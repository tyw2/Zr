#!/usr/bin/env python3

import json
from decimal import Decimal
from pathlib import Path
from zr_validation.settings import *
from zr_validation.core.common import pair_block, write_file, write_file_if_missing
from zr_validation.core.structures import crystal,write_lmp,append_cfg


def eos_values():
    start,stop,step=map(lambda x:Decimal(str(x)),(EOS_A_START,EOS_A_END,EOS_A_STEP))
    if step<=0 or stop<start: raise ValueError('Invalid EOS a range')
    v=start
    while v<=stop:
        yield v
        v+=step


def eos_case_name(phase, a):
    name=f'a_{a:.6f}'
    if phase=='HCP':
        c=a*Decimal(str(EOS_HCP_CA))
        name+=f'_c_{c:.6f}'
    return name


def eos_input(phase,pot):
    relax='''min_style cg
min_modify line quadratic
minimize 1.0e-12 1.0e-12 2000 20000
''' if EOS_RELAX_ATOMS[phase] else ''
    return f'''units metal
atom_style atomic
boundary p p p
read_data data.lmp
{pair_block(str(pot))}
mass 1 91.224
neighbor 1.0 nsq
neigh_modify delay 0 every 1 check yes
thermo 1
thermo_style custom step pe ke etotal press vol lx ly lz
run 0
{relax}run 0
variable E equal pe
variable V equal vol
print "EOS_RESULT E=${{E}} V=${{V}}" file result.txt screen yes
'''


def eos_execute_sh(cores):
    return f'''#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
LMP="{LMP_EXE}"
pids=()
failed=0
for dir in structures/a_*/; do
  [[ -d "$dir" ]] || continue
  if [[ -s "$dir/result.txt" && -s "$dir/log.lammps" ]] && grep -q '^EOS_RESULT E=' "$dir/result.txt"; then
    echo "[skip] $dir"
    continue
  fi
  ( cd "$dir" && srun --exclusive --ntasks=1 --cpus-per-task=1 "$LMP" -in in.lmp > stdout.txt 2>&1 ) &
  pids+=("$!")
  if (( ${{#pids[@]}} >= {cores} )); then
    for pid in "${{pids[@]}}"; do wait "$pid" || failed=1; done
    pids=()
  fi
done
for pid in "${{pids[@]}}"; do wait "$pid" || failed=1; done
if (( failed )); then echo 'ERROR: Some EOS configurations failed' >&2; exit 1; fi
python3 collect_energy.py
'''


def eos_collect_py():
    return '''#!/usr/bin/env python3
import csv,json,re,sys
from pathlib import Path
root=Path(__file__).resolve().parent
rows=[]
for d in sorted(root.glob('structures/a_*')):
    if not (d/'data.lmp').exists(): continue
    meta=json.loads((d/'meta.json').read_text())
    fp=d/'result.txt'
    if not fp.is_file() or not (d/'log.lammps').is_file():
        sys.exit(f'Missing result/log at {d}')
    m=re.search(r'EOS_RESULT E=([^\\s]+) V=([^\\s]+)',fp.read_text())
    if not m: sys.exit(f'Invalid EOS result at {d}')
    E,V=map(float,m.groups())
    n=meta['natoms']
    rows.append(dict(a_A=meta['a'],c_A=meta.get('c',''),natoms=n,volume_A3=V,
                     volume_per_atom_A3=V/n,E_total_eV=E,E_per_atom_eV=E/n))
if not rows or len(rows)!=len(list(root.glob('structures/a_*/meta.json'))):
    sys.exit('Incomplete EOS scan')
rows.sort(key=lambda r:r['a_A'])
with (root/'eos_result.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]))
    writer.writeheader(); writer.writerows(rows)
print(f'[OK] {len(rows)} EOS points -> {root / "eos_result.csv"}')
'''


def eos_plot_py():
    return '''#!/usr/bin/env python3
"""Run inside EOS-BCC or EOS-HCP after collecting data (requires matplotlib)."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path(__file__).resolve().parent
with (p/'eos_result.csv').open() as f: rows=list(csv.DictReader(f))
v=[float(r['volume_per_atom_A3']) for r in rows]
e=[float(r['E_per_atom_eV']) for r in rows]
plt.figure(figsize=(7,5))
plt.plot(v,e,'o-',markersize=3)
plt.xlabel('Volume / atom (Angstrom^3)')
plt.ylabel('Energy / atom (eV)')
plt.title(p.name + ' Zr EOS')
plt.grid(alpha=.3)
plt.tight_layout()
plt.savefig(p/'eos_curve.png',dpi=180)
print('Saved',p/'eos_curve.png')
'''

def eos_run_sh(test_dir,cores):
    return f'''#!/bin/bash
#SBATCH --account={SLURM_ACCOUNT}
#SBATCH --nodes=1
#SBATCH --ntasks={cores}
#SBATCH --cpus-per-task=1
#SBATCH --time={BATCH_TIME}
module purge
module load StdEnv/2023
module load gcc/12.3
module load openmpi/4.1.5
cd "{test_dir}"
bash execute_eos.sh
'''


def setup_eos(test_name,pot,base_dir):
    phase=test_name.split('-',1)[1]
    test_dir=base_dir/test_name
    test_dir.mkdir(parents=True,exist_ok=True)
    config={'potential':str(pot),'phase':phase,'a_range':[str(EOS_A_START),str(EOS_A_END),str(EOS_A_STEP)],
            'ca':str(EOS_HCP_CA) if phase=='HCP' else None,'repeat':list(EOS_REPEAT),
            'relax_atoms':EOS_RELAX_ATOMS[phase], 'nlh':TABLE_POTENTIAL}
    fp=test_dir/'simulation_config.json'
    if fp.exists() and json.loads(fp.read_text())!=config:
        raise RuntimeError(f'{test_dir}: EOS settings/potential changed; use another OUTPUT_ROOT or archive directory')
    write_file_if_missing(fp,json.dumps(config,indent=2,ensure_ascii=False)+'\n')
    summary=test_dir/'structures'/f'all_{phase.lower()}.cfg'
    if not summary.exists():
        summary.parent.mkdir(parents=True,exist_ok=True)
        summary.write_text('')

    all_complete=(summary.exists() and summary.stat().st_size>0 and
                  all((test_dir/'structures'/eos_case_name(phase,a)/'meta.json').exists() for a in eos_values()))
    if not all_complete:
        summary.write_text('')
    for a in eos_values():
        aval=float(a)
        cval=float(a*Decimal(str(EOS_HCP_CA))) if phase=='HCP' else None
        d=test_dir/'structures'/eos_case_name(phase,a)
        cell,pos=crystal(phase,aval,c=cval,repeat=EOS_REPEAT)
        if not all_complete or not (d/'data.lmp').exists():
            write_lmp(d/'data.lmp',cell,pos)
        meta={'element':'Zr','phase':phase,'a':aval,'c':cval,'natoms':len(pos),'repeat':EOS_REPEAT}
        write_file_if_missing(d/'meta.json',json.dumps(meta,indent=2)+'\n')
        write_file_if_missing(d/'in.lmp',eos_input(phase,pot))
        if not all_complete:
            append_cfg(summary,cell,pos, f'a={aval:.6f} structure={phase.lower()}'+
                       (f' c={cval:.6f} c_over_a={EOS_HCP_CA}' if cval else ''))
    write_file_if_missing(test_dir/'execute_eos.sh',eos_execute_sh(EOS_CORES),executable=True)
    write_file_if_missing(test_dir/'run.sh',eos_run_sh(test_dir,EOS_CORES),executable=True)
    write_file_if_missing(test_dir/'collect_energy.py',eos_collect_py(),executable=True)
    write_file_if_missing(test_dir/'plot_eos.py',eos_plot_py(),executable=True)
    write_file_if_missing(test_dir/'collect_energy.sh',
                          '#!/bin/bash\nset -euo pipefail\ncd "$(dirname "$0")"\npython3 collect_energy.py\n',executable=True)
