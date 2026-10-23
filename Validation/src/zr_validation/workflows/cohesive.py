#!/usr/bin/env python3

import json
from zr_validation.settings import *
from zr_validation.core.common import pair_block,write_file_if_missing
from zr_validation.core.structures import crystal,write_lmp,write_isolated_atom


def cohesive_input(pot,which,bulk_name=''):
    is_bulk=(which=='bulk')
    structure=f'../configuration/{bulk_name}' if is_bulk else '../configuration/single_atom.data'
    relax='''min_style cg
min_modify line quadratic
minimize 1.0e-12 1.0e-12 100000 100000
''' if is_bulk else ''
    return f'''units metal
atom_style atomic
boundary p p p
read_data {structure}
{pair_block(str(pot))}
mass 1 91.224
neighbor 1.0 nsq
neigh_modify delay 0 every 1 check yes
thermo 10
thermo_style custom step pe etotal press vol
run 0
{relax}run 0
variable E equal pe
variable N equal count(all)
print "COH_RESULT E=${{E}} N=${{N}}" file result.txt screen yes
'''


def cohesive_execute_sh():
    return f'''#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
LMP="{LMP_EXE}"
pids=()
for d in bulk single; do
  if [[ -s "$d/result.txt" && -s "$d/log.lammps" ]] && grep -q '^COH_RESULT E=' "$d/result.txt"; then
    echo "[skip] $d"
    continue
  fi
  ( cd "$d" && srun --exclusive --ntasks=1 --cpus-per-task=1 "$LMP" -in in.lmp > stdout.txt 2>&1 ) &
  pids+=("$!")
done
failed=0
for pid in "${{pids[@]}}"; do wait "$pid" || failed=1; done
if (( failed )); then echo 'ERROR: cohesive LAMMPS calculation failed' >&2; exit 1; fi
python3 collect_energy.py
'''


def cohesive_collect_py():
    return '''#!/usr/bin/env python3
import csv,json,re,sys
from pathlib import Path
root=Path(__file__).resolve().parent
config=json.loads((root/'simulation_config.json').read_text())
def get(which):
    fp=root/which/'result.txt'
    if not fp.is_file() or not (root/which/'log.lammps').is_file():
        sys.exit(f'Missing {which} result/log')
    m=re.search(r'COH_RESULT E=([^\\s]+) N=([^\\s]+)',fp.read_text())
    if not m: sys.exit(f'Malformed {which} result')
    return float(m.group(1)),int(float(m.group(2)))
Eb,N=get('bulk')
Ea,Na=get('single')
if Na!=1 or N!=config['natoms']:
    sys.exit('Wrong number of atoms')
binding=Ea-Eb/N
row=dict(phase=config['phase'],natoms=N,E_bulk_total_eV=Eb,
         E_bulk_per_atom_eV=Eb/N,E_atom_eV=Ea,cohesive_energy_positive_eV_atom=binding,
         bulk_minus_atom_eV_atom=-binding)
(root/'cohesive_result.json').write_text(json.dumps(row,indent=2)+'\\n')
with (root/'cohesive_result.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(row))
    w.writeheader(); w.writerow(row)
print(f'[OK] {config["phase"]} binding cohesive energy: {binding:.8f} eV/atom')
'''


def cohesive_run_sh(test_dir):
    return f'''#!/bin/bash
#SBATCH --account={SLURM_ACCOUNT}
#SBATCH --nodes=1
#SBATCH --ntasks={COHESIVE_CORES}
#SBATCH --cpus-per-task=1
#SBATCH --time={BATCH_TIME}
module purge
module load StdEnv/2023
module load gcc/12.3
module load openmpi/4.1.5
cd "{test_dir}"
bash execute_cohesive.sh
'''


def setup_cohesive(test_name,pot,base_dir):
    phase=test_name.split('-',1)[1]
    d=base_dir/test_name
    d.mkdir(parents=True,exist_ok=True)
    cfg=COHESIVE_CONFIG[phase]
    cell,pos=crystal(phase,cfg['a'],c=cfg.get('c'),repeat=cfg['repeat'],conventional=(phase=='BCC'))
    cell_prim,pos_prim=crystal(phase,cfg['a'],c=cfg.get('c'),conventional=(phase=='BCC'))
    bulk_name=cfg['filename']
    meta={'phase':phase,'a':cfg['a'],'c':cfg.get('c'), 'repeat':cfg['repeat'],
          'bulk_data':bulk_name,'natoms':len(pos),'potential':str(pot),'nlh':TABLE_POTENTIAL}
    path=d/'simulation_config.json'
    if path.exists() and json.loads(path.read_text())!=meta:
        raise RuntimeError(f'{d}: cohesive configuration changed; use a new directory')
    write_file_if_missing(path,json.dumps(meta,indent=2,ensure_ascii=False)+'\n')
    structure=d/'configuration'
    if not (structure/bulk_name).exists(): write_lmp(structure/bulk_name,cell,pos)
    if not (structure/cfg['primitive']).exists(): write_lmp(structure/cfg['primitive'],cell_prim,pos_prim)
    if phase=='FCC':

        extra_cell,extra_pos=crystal('FCC',cfg['a'],repeat=(3,3,3))
        if not (structure/'fcc_3x3x3.data').exists():
            write_lmp(structure/'fcc_3x3x3.data',extra_cell,extra_pos)
    if not (structure/'single_atom.data').exists(): write_isolated_atom(structure/'single_atom.data')
    for which in ('bulk','single'):
        write_file_if_missing(d/which/'in.lmp',cohesive_input(pot,which,bulk_name))
    write_file_if_missing(d/'execute_cohesive.sh',cohesive_execute_sh(),executable=True)
    write_file_if_missing(d/'run.sh',cohesive_run_sh(d),executable=True)
    write_file_if_missing(d/'collect_energy.py',cohesive_collect_py(),executable=True)
    write_file_if_missing(d/'collect_energy.sh',
                          '#!/bin/bash\nset -euo pipefail\ncd "$(dirname "$0")"\npython3 collect_energy.py\n',executable=True)
