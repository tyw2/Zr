#!/usr/bin/env python3

import json
import shutil
from pathlib import Path
from zr_validation.settings import (FORMATION_CONFIG, FORMATION_CORES, LMP_EXE,
                      SLURM_ACCOUNT, BATCH_TIME, TABLE_POTENTIAL)
from zr_validation.core.common import pair_block,write_file_if_missing
from zr_validation.workflows.formation.structures import generate_structures


def lmp_input(pot, conf, label):

    relax = (f'fix cell all box/relax aniso 0.0 vmax 0.001\n'
             if conf['relaxation'] == 'aniso' else '')
    if conf['relaxation'] not in ('aniso', 'fixed'):
        raise ValueError('Formation relaxation must be aniso or fixed')
    unfix = 'unfix cell\n' if relax else ''
    return f'''units metal
atom_style atomic
boundary p p p
read_data data.lmp
mass 1 91.224
{pair_block(str(pot))}
neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes
thermo 50
thermo_style custom step pe etotal press pxx pyy pzz lx ly lz fmax fnorm
min_style cg
min_modify line quadratic dmax 0.1
{relax}minimize 1.0e-12 1.0e-12 {conf['maxiter']} {conf['maxeval']}
{unfix}run 0
variable E equal pe
variable N equal count(all)
variable FM equal fmax
print "FORM_RESULT E=${{E}} N=${{N}} FMAX=${{FM}}" file result.txt screen yes
write_data relaxed.data
write_dump all custom relaxed.dump id type x y z fx fy fz modify sort id
'''


def execute_script(subs):
    sub_block = ' '.join(subs)
    return f'''#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
LMP="{LMP_EXE}"
pids=()
for d in {sub_block}; do
  if [[ -s "$d/result.txt" && -s "$d/log.lammps" ]] && \\
     grep -q '^FORM_RESULT E=' "$d/result.txt" && \\
     grep -q 'Total wall time:' "$d/log.lammps"; then
    echo "[skip] $d"
    continue
  fi
  (
    cd "$d"
    rm -f result.txt
    srun --exclusive --ntasks=1 --cpus-per-task=1 \\
       "$LMP" -in in.lmp -log log.lammps > stdout.txt 2>&1
  ) &
  pids+=("$!")
done
failed=0
for pid in "${{pids[@]}}"; do
  wait "$pid" || failed=1
done
if (( failed )); then
  echo 'ERROR: at least one formation-energy LAMMPS job failed' >&2
  exit 1
fi
python3 collect_energy.py
'''


def run_script(d, cores):
    return f'''#!/bin/bash
#SBATCH --account={SLURM_ACCOUNT}
#SBATCH --nodes=1
#SBATCH --ntasks={cores}
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time={BATCH_TIME}
module purge
module load StdEnv/2023
module load gcc/12.3
module load openmpi/4.1.5
export OMP_NUM_THREADS=1
cd "{d}"
bash execute_formation.sh
'''


def setup_formation(test_name, potential, base_dir):
    conf = dict(FORMATION_CONFIG[test_name])
    conf['test'] = test_name
    conf['repeat'] = list(conf['repeat'])
    conf['potential'] = str(potential)
    conf['nlh_table'] = TABLE_POTENTIAL
    d = base_dir / test_name
    d.mkdir(parents=True,exist_ok=True)
    simulation_file = d / 'simulation_config.json'
    if simulation_file.exists() and json.loads(simulation_file.read_text()) != conf:
        raise RuntimeError(f'{d}: formation settings changed; use new OUTPUT_ROOT or archive old folder')
    write_file_if_missing(simulation_file, json.dumps(conf,indent=2,ensure_ascii=False)+'\n')
    details, ncases = generate_structures(d, conf) if not (d/'structure_info.json').exists() else (json.loads((d/'structure_info.json').read_text()), 3 if conf['kind'] in ('basal','out_of_plane') else 2)

    subdirs = ['bulk', 'defect'] + (['monovac'] if ncases == 3 else [])
    if any(not (d/sub/'data.lmp').is_file() for sub in subdirs):
        details, ncases = generate_structures(d, conf)
    write_file_if_missing(d/'structure_info.json', json.dumps(details,indent=2,ensure_ascii=False)+'\n')
    for label in subdirs:
        write_file_if_missing(d/label/'in.lmp', lmp_input(potential,conf,label))
    local = d/'collect_energy.py'
    if not local.is_file():
        shutil.copy2(Path(__file__).with_name('post.py'), local)
    write_file_if_missing(d/'execute_formation.sh', execute_script(subdirs), executable=True)
    write_file_if_missing(d/'run.sh',run_script(d, FORMATION_CORES[test_name]),executable=True)
    write_file_if_missing(d/'collect_energy.sh',
           '#!/bin/bash\nset -euo pipefail\ncd "$(dirname "$0")"\npython3 collect_energy.py\n',executable=True)
    print(f'  [formation] {d}: N={details["N_bulk"]}, defect N={details["N_defect"]}, {ncases} subruns')
