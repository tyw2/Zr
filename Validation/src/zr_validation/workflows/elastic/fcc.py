#!/usr/bin/env python3


import json
import shutil
from pathlib import Path
from zr_validation.settings import (LMP_EXE, SLURM_ACCOUNT, BATCH_TIME, ELASTIC_FCC_A,
                      ELASTIC_FCC_REPEAT, ELASTIC_FCC_STRAIN,
                      ELASTIC_FCC_BOX_RELAX, ELASTIC_FCC_ETOL,
                      ELASTIC_FCC_FTOL, ELASTIC_FCC_MAXITER,
                      ELASTIC_FCC_MAXEVAL, ELASTIC_FCC_DMAX,
                      TABLE_POTENTIAL, TABLE_NAME)
from zr_validation.core.common import pair_block, write_file_if_missing
from zr_validation.core.structures import crystal, write_lmp, append_cfg


DIRECTIONS = {1: 'x', 2: 'y', 3: 'z', 4: 'yz', 5: 'xz', 6: 'xy'}


def lammps_prefix(potential: str, reference_data: str) -> str:
    return f'''units metal
atom_style atomic
boundary p p p
read_data {reference_data}
mass 1 91.224
{pair_block(potential)}
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
thermo 50
thermo_style custom step pe pxx pyy pzz pyz pxz pxy press vol fmax fnorm lx ly lz
min_style cg
min_modify line quadratic dmax {ELASTIC_FCC_DMAX}
'''


def minimize_command() -> str:
    return (f'minimize {ELASTIC_FCC_ETOL:.12g} {ELASTIC_FCC_FTOL:.12g} '
            f'{ELASTIC_FCC_MAXITER} {ELASTIC_FCC_MAXEVAL}')


def stress_print(path: str, marker: str):
    return f'''run 0
variable E equal pe
variable V equal vol
variable N equal count(all)
variable F equal fmax
variable P1 equal pxx
variable P2 equal pyy
variable P3 equal pzz
variable P4 equal pyz
variable P5 equal pxz
variable P6 equal pxy
print "{marker} E=${{E}} V=${{V}} N=${{N}} FMAX=${{F}} PXX=${{P1}} PYY=${{P2}} PZZ=${{P3}} PYZ=${{P4}} PXZ=${{P5}} PXY=${{P6}}" file {path} screen yes
'''


def equil_input(potential: str) -> str:

    return (lammps_prefix(potential, 'init.lammps')
            + ''
            + minimize_command() + '\n'
            + ''
            + f'fix boxeq all box/relax {ELASTIC_FCC_BOX_RELAX} 0.0 vmax 0.001\n'
            + minimize_command() + '\n'
            + 'unfix boxeq\n'
            + 'write_data equilibrated.data\n'
            + 'write_dump all custom relaxed.dump id type xu yu zu fx fy fz modify sort id\n'
            + stress_print('reference.txt', 'ELASTIC_REF'))


def strain_input(potential: str, direction: int, sign: int) -> str:
    eps = ELASTIC_FCC_STRAIN * sign
    field = DIRECTIONS[direction]
    if direction <= 3:
        change = f'change_box all {field} scale {1.0 + eps:.12f} remap units box'
    else:

        length = 'lz' if direction in (4, 5) else 'ly'
        change = (f'variable dtilt equal ({eps:.12g})*{length}\n'
                  f'change_box all {field} delta ${{dtilt}} remap units box')
    name = f'strain_{direction}_{"pos" if sign > 0 else "neg"}'
    return (lammps_prefix(potential, 'equilibrated.data')
            + 'change_box all triclinic\n'
            + change + '\n'
            + ''
            + minimize_command() + '\n'
            + stress_print(f'{name}.txt', 'ELASTIC_STRESS'))


def execute_script():
    pairs = '\n'.join(
        f'run_one in.strain_{j}_{s}.lmp strain_{j}_{s}.txt'
        for j in range(1, 7) for s in ('neg', 'pos'))
    return f'''#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
LMP="{LMP_EXE}"
if [[ -s elastic_result.json && -s elastic_result.csv ]]; then
    echo '[skip] Completed elastic tensor'
    exit 0
fi
run_one() {{
    local input="$1" result="$2"
    local tag="${{input#in.}}"
    rm -f "$result" "log.${{tag}}"
    srun --exclusive --ntasks=1 --cpus-per-task=1 \\
         "$LMP" -in "$input" -log "log.${{tag}}" > "stdout.${{tag}}" 2>&1
    [[ -s "$result" ]] || {{ echo "ERROR: missing $result" >&2; return 1; }}
}}
if [[ ! -s reference.txt || ! -s equilibrated.data ]]; then
    rm -f reference.txt equilibrated.data
    run_one in.equilibrate.lmp reference.txt
fi
{pairs}
python3 collect_energy.py
'''


def single_run_script(test_dir: Path):
    return f'''#!/bin/bash
#SBATCH --account={SLURM_ACCOUNT}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time={BATCH_TIME}
module purge
module load StdEnv/2023
module load gcc/12.3
module load openmpi/4.1.5
export OMP_NUM_THREADS=1
cd "{test_dir}"
bash execute_elastic.sh
'''


def setup_elastic_fcc(test_name, pot, base_dir):
    dest = base_dir / test_name
    dest.mkdir(parents=True, exist_ok=True)
    parameters = dict(phase='FCC', potential=str(pot), nlh=TABLE_POTENTIAL,
                      a=ELASTIC_FCC_A, repeat=list(ELASTIC_FCC_REPEAT),
                      strain=ELASTIC_FCC_STRAIN,
                      box_relax=ELASTIC_FCC_BOX_RELAX,
                      etol=ELASTIC_FCC_ETOL, ftol=ELASTIC_FCC_FTOL,
                      maxiter=ELASTIC_FCC_MAXITER,
                      maxeval=ELASTIC_FCC_MAXEVAL, dmax=ELASTIC_FCC_DMAX)
    cfg_path = dest / 'simulation_config.json'
    if cfg_path.exists() and json.loads(cfg_path.read_text()) != parameters:
        raise RuntimeError(f'{dest} settings changed: choose a new OUTPUT_ROOT or archive this phase folder')
    write_file_if_missing(cfg_path, json.dumps(parameters, indent=2, ensure_ascii=False)+'\n')
    cell, positions = crystal('FCC', ELASTIC_FCC_A, conventional=True,
                              repeat=ELASTIC_FCC_REPEAT)
    if not (dest / 'init.lammps').exists():
        write_lmp(dest / 'init.lammps', cell, positions)
    if not (dest / 'init.cfg').exists():
        append_cfg(dest / 'init.cfg', cell, positions, 'initial_FCC_elastic')
    write_file_if_missing(dest/'in.equilibrate.lmp', equil_input(str(pot)))
    for j in range(1, 7):
        for s in (-1, 1):
            name = f'in.strain_{j}_{"pos" if s > 0 else "neg"}.lmp'
            write_file_if_missing(dest / name, strain_input(str(pot), j, s))
    write_file_if_missing(dest/'execute_elastic.sh', execute_script(), executable=True)
    write_file_if_missing(dest/'run.sh', single_run_script(dest), executable=True)
    local_collector = dest/'collect_energy.py'
    if not local_collector.exists():
        shutil.copy2(Path(__file__).with_name('post.py'), local_collector)
    write_file_if_missing(dest/'collect_energy.sh',
                          '#!/bin/bash\nset -euo pipefail\ncd "$(dirname "$0")"\npython3 collect_energy.py\n',
                          executable=True)
    print(f'  [elastic] {dest}: {len(positions)} atoms, ±{ELASTIC_FCC_STRAIN:g} strain')
