#!/usr/bin/env python3


from zr_validation.settings import *
from zr_validation.core.common import pair_block

def neb_in(pot_path: str, start_idx: int = 0, title: str = "defect migration") -> str:


    ivar = "part" if start_idx == 0 else f"part+{start_idx}"

    return f"""
units           metal
atom_style      atomic
boundary        p p p
atom_modify     map array

read_data       replica/replica.{start_idx}.lmp

{pair_block(pot_path)}

mass            1 91.224

min_style       fire

variable        i equal {ivar}

thermo          50
thermo_style    custom step pe press fmax fnorm

dump            D all custom 100 dump.neb.${{i}}.* id type x y z fx fy fz
dump_modify     D sort id

fix             1 all neb 1.0

neb             0.0 0.001 4000 4000 20 each replica/replica.${{i}}.coords verbosity verbose
"""

def relax_endpoint_in(pot_path: str, idx: int, which: str) -> str:

    return f"""
units           metal
atom_style      atomic
boundary        p p p

read_data       replica/replica.{idx}.lmp

{pair_block(pot_path)}

mass            1 91.224
thermo          50
thermo_style    custom step pe press fmax fnorm

min_style       fire
minimize        1.0e-10 1.0e-10 10000 10000

write_data      replica/replica.{idx}.lmp
"""

_LMP_TO_COORDS_FUNC = r'''
def lmp_to_coords(lmp_file, coords_file):
    """把 lmp_file 第 LMP_START_LINE–LMP_END_LINE 行写入 coords_file,
    然后按原子 ID 排序重写(格式: id x y z)。"""

    with open(lmp_file, "r") as f:
        lmp_lines = f.readlines()

    replacement_lines = lmp_lines[LMP_START_LINE - 1:LMP_END_LINE]

    if len(replacement_lines) != N_ATOMS:
        raise ValueError(
            "{} 中第 {}-{} 行数量不正确:预期 {} 行,实际得到 {} 行".format(
                lmp_file, LMP_START_LINE, LMP_END_LINE,
                N_ATOMS, len(replacement_lines)))

    with open(coords_file, "w") as f:
        f.write("{}\n".format(N_ATOMS))
        f.writelines(replacement_lines)

    print("[OK] wrote {} atoms from lines {}-{} of {} into {}".format(
        N_ATOMS, LMP_START_LINE, LMP_END_LINE, lmp_file, coords_file))

    with open(coords_file, "r") as f:
        lines = [line.split() for line in f if line.strip()]

    if not lines:
        raise ValueError("{} 是空文件".format(coords_file))

    N = int(lines[0][0])
    atoms = []

    for line_number, parts in enumerate(lines[1:], start=2):
        if len(parts) < 5:
            raise ValueError(
                "{} 第 {} 行格式不正确:至少需要 5 列,实际只有 {} 列".format(
                    coords_file, line_number, len(parts)))

        atom_id = int(parts[0])
        x = float(parts[2])
        y = float(parts[3])
        z = float(parts[4])
        atoms.append((atom_id, x, y, z))

    if len(atoms) != N:
        raise ValueError(
            "文件第一行声明有 {} 个原子,但实际读取到 {} 个原子".format(
                N, len(atoms)))

    atom_ids = [atom[0] for atom in atoms]
    if len(set(atom_ids)) != len(atom_ids):
        duplicate_ids = sorted(
            atom_id for atom_id in set(atom_ids)
            if atom_ids.count(atom_id) > 1)
        raise ValueError("发现重复的原子 ID:{}".format(duplicate_ids))

    atoms.sort(key=lambda atom: atom[0])

    with open(coords_file, "w") as f:
        f.write("{}\n".format(N))
        for atom_id, x, y, z in atoms:
            f.write("{} {:.16f} {:.16f} {:.16f}\n".format(atom_id, x, y, z))

    print("[OK] sorted {} atoms by atom ID and wrote {}".format(N, coords_file))
'''

def endpoint_helper_py(nrep: int = NEB_NREP,
                       start_line: int = 17,
                       end_line: int = 161,
                       start_idx: int = 0) -> str:


    first_idx = start_idx
    final_idx = start_idx + nrep - 1
    header = f'''#!/usr/bin/env python3
"""Update NEB endpoint coords from relaxed LAMMPS data files.

  replica/replica.{first_idx}.lmp -> replica/replica.{first_idx}.coords
  replica/replica.{final_idx}.lmp -> replica/replica.{final_idx}.coords

从 .lmp 第 {start_line}-{end_line} 行截取 Atoms 数据({end_line - start_line + 1} 原子),
写入 coords 并按原子 ID 排序(格式: id x y z)。
"""
from pathlib import Path
import sys

LMP_START_LINE = {start_line}
LMP_END_LINE = {end_line}
N_ATOMS = LMP_END_LINE - LMP_START_LINE + 1

PAIRS = [
    ("replica/replica.{first_idx}.lmp", "replica/replica.{first_idx}.coords"),
    ("replica/replica.{final_idx}.lmp", "replica/replica.{final_idx}.coords"),
]

'''
    footer = '''

def update():
    for lmp_file, coords_file in PAIRS:
        if not Path(lmp_file).exists():
            raise FileNotFoundError(lmp_file)
        lmp_to_coords(lmp_file, coords_file)


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] != "update":
        raise SystemExit("Usage: python3 endpoint_helper.py update")
    update()
'''
    return header + _LMP_TO_COORDS_FUNC + footer

def single_neb_test_run_sh(np: int, input_file: str) -> str:


    return f"""#!/bin/bash
#SBATCH --account={SLURM_ACCOUNT}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={np}
#SBATCH --cpus-per-task=1
#SBATCH --time=0-00:59

module purge
module load StdEnv/2023
module load gcc/12.3
module load openmpi/4.1.5

set -e
mpirun -np 1 {LMP_EXE} -in in_relax_initial.lmp > relax_initial.out 2>&1
mpirun -np 1 {LMP_EXE} -in in_relax_final.lmp   > relax_final.out 2>&1
python3 endpoint_helper.py update
mpirun -np {np} {LMP_EXE} -partition {np}x1 -in {input_file}
"""

def neb_collect_sh() -> str:
    return r"""#!/bin/bash
output_file="./neb_energy_output.txt"
> "$output_file"

shopt -s nullglob
logs=( log.lammps.* )
if [ ${#logs[@]} -eq 0 ]; then
    logs=( log.lammps )
fi

for logf in "${logs[@]}"; do
    if [ -f "$logf" ]; then
        last_step=$(grep -E '^[[:space:]]*[0-9]+[[:space:]]' "$logf" | tail -n 1)
        if [ -n "$last_step" ]; then
            echo "$logf : $last_step" >> "$output_file"
        else
            echo "$logf : no step line found" >> "$output_file"
        fi
    fi
done
echo "Done. Saved to $output_file"
"""

def detect_start_index(src: Path) -> int:


    if not src.exists():
        print(f"    [warn] {src} not found; assuming start index 0.")
        return 0
    for start in (0, 1):
        if (src / f"replica.{start}.coords").exists() or \
           (src / f"replica.{start}.lmp").exists():
            missing = [i for i in range(start, start + NEB_NREP)
                       if not (src / f"replica.{i}.coords").exists()]
            if missing:
                print(f"    [warn] {src}: start={start}, "
                      f"but missing replica.<i>.coords for images {missing}")
            return start
    print(f"    [warn] cannot detect image start index in {src}; defaulting to 0.")
    return 0

