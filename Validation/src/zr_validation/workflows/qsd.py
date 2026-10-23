#!/usr/bin/env python3


from zr_validation.settings import *
from zr_validation.core.common import pair_block

def qsd_steps(test_name: str) -> int:

    return round(QSD_CONFIG[test_name]["dmax"] / QSD_D_STEP)

def qsd_lattice_block() -> str:
    if QSD_LATTICE_MODE == "legacy_hcp":

        return f"lattice         hcp {QSD_A}"
    if QSD_LATTICE_MODE == "custom_ca":

        return f"""lattice         custom {QSD_A} &
                a1 1.0 0.0 0.0 &
                a2 0.0 1.7320508075688772 0.0 &
                a3 0.0 0.0 {QSD_CA} &
                basis 0.0 0.0 0.0 &
                basis 0.5 0.5 0.0 &
                basis 0.5 0.8333333333333333 0.5 &
                basis 0.0 0.3333333333333333 0.5"""
    raise ValueError(f"Bad QSD_LATTICE_MODE: {QSD_LATTICE_MODE}")

def qsd_input(test_name: str, pot_path: str) -> str:

    cfg = QSD_CONFIG[test_name]
    drag = cfg["drag"]
    header = f"""
units           metal
dimension       3
boundary        p p p
atom_style      atomic
atom_modify     map array

{qsd_lattice_block()}
region          box block 0 4 0 2 0 2 units lattice
create_box      1 box
create_atoms    1 box
mass            1 91.224

{pair_block(pot_path, table_points=4000 if test_name == "QSD-AC" else 1000)}
neighbor        2.0 bin
neigh_modify    every 1 delay 0 check yes
thermo          1
thermo_style    custom step atoms pe

variable        id_drag equal {drag}
group           gdrag id ${{id_drag}}
"""
    if test_name == "QSD-A":
        movement = """
variable ex equal -1.0
variable ey equal  0.0
variable ez equal  0.0
"""
    elif test_name == "QSD-C":
        movement = """
variable ex equal  0.0
variable ey equal  0.0
variable ez equal -1.0
"""
    elif test_name == "QSD-AC":
        target = cfg["target"]
        movement = f"""
variable id_target equal {target}
variable x1 equal x[${{id_drag}}]
variable y1 equal y[${{id_drag}}]
variable z1 equal z[${{id_drag}}]
variable x2 equal x[${{id_target}}]
variable y2 equal y[${{id_target}}]
variable z2 equal z[${{id_target}}]
variable dx0 equal v_x2-v_x1
variable dy0 equal v_y2-v_y1
variable dz0 equal v_z2-v_z1
variable r0 equal sqrt(v_dx0*v_dx0+v_dy0*v_dy0+v_dz0*v_dz0)
variable ex equal v_dx0/v_r0
variable ey equal v_dy0/v_r0
variable ez equal v_dz0/v_r0
print \"QSD-AC: drag=${{id_drag}} target=${{id_target}} initial_AC_distance_A=${{r0}}\"
print \"QSD-AC: direction=${{ex}} ${{ey}} ${{ez}}\"
"""
    else:
        raise ValueError(test_name)

    tail = """
variable        dx equal ${d}*${ex}
variable        dy equal ${d}*${ey}
variable        dz equal ${d}*${ez}
displace_atoms  gdrag move ${dx} ${dy} ${dz} units box
fix             frozen gdrag setforce 0.0 0.0 0.0

"""
    if QSD_RELAX_OTHERS:
        tail += """
min_style       cg
min_modify      line quadratic
minimize        0.0 1.0e-8 10000 100000
"""
    else:
        tail += """min_style       cg
min_modify      line quadratic
run             0
"""
    tail += """
variable        E equal pe
print \"d=${d} E=${E}\" file energy.txt screen yes
"""
    return header + movement + tail

def qsd_execute_sh(test_name: str) -> str:

    nsteps = qsd_steps(test_name)
    return f'''#!/bin/bash
set -euo pipefail
LMP_BIN="${{LMP:-{LMP_EXE}}}"
MAX_PARALLEL={QSD_PARALLEL}
LAST_STEP={nsteps}
mkdir -p .
failed=0
pids=()

wait_wave() {{
    local pid
    for pid in "${{pids[@]}}"; do
        if ! wait "$pid"; then
            failed=1
        fi
    done
    pids=()
    if (( failed != 0 )); then
        echo "ERROR: one or more {test_name} displacement jobs failed" >&2
        exit 1
    fi
}}

for (( i=0; i<=LAST_STEP; i++ )); do
    d="$((i / 10)).$((i % 10))"
    step_dir="Step_${{i}}"
    mkdir -p "$step_dir"
    if [[ -s "$step_dir/energy.txt" ]] && \\
       grep -Eq '^d=[0-9]+[.][0-9]+ E=[+-]?[0-9.eE+-]+' "$step_dir/energy.txt"; then
        continue
    fi
    cp qsd.in "$step_dir/qsd.in"
    (
        cd "$step_dir"
        srun --exclusive --ntasks=1 --cpus-per-task=1 \\
             "$LMP_BIN" -in qsd.in -var d "$d" > stdout.txt 2>&1
    ) &
    pids+=("$!")
    if (( ${{#pids[@]}} >= MAX_PARALLEL )); then
        wait_wave
    fi
done
wait_wave

echo "All {test_name} displacement points are ready."
'''

def qsd_single_run_sh(test_name: str) -> str:
    return f'''#!/bin/bash
#SBATCH --account={SLURM_ACCOUNT}
#SBATCH --job-name={test_name.replace('-', '_')}
#SBATCH --nodes=1
#SBATCH --ntasks={min(QSD_PARALLEL, qsd_steps(test_name) + 1)}
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time={BATCH_TIME}

module purge
module load StdEnv/2023
module load gcc/12.3
module load openmpi/4.1.5

set -euo pipefail
cd "$(dirname "$0")"
export LMP={LMP_EXE}
bash execute_qsd.sh
bash collect_energy.sh
'''

def qsd_collect_py(test_name: str) -> str:
    cfg = QSD_CONFIG[test_name]
    nsteps = qsd_steps(test_name)
    legacy = cfg["legacy_txt"]

    return r'''#!/usr/bin/env python3
"""Collect QSD data: d, E(d), and delta E = E(d)-E(0) (eV)."""
import csv
import math
import re
from pathlib import Path

LAST_STEP = ''' + str(nsteps) + r'''
LEGACY_FILE = ''' + repr(legacy) + r'''
D_STEP = ''' + repr(QSD_D_STEP) + r'''
RX = re.compile(r"^d=([^ ]+)\s+E=([^ ]+)\s*$")
rows = []
for i in range(LAST_STEP + 1):
    p = Path(f"Step_{i}/energy.txt")
    if not p.is_file():
        raise SystemExit(f"ERROR: missing {p}")
    text = p.read_text().strip().splitlines()
    if len(text) != 1:
        raise SystemExit(f"ERROR: {p} has {len(text)} lines, expected exactly 1")
    m = RX.fullmatch(text[0].strip())
    if m is None:
        raise SystemExit(f"ERROR: malformed {p}: {text[0]}")
    d, energy = float(m.group(1)), float(m.group(2))
    if not (math.isfinite(d) and math.isfinite(energy)):
        raise SystemExit(f"ERROR: NaN or Inf in {p}")
    if abs(d - i*D_STEP) > 1e-5:
        raise SystemExit(f"ERROR: wrong d in {p}: {d}, expected {i*D_STEP}")
    rows.append((d,energy))

E0 = rows[0][1]
with open("qsd_result.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["d_A", "E_eV", "DeltaE_eV"])
    for d, energy in rows:
        writer.writerow([f"{d:.1f}", f"{energy:.12f}", f"{energy-E0:.12f}"])
with open(LEGACY_FILE, "w") as f:
    for d, energy in rows:
        f.write(f"d={d:.1f} E={energy:.12f}\n")
print(f"Collected {len(rows)} points: qsd_result.csv / {LEGACY_FILE}")
'''

def qsd_compare_py() -> str:
    return r'''#!/usr/bin/env python3
"""Optional comparison of QSD-A, QSD-AC, QSD-C curves for one seed."""
import csv
from pathlib import Path
try:
    import matplotlib.pyplot as plt
except ImportError:
    raise SystemExit("matplotlib is required for plotting; CSV results need no extra packages")
for test_name in ("QSD-A", "QSD-AC", "QSD-C"):
    p = Path(test_name) / "qsd_result.csv"
    if not p.exists():
        print(f"Skip missing {p}")
        continue
    with p.open() as f:
        rows = list(csv.DictReader(f))
    plt.plot([float(r["d_A"]) for r in rows],
             [float(r["DeltaE_eV"]) for r in rows], label=test_name)
plt.xlabel("Displacement d (Angstrom)")
plt.ylabel("E(d) - E(0) (eV)")
plt.legend()
plt.tight_layout()
plt.savefig("qsd_comparison.png", dpi=250)
print("Wrote qsd_comparison.png")
'''

