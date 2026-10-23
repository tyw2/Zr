#!/usr/bin/env python3


from zr_validation.settings import *
from zr_validation.core.common import write_file, write_file_if_missing
from zr_validation.workflows.melting.inputs import melt_temps,melt_interface_name,melt_prep_input,melt_scan_input
from zr_validation.workflows.melting.fit_template import MELT_ANALYSIS_FUNCTIONS

def melt_slurm_header(test_name: str, d: Path, step: str, nproc: int, walltime: str):
    job_label = "%A_%a" if step == "scan" else "%j"
    return f'''#!/bin/bash
#SBATCH --account={SLURM_ACCOUNT}
#SBATCH --job-name={test_name.replace("Melt-", "m")}_{step}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={nproc}
#SBATCH --cpus-per-task=1
#SBATCH --time={walltime}
#SBATCH --output={d}/{step}-{job_label}.out
#SBATCH --error={d}/{step}-{job_label}.err
'''

def melt_prep_sh(test_name: str, d: Path) -> str:
    return (melt_slurm_header(test_name, d, "prep", BATCH_CORES, MELT_WALLTIME)
        + f'''
set -euo pipefail
module purge
module load StdEnv/2023
module load gcc/12.3
module load openmpi/4.1.5
export OMP_NUM_THREADS=1
cd "{d}"
if [ -s "{melt_interface_name(test_name)}" ]; then
    echo "[skip] Existing interface data: {melt_interface_name(test_name)}"
    exit 0
fi
mpirun -np "${{SLURM_NTASKS}}" "{LMP_EXE}" -in in_prepare.lmp -log prepare.log
''')

def melt_scan_sh(test_name: str, d: Path) -> str:
    n = len(melt_temps())
    return (melt_slurm_header(test_name, d, "scan", BATCH_CORES, MELT_WALLTIME)
        + f'''

set -euo pipefail
module purge
module load StdEnv/2023
module load gcc/12.3
module load openmpi/4.1.5
export OMP_NUM_THREADS=1
cd "{d}"
T=$(( {MELT_T_MIN} + SLURM_ARRAY_TASK_ID * {MELT_T_STEP} ))
OUTDIR="scan/${{T}}K"
mkdir -p "$OUTDIR"
if [[ -s "$OUTDIR/rdf.txt" && -s "$OUTDIR/coexist.dat" ]]; then
    echo "[skip] $T K already complete"
    exit 0
fi
if [[ ! -s "{melt_interface_name(test_name)}" ]]; then
    echo "[ERROR] interface data missing, run prepare first" >&2
    exit 2
fi
mpirun -np "${{SLURM_NTASKS}}" "{LMP_EXE}" -var T "$T" -var outdir "$OUTDIR" -in in_scan.lmp -log "$OUTDIR/log.lammps"
''')

def melt_post_sh(test_name: str, d: Path) -> str:
    return (melt_slurm_header(test_name, d, "post", 1, "0-00:30")
      + f'''
set -euo pipefail
module purge
module load StdEnv/2023
module load scipy-stack
cd "{d}"
python3 melt_postprocess.py
''')

def melt_postprocess_py() -> str:

    head = f'''#!/usr/bin/env python3
import os
import sys
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

T_MIN, T_MAX, T_STEP = {float(MELT_T_MIN)}, {float(MELT_T_MAX)}, {float(MELT_T_STEP)}
LOW_T_WINDOW = {MELT_LOW_T_WINDOW}
HIGH_T_WINDOW = {MELT_HIGH_T_WINDOW}
SMOOTH_WINDOW = {MELT_SMOOTH_WINDOW}
PEAK_SEARCH_MAX = {MELT_PEAK_SEARCH_MAX}
QUARTER_FRAC = 0.25
SCAN_DIR = "scan"
OUT_CSV = "fwqm_vs_T.csv"
OUT_FIG = "melt_point_fit.png"

def temp_list():
    return np.arange(T_MIN, T_MAX + T_STEP * 0.5, T_STEP)

def rdf_path(T):
    return os.path.join(SCAN_DIR, f"{{T:.0f}}K", "rdf.txt")

'''
    return head + MELT_ANALYSIS_FUNCTIONS + '\n\nif __name__ == "__main__":\n    analyze()\n'

def melt_metadata(test_name: str, pot_file: Path):
    return {
        "test": test_name, "pot_path": str(pot_file),
        "nlh_path": TABLE_POTENTIAL, "configuration": MELT_CONFIG[test_name],
        "prep_liquid_ensemble": MELT_PREP_LIQUID_ENSEMBLE,
        "T_solid_K": MELT_T_SOLID, "T_liquid_K": MELT_T_LIQUID,
        "pressure_bar": MELT_PEXT, "timestep_ps": MELT_TIMESTEP,
        "prep_solid_steps": MELT_PREP_SOLID_STEPS,
        "prep_liquid_steps": MELT_PREP_LIQUID_STEPS,
        "equil_steps": MELT_EQUIL_STEPS, "prod_steps": MELT_PROD_STEPS,
        "rdf_bins": MELT_RDF_BINS, "rdf_cutoff": MELT_RDF_CUTOFF,
        "rdf_comm_cutoff": MELT_COMM_CUTOFF,
        "temp_min": MELT_T_MIN, "temp_max": MELT_T_MAX, "temp_step": MELT_T_STEP,
        "smooth_window": MELT_SMOOTH_WINDOW,
        "low_window": MELT_LOW_T_WINDOW, "high_window": MELT_HIGH_T_WINDOW,
    }

def setup_melt_test(test_name: str, pot_file: Path, base_dir: Path):
    import json
    d = base_dir / test_name
    d.mkdir(parents=True, exist_ok=True)
    metadata = melt_metadata(test_name, pot_file)
    metadata_file = d / "simulation_config.json"
    if metadata_file.exists():
        old = json.loads(metadata_file.read_text(encoding="utf-8"))
        if old != json.loads(json.dumps(metadata)):
            raise RuntimeError(f"{d}: melt settings/potential changed since setup. "
                               "Choose a new OUTPUT_ROOT or archive/rename the previous Melt folder.")
    else:
        write_file(metadata_file, json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")
    write_file_if_missing(d / "in_prepare.lmp", melt_prep_input(test_name, str(pot_file)))
    write_file_if_missing(d / "in_scan.lmp", melt_scan_input(test_name, str(pot_file)))
    write_file_if_missing(d / "run_prepare.sh", melt_prep_sh(test_name, d), executable=True)
    write_file_if_missing(d / "run_scan_array.sh", melt_scan_sh(test_name, d), executable=True)
    write_file_if_missing(d / "run_post.sh", melt_post_sh(test_name, d), executable=True)
    write_file_if_missing(d / "melt_postprocess.py", melt_postprocess_py(), executable=True)
    (d / "scan").mkdir(exist_ok=True)

