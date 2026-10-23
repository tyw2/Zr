#!/usr/bin/env python3


import os
import shutil
from pathlib import Path
from zr_validation.settings import *
from zr_validation.workflows.qsd import qsd_steps
from zr_validation.core.setup_tests import is_task_done,detect_all_gsfe_counts,cores_for_gsfe
from zr_validation.workflows.melting.inputs import melt_temps
from zr_validation.workflows.melting.manage import melt_submit_manager_py,melt_collect_manager_py
from zr_validation.core.common import write_file

GSFE_CONFIG_COUNTS = {}

GSFE_CORES = {}

CORES_NEB = NEB_NREP

CORES_PER_POT = 0

TEST_SPECS = []

def configure_task_specs(selected_tests: set):

    global GSFE_CONFIG_COUNTS, GSFE_CORES, CORES_PER_POT, TEST_SPECS
    GSFE_CONFIG_COUNTS = detect_all_gsfe_counts() if selected_tests & GSFE_TESTS else {}
    GSFE_CORES = {name: cores_for_gsfe(n) for name,n in GSFE_CONFIG_COUNTS.items()}
    TEST_SPECS = (
        [(name, GSFE_CORES[name]) for name in GSFE_TEST_CONFIG if name in selected_tests]
        + [(name, CORES_NEB) for name in NEB_TEST_CONFIG if name in selected_tests]
        + [(name, min(QSD_PARALLEL, qsd_steps(name)+1))
           for name in QSD_CONFIG if name in selected_tests]
        + [(name, EOS_CORES) for name in sorted(EOS_TESTS & selected_tests)]
        + [(name, COHESIVE_CORES) for name in sorted(COHESIVE_TESTS & selected_tests)]
        + [(name, ELASTIC_CORES) for name in sorted(ELASTIC_TESTS & selected_tests)]
        + [(name, FORMATION_CORES[name]) for name in sorted(FORMATION_TESTS & selected_tests)]
    )
    if not (1 <= QSD_PARALLEL <= BATCH_CORES):
        raise ValueError("QSD_PARALLEL must be in range 1..BATCH_CORES")
    if not (1 <= EOS_CORES <= BATCH_CORES):
        raise ValueError("EOS_CORES must be in 1..BATCH_CORES")
    if COHESIVE_CORES != 2:
        raise ValueError("COHESIVE_CORES must be 2 (bulk and single run concurrently)")
    if abs(QSD_D_STEP - 0.1) > 1e-12:
        raise ValueError("Current QSD execute_qsd.sh uses 0.1-A increments: set QSD_D_STEP=0.1")
    if ELASTIC_BCC_BOX_RELAX not in ("iso", "aniso"):
        raise ValueError("ELASTIC_BCC_BOX_RELAX must be iso or aniso")
    if not (0 < ELASTIC_BCC_STRAIN < 0.05):
        raise ValueError("ELASTIC_BCC_STRAIN must be in (0,0.05)")
    if ELASTIC_FCC_BOX_RELAX not in ("iso", "aniso"):
        raise ValueError("ELASTIC_FCC_BOX_RELAX must be iso or aniso")
    if not (0 < ELASTIC_FCC_STRAIN < 0.05):
        raise ValueError("ELASTIC_FCC_STRAIN must be in (0,0.05)")
    if ELASTIC_HCP_BOX_RELAX not in ("coupled_xy", "aniso"):
        raise ValueError("ELASTIC_HCP_BOX_RELAX must be coupled_xy or aniso")
    if not (0 < ELASTIC_HCP_STRAIN < 0.05):
        raise ValueError("ELASTIC_HCP_STRAIN must be in (0,0.05)")
    if ELASTIC_CORES != 1:
        raise ValueError("Elastic-BCC/FCC/HCP tasks use exactly 1 core per test")
    if any(cores < 1 or cores > BATCH_CORES for cores in FORMATION_CORES.values()):
        raise ValueError("Invalid formation task cores")
    CORES_PER_POT = sum(cores for _, cores in TEST_SPECS)

def ask_yes_no(question: str) -> bool:
    while True:
        ans = input(f"{question} [y/n]: ").strip().lower()
        if ans in ("y", "yes"):
            return True
        if ans in ("n", "no"):
            return False
        print("  Please answer y or n.")

def ask_recalculate_mode() -> set:

    print()
    selected = set()
    if ask_yes_no("Run all GSFE tests (auto-detect config counts and cores)?"):
        selected |= GSFE_TESTS
    if ask_yes_no("Run NEB tests (SIA_basal/nonbasal + vacancy_basal/nonbasal)?"):
        selected |= NEB_TESTS
    if ask_yes_no("Run QSD tests (QSD-A / QSD-AC / QSD-C)?"):
        selected |= QSD_TESTS
    if ask_yes_no("Prepare melting tests (BCC / FCC / HCP; full-node 201-point scans)?"):
        selected |= MELT_TESTS
    if ask_yes_no("Run EOS tests (BCC / HCP; a=2.0..5.0 A)?"):
        selected |= EOS_TESTS
    if ask_yes_no("Run cohesive energy tests (BCC / FCC / HCP)?"):
        selected |= COHESIVE_TESTS
    if ask_yes_no("Run BCC, FCC and HCP elastic constants (0 K)?"):
        selected |= ELASTIC_TESTS
    if ask_yes_no("Run HCP defect formation energies (1_vac / BO / 2V basal / 2V out-of-plane)?"):
        selected |= FORMATION_TESTS
    return selected

def make_all_tasks(seed_indices: list, selected_tests: set):
    tasks   = []
    skipped = []

    for seed_idx in seed_indices:
        for test_name, cores in TEST_SPECS:

            if test_name not in selected_tests:
                skipped.append(f"seed_{seed_idx}/{test_name}  [not selected]")
                continue

            test_dir = OUTPUT_ROOT / f"seed_{seed_idx}" / test_name
            if test_name in (QSD_TESTS | EOS_TESTS | COHESIVE_TESTS | ELASTIC_TESTS | FORMATION_TESTS) and is_task_done(test_dir, test_name):
                skipped.append(f"seed_{seed_idx}/{test_name}  [already collected]")
                continue
            tasks.append({"seed": seed_idx, "test": test_name, "cores": cores})

    return tasks, skipped

def pack_tasks(tasks: list, capacity: int = BATCH_CORES) -> list:

    batches = []
    for task in sorted(tasks, key=lambda x: x["cores"], reverse=True):
        placed = False
        for batch in batches:
            used = sum(t["cores"] for t in batch)
            if used + task["cores"] <= capacity:
                batch.append(task)
                placed = True
                break
        if not placed:
            batches.append([task])
    return batches

def task_command(task: dict) -> str:
    seed_idx = task["seed"]
    test = task["test"]
    work_dir = OUTPUT_ROOT / f"seed_{seed_idx}" / test
    nrep = NEB_NREP

    if test in GSFE_TESTS:
        max_parallel = task["cores"]

        return f"""
(
    cd {work_dir}
    set -e
    MAX_PARALLEL={max_parallel}

    mapfile -t CONFIGS < <(find -L configuration -maxdepth 1 -type f \
        -name 'input.pos*' -printf '%p\\n' | sort -V)

    NCONFIG=${{#CONFIGS[@]}}

    if (( NCONFIG == 0 )); then
        echo "ERROR: seed_{seed_idx}/{test} found no input.pos*" >&2
        exit 1
    fi

    echo "seed_{seed_idx}/{test}: detected $NCONFIG configs"
    echo "seed_{seed_idx}/{test}: allocated $MAX_PARALLEL cores"

    count=0

    for cfg_file in "${{CONFIGS[@]}}"; do
        fname=$(basename "$cfg_file")
        step_index="${{fname#input.pos}}"

        step_dir="{work_dir}/Step_${{step_index}}"
        mkdir -p "$step_dir"

        cp in.gsfe "$step_dir/"
        cp "$cfg_file" "$step_dir/input.pos${{step_index}}"
        sed -i "s/\\${{step_index}}/$step_index/g" "$step_dir/in.gsfe"

        (
            cd "$step_dir"
            srun --exclusive --ntasks=1 --cpus-per-task=1 \
                 $LMP -in in.gsfe > lammps.out 2>&1
        ) &

        count=$((count + 1))

        if (( count % MAX_PARALLEL == 0 )); then
            wait
        fi
    done

    wait
    echo "seed_{seed_idx} {test} done ($NCONFIG configs)"
) &
"""

    if test in NEB_TESTS:
        input_file = NEB_TEST_CONFIG[test]["input"]
        return f'''
(
    cd {work_dir}
    set -e
    srun --exclusive --ntasks=1 --cpus-per-task=1 \\
         $LMP -in in_relax_initial.lmp > relax_initial.out 2>&1
    srun --exclusive --ntasks=1 --cpus-per-task=1 \\
         $LMP -in in_relax_final.lmp > relax_final.out 2>&1
    python3 endpoint_helper.py update
    srun --exclusive --ntasks={nrep} --cpus-per-task=1 \\
         $LMP -partition {nrep}x1 -in {input_file} > lammps.out 2>&1
    echo "seed_{seed_idx} {test} endpoint in-place relax + CI/CL-NEB done"
) &
'''

    if test in QSD_TESTS:
        return f'''
(
    cd {work_dir}
    set -e
    bash execute_qsd.sh
    echo "seed_{seed_idx} {test} done"
) &
'''

    if test in EOS_TESTS:
        return f'''
(
    cd "{work_dir}"
    set -e
    bash execute_eos.sh
    echo "seed_{seed_idx} {test} done"
) &
'''

    if test in COHESIVE_TESTS:
        return f'''
(
    cd "{work_dir}"
    set -e
    bash execute_cohesive.sh
    echo "seed_{seed_idx} {test} done"
) &
'''

    if test in FORMATION_TESTS:
        return f'''
(
    cd "{work_dir}"
    set -e
    bash execute_formation.sh
    echo "seed_{seed_idx} {test} done"
) &
'''

    if test in ELASTIC_TESTS:
        return f'''
(
    cd "{work_dir}"
    set -e
    bash execute_elastic.sh
    echo "seed_{seed_idx} {test} done"
) &
'''

    raise ValueError(f"Unknown test: {test}")

def make_batch_script(batch_idx: int, tasks: list, batch_dir: Path) -> str:
    used = sum(t["cores"] for t in tasks)
    task_summary = ", ".join(
        f"seed_{t['seed']}/{t['test']}({t['cores']})" for t in tasks
    )

    header = f'''#!/bin/bash
#SBATCH --account={SLURM_ACCOUNT}
#SBATCH --job-name=batch_{batch_idx:02d}
#SBATCH --nodes=1
#SBATCH --ntasks={BATCH_CORES}
#SBATCH --cpus-per-task=1
#SBATCH --time={BATCH_TIME}
#SBATCH --output={batch_dir}/slurm-%j.out
#SBATCH --error={batch_dir}/slurm-%j.err

module purge
module load StdEnv/2023
module load gcc/12.3
module load openmpi/4.1.5

set -e
LMP={LMP_EXE}


echo "=== Batch job start: $(date) ==="
echo "Active cores: {used}/{BATCH_CORES}"

'''

    body = 'task_pids=()\n' + "".join(task_command(t) + 'task_pids+=("$!")\n' for t in tasks)

    footer_lines = [
        '',
        '',
        '',
        'echo "Waiting for all tasks in this batch to complete..."',
        "failed=0",
        'for pid in "${task_pids[@]}"; do',
        '  if ! wait "$pid"; then failed=1; fi',
        'done',
        'if (( failed )); then echo "ERROR: at least one packed task failed" >&2; exit 1; fi',
        'echo "=== Batch tasks done: $(date) ==="',
        "",
        '',
        'echo "--- Collecting energies for this batch ---"',
    ]
    for t in tasks:
        collect_dir = OUTPUT_ROOT / f"seed_{t['seed']}" / t["test"]
        footer_lines.append(f"( cd {collect_dir} && bash collect_energy.sh )")
    footer_lines.append('echo "=== Batch job end: $(date) ==="')

    return header + body + "\n".join(footer_lines) + "\n"

def make_submit_all(batch_dirs: list) -> str:
    lines = ["#!/bin/bash", '', "set -e", ""]
    for bd in batch_dirs:
        lines.append(f"echo 'Submitting {bd.name}...'")
        lines.append(f"sbatch {bd}/run_batch.sh")
    lines.append("echo 'All batches submitted.'")
    return "\n".join(lines) + "\n"

def make_collect_all(seed_indices: list, selected_tests: set) -> str:
    lines = [
        "#!/bin/bash",
        '',
        "",
        f"TESTS_ROOT={OUTPUT_ROOT}",
        "",
    ]
    for idx in seed_indices:
        for sub in [
            "prism_I_1", "basal_1", "prism_I_2", "basal_2",
            "prism_II_1", "prism_II_2",
            "pyram_I_1_a", "pyram_I_1_a+c", "pyram_I_2",
            "pyram_II_1_a+c", "pyram_II_2",
            "SIA_basal", "SIA_nonbasal",
            "vacancy_basal", "vacancy_nonbasal",
            "QSD-A", "QSD-AC", "QSD-C",
            "EOS-BCC", "EOS-HCP",
            "Cohesive-BCC", "Cohesive-FCC", "Cohesive-HCP", "Elastic-BCC", "Elastic-FCC", "Elastic-HCP",
            "Formation-1_vac", "Formation-BO", "Formation-di_vac_inplane", "Formation-di_vac_outplane",
        ]:
            if sub not in selected_tests:
                continue
            lines.append(
                f"( cd $TESTS_ROOT/seed_{idx}/{sub} && bash collect_energy.sh )"
                f" 2>/dev/null || true"
            )
    lines.append("echo 'Collection complete.'")
    return "\n".join(lines) + "\n"

def find_potentials(seed_filter=None):
    results = []
    if not POT_ROOT.exists():
        print(f"[error] Potential root does not exist: {POT_ROOT}")
        return results

    for idx in (sorted(seed_filter) if seed_filter is not None else range(0, 201)):
        fname = f"seed_{idx}.almtp"
        pot_file = POT_ROOT / fname
        if pot_file.exists():
            results.append((idx, pot_file))
            continue
        pot_file = POT_ROOT / f"seed_run_{idx}" / fname
        if pot_file.exists():
            results.append((idx, pot_file))
            continue
        print(f"    [warn] {fname} not found in either location, skipping.")

    return results

