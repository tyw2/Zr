#!/usr/bin/env python3


import os
import shutil
from pathlib import Path
from zr_validation.settings import *
from zr_validation.config_check import ensure_configured
from zr_validation.core.common import write_file
from zr_validation.core.setup_tests import setup_one_potential
from zr_validation.workflows.neb import detect_start_index
from zr_validation.workflows.melting.inputs import melt_temps
from zr_validation.workflows.melting.manage import melt_submit_manager_py, melt_collect_manager_py
from zr_validation.core.static_summary import generate_summary_py
from zr_validation.core import scheduler
from zr_validation.core.scheduler import (configure_task_specs,ask_recalculate_mode,find_potentials,make_all_tasks,
                       pack_tasks,make_batch_script,make_submit_all,make_collect_all)

def category_label(selected_tests):
    if selected_tests == GSFE_TESTS: return "gsfe"
    if selected_tests == NEB_TESTS: return "neb"
    if selected_tests == QSD_TESTS: return "qsd"
    if selected_tests == MELT_TESTS: return "melting"
    if selected_tests == EOS_TESTS: return "eos"
    if selected_tests == COHESIVE_TESTS: return "cohesive"
    if selected_tests == ELASTIC_TESTS: return "elastic"
    if selected_tests == FORMATION_TESTS: return "formation"
    return "all"

def main(selected_tests=None, seed_filter=None, label=None):
    ensure_configured(selected_tests)
    print(f"Potential root : {POT_ROOT}")
    print(f"Output root    : {OUTPUT_ROOT}")
    print(f"GSFE tests     : {len(GSFE_TESTS)}")
    print(f"NEB tests      : {len(NEB_TESTS)}")
    print(f"QSD tests      : {len(QSD_TESTS)}")
    print(f"Melting tests  : {len(MELT_TESTS)}")
    print(f"EOS tests      : {len(EOS_TESTS)}")
    print(f"Cohesive tests : {len(COHESIVE_TESTS)}")
    print(f"Elastic tests  : {len(ELASTIC_TESTS)}")
    print(f"Formation tests: {len(FORMATION_TESTS)}")
    print(f"NEB images     : {NEB_NREP}")
    print(f"Batch packing  : task-level first-fit packing ({BATCH_CORES} cores/job)")
    print()


    if OUTPUT_ROOT.exists():
        print(f"[incremental] {OUTPUT_ROOT} already exists — keeping existing data.")
    else:
        OUTPUT_ROOT.mkdir(parents=True)

    pots = find_potentials(seed_filter)
    if not pots:
        print("No potentials found. Nothing to do.")
        return

    print(f"Found {len(pots)} potential(s).")


    if selected_tests is None:
        selected_tests = ask_recalculate_mode()
    if not selected_tests:
        print("\nNothing selected. Exiting.")
        return
    configure_task_specs(selected_tests)
    print(f"Selected tests: {sorted(selected_tests)}")
    print(f"Task core reservations/potential: {scheduler.CORES_PER_POT}")
    if selected_tests & GSFE_TESTS:
        print("GSFE auto-detected configuration counts / allocated cores:")
        for name in GSFE_TEST_CONFIG:
            print(f"  {name:18s}: {scheduler.GSFE_CONFIG_COUNTS[name]:4d} configs -> {scheduler.GSFE_CORES[name]:3d} cores")
    if selected_tests & QSD_TESTS:
        print("QSD: A=51, AC=31, C=51 points; no relaxation" if not QSD_RELAX_OTHERS
              else "QSD: A=51, AC=31, C=51 points; constrain-and-relax others")
        print(f"QSD lattice mode: {QSD_LATTICE_MODE}, a={QSD_A}, parallel={QSD_PARALLEL}")
    if selected_tests & EOS_TESTS:
        print(f"EOS: BCC/HCP, a={EOS_A_START}..{EOS_A_END} step {EOS_A_STEP}; {EOS_CORES} cores/test")
    if selected_tests & COHESIVE_TESTS:
        print("Cohesive: BCC/FCC/HCP, fixed-box CG + isolated atom")
    if selected_tests & ELASTIC_TESTS:
        if 'Elastic-BCC' in selected_tests:
            print(f"BCC Elastic: a={ELASTIC_BCC_A}, repeat={ELASTIC_BCC_REPEAT}, "
                  f"±{ELASTIC_BCC_STRAIN:g} strain, box/relax={ELASTIC_BCC_BOX_RELAX}")
        if 'Elastic-FCC' in selected_tests:
            print(f"FCC Elastic: a={ELASTIC_FCC_A}, repeat={ELASTIC_FCC_REPEAT}, "
                  f"±{ELASTIC_FCC_STRAIN:g} strain, box/relax={ELASTIC_FCC_BOX_RELAX}")
        if 'Elastic-HCP' in selected_tests:
            print(f"HCP Elastic: a={ELASTIC_HCP_A}, c={ELASTIC_HCP_C}, repeat={ELASTIC_HCP_REPEAT}, "
                  f"±{ELASTIC_HCP_STRAIN:g} strain, box/relax={ELASTIC_HCP_BOX_RELAX}")
    if selected_tests & FORMATION_TESTS:
        print("Formation: HCP monovac / BO SIA / basal di-vac / out-of-plane di-vac")
        print("Di-vac binding computed with same-cell monovac reference; see README for size warnings")
    if selected_tests & MELT_TESTS:
        print(f"Melting: {len(melt_temps())} temps ({MELT_T_MIN}..{MELT_T_MAX} K, {MELT_T_STEP} K), "
              f"prep ensemble={MELT_PREP_LIQUID_ENSEMBLE}")
        print("Melting resources: 192 MPI ranks / temperature, independent job arrays; "
              "NOT packed with GSFE/NEB/QSD.")


    starts = {name: detect_start_index(cfg["src"])
              for name, cfg in NEB_TEST_CONFIG.items() if name in selected_tests}
    if selected_tests & NEB_TESTS:
        print("\nNEB image start indices (auto-detected):")
        for k, v in starts.items():
            print(f"  {k:18s}: start={v}  (images {v}..{v + NEB_NREP - 1})")


    print("\n--- Setting up test directories ---")
    for idx, pot_file in pots:
        pot_work_dir = OUTPUT_ROOT / f"seed_{idx}"
        pot_work_dir.mkdir(parents=True, exist_ok=True)
        setup_one_potential(pot_file, idx, pot_work_dir, starts, selected_tests)

    if label is None:
        label = category_label(selected_tests)


    seed_indices = [idx for idx, _ in pots]
    all_tasks, skipped = make_all_tasks(seed_indices, selected_tests)

    if skipped:
        print(f"\n--- Skipping {len(skipped)} task(s) ---")
        for s in skipped:
            print(f"  [skip] {s}")


    if selected_tests & MELT_TESTS:
        write_file(OUTPUT_ROOT / "submit_melting.py", melt_submit_manager_py(), executable=True)
        write_file(OUTPUT_ROOT / "collect_melting.py", melt_collect_manager_py(), executable=True)

    if selected_tests & (EOS_TESTS | COHESIVE_TESTS | ELASTIC_TESTS | FORMATION_TESTS):
        write_file(OUTPUT_ROOT / "collect_static.py", generate_summary_py(), executable=True)

    batch_dirs = []
    if all_tasks:

        batch_root = OUTPUT_ROOT / "jobs" / label
        batch_root.mkdir(parents=True, exist_ok=True)
        for old_batch in sorted(batch_root.glob("batch_*")):
            if old_batch.is_dir():
                shutil.rmtree(old_batch)
                print(f"  [removed old batch] {old_batch.name}")

        batches = pack_tasks(all_tasks, BATCH_CORES)
        print(f"\n=== Creating {len(batches)} batch job(s) for {len(all_tasks)} pending task(s) ===")
        for bidx, batch_tasks in enumerate(batches, start=1):
            batch_dir = batch_root / f"batch_{bidx:02d}"
            batch_dir.mkdir(parents=True, exist_ok=True)
            script = make_batch_script(bidx, batch_tasks, batch_dir)
            write_file(batch_dir / "run_batch.sh", script, executable=True)
            batch_dirs.append(batch_dir)
            n_active = sum(t["cores"] for t in batch_tasks)
            task_summary = ", ".join(
                f"seed_{t['seed']}/{t['test']}({t['cores']})" for t in batch_tasks
            )
            print(f"  batch_{bidx:02d}: {n_active}/{BATCH_CORES} cores active")
            print(f"    {task_summary}")
        write_file(OUTPUT_ROOT / f"submit_{label}.sh", make_submit_all(batch_dirs), executable=True)
        write_file(OUTPUT_ROOT / f"collect_{label}.sh", make_collect_all(seed_indices, selected_tests), executable=True)
    elif not (selected_tests & MELT_TESTS):
        print("\nNo tasks to run. Nothing to submit.")
        return

    print("\nDone.")
    if batch_dirs:
        print(f"Submit short tests: cd {OUTPUT_ROOT} && ./submit_{label}.sh")
        print(f"Collect results:    cd {OUTPUT_ROOT} && ./collect_{label}.sh")
        if selected_tests & (EOS_TESTS | COHESIVE_TESTS | ELASTIC_TESTS | FORMATION_TESTS):
            print(f"EOS/Cohesive/Elastic/Formation summary: cd {OUTPUT_ROOT} && python3 collect_static.py")
    if selected_tests & MELT_TESTS:
        print(f"Melting preview:             cd {OUTPUT_ROOT} && python3 submit_melting.py --seed {seed_indices[0]}")
        print(f"Submit one seed melting:     cd {OUTPUT_ROOT} && python3 submit_melting.py --seed {seed_indices[0]} --execute")
        print(f"Submit ALL melting:         cd {OUTPUT_ROOT} && python3 submit_melting.py --all --execute")
        print(f"After completion:           cd {OUTPUT_ROOT} && python3 collect_melting.py")

