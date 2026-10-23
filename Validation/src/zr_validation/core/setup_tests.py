#!/usr/bin/env python3


import os
import shutil
from pathlib import Path
from zr_validation.settings import *
from zr_validation.core.common import write_file,write_file_if_missing,make_symlink,copy_dir_replace
from zr_validation.workflows.gsfe import gsfe_input_for_test,single_auto_gsfe_run_sh,auto_gsfe_collect_sh
from zr_validation.workflows.neb import detect_start_index,neb_in,relax_endpoint_in,endpoint_helper_py,single_neb_test_run_sh,neb_collect_sh
from zr_validation.workflows.qsd import qsd_input,qsd_execute_sh,qsd_single_run_sh,qsd_collect_py,qsd_compare_py,qsd_steps
from zr_validation.workflows.melting.jobs import setup_melt_test
from zr_validation.workflows.eos import setup_eos
from zr_validation.workflows.cohesive import setup_cohesive
from zr_validation.workflows.elastic.bcc import setup_elastic_bcc
from zr_validation.workflows.elastic.fcc import setup_elastic_fcc
from zr_validation.workflows.elastic.hcp import setup_elastic_hcp
from zr_validation.workflows.formation.workflow import setup_formation

def is_task_done(test_dir: Path, test_name: str) -> bool:


    if test_name in NEB_TESTS:
        return False
    if test_name in FORMATION_TESTS:
        return ((test_dir / "formation_result.json").is_file() and
                (test_dir / "formation_result.csv").is_file())
    if test_name in EOS_TESTS:
        return (test_dir / "eos_result.csv").is_file() and (test_dir / "eos_result.csv").stat().st_size>0
    if test_name in ELASTIC_TESTS:
        return (test_dir / "elastic_result.json").is_file() and (test_dir / "elastic_result.csv").is_file()
    if test_name in COHESIVE_TESTS:
        return (test_dir / "cohesive_result.json").is_file() and (test_dir / "cohesive_result.json").stat().st_size>0
    if test_name in QSD_TESTS:
        result = test_dir / "qsd_result.csv"
        if not result.is_file():
            return False
        with result.open() as f:
            lines = [line.strip() for line in f if line.strip()]
        return (len(lines) == qsd_steps(test_name)+2
                and lines[0].startswith("d_A,E_eV,DeltaE_eV"))
    result_file = test_dir / "gsfe_result.txt"
    return result_file.exists() and result_file.stat().st_size > 0

def setup_one_potential(pot_file: Path, seed_idx: int, base_dir: Path, starts: dict, selected_tests: set):
    pot_path_str = str(pot_file)
    print(f"  seed_{seed_idx}: {pot_file.name}")


    for test_name, cfg in GSFE_TEST_CONFIG.items():
        if test_name not in selected_tests:
            continue
        d = base_dir / test_name

        write_file_if_missing(
            d / "in.gsfe",
            gsfe_input_for_test(test_name, pot_path_str)
        )

        write_file_if_missing(
            d / "run.sh",
            single_auto_gsfe_run_sh(str(d), test_name)
        )

        write_file_if_missing(
            d / "collect_energy.sh",
            auto_gsfe_collect_sh(test_name)
        )

        for f in ["run.sh", "collect_energy.sh"]:
            (d / f).chmod(0o755)

        if cfg["src"].exists():
            make_symlink(cfg["src"], d / "configuration")
        else:
            print(f"    [warn] {cfg['src']} not found, configuration symlink skipped.")


    for test_name, cfg in NEB_TEST_CONFIG.items():
        if test_name not in selected_tests:
            continue
        d = base_dir / test_name
        s = starts[test_name]
        final_idx = s + NEB_NREP - 1

        copy_dir_replace(cfg["src"], d / "replica")
        write_file(d / "endpoint_helper.py",
                   endpoint_helper_py(NEB_NREP, cfg["start_line"], cfg["end_line"], s),
                   executable=True)
        write_file(d / "in_relax_initial.lmp",
                   relax_endpoint_in(pot_path_str, s, "initial"))
        write_file(d / "in_relax_final.lmp",
                   relax_endpoint_in(pot_path_str, final_idx, "final"))
        write_file(d / cfg["input"], neb_in(pot_path_str, s, cfg["title"]))
        write_file(d / "run.sh",
                   single_neb_test_run_sh(NEB_NREP, cfg["input"]), executable=True)
        write_file_if_missing(d / "collect_energy.sh", neb_collect_sh())
        (d / "collect_energy.sh").chmod(0o755)


    for test_name in QSD_CONFIG:
        if test_name not in selected_tests:
            continue
        d = base_dir / test_name
        write_file_if_missing(d / "qsd.in", qsd_input(test_name, pot_path_str))
        write_file_if_missing(d / "execute_qsd.sh", qsd_execute_sh(test_name), executable=True)
        write_file_if_missing(d / "run.sh", qsd_single_run_sh(test_name), executable=True)
        write_file_if_missing(d / "collect_energy.py", qsd_collect_py(test_name), executable=True)
        write_file_if_missing(
            d / "collect_energy.sh",
            "#!/bin/bash\nset -euo pipefail\ncd \"$(dirname \"$0\")\"\npython3 collect_energy.py\n",
            executable=True,
        )
        for name in ("run.sh", "execute_qsd.sh", "collect_energy.sh", "collect_energy.py"):
            (d / name).chmod(0o755)
    if selected_tests & QSD_TESTS:
        write_file_if_missing(base_dir / "plot_qsd.py", qsd_compare_py(), executable=True)


    for test_name in sorted(EOS_TESTS & selected_tests):
        setup_eos(test_name, pot_file, base_dir)


    for test_name in sorted(COHESIVE_TESTS & selected_tests):
        setup_cohesive(test_name, pot_file, base_dir)


    for test_name in sorted(ELASTIC_TESTS & selected_tests):
        if test_name == 'Elastic-BCC':
            setup_elastic_bcc(test_name, pot_file, base_dir)
        elif test_name == 'Elastic-FCC':
            setup_elastic_fcc(test_name, pot_file, base_dir)
        elif test_name == 'Elastic-HCP':
            setup_elastic_hcp(test_name, pot_file, base_dir)
        else:
            raise ValueError(f'Unknown elastic phase: {test_name}')


    for test_name in sorted(FORMATION_TESTS & selected_tests):
        setup_formation(test_name, pot_file, base_dir)


    for test_name in MELT_CONFIG:
        if test_name in selected_tests:
            setup_melt_test(test_name, pot_file, base_dir)

def count_input_configs(config_dir: Path) -> int:

    if not config_dir.exists():
        return 0
    return sum(1 for p in config_dir.glob("input.pos*") if p.is_file())

def detect_all_gsfe_counts() -> dict:


    counts = {}

    for test_name, cfg in GSFE_TEST_CONFIG.items():
        n = count_input_configs(cfg["src"])

        if n <= 0:
            raise RuntimeError(
                f"GSFE config detection failed for {test_name}: "
                f"no input.pos* found in {cfg['src']}"
            )

        counts[test_name] = n

    return counts

def cores_for_gsfe(nconfig: int) -> int:

    return min(nconfig, BATCH_CORES)

