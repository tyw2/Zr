import argparse
import sys
from zr_validation.settings import (
    GSFE_TESTS, NEB_TESTS, QSD_TESTS, MELT_TESTS,
    EOS_TESTS, COHESIVE_TESTS, ELASTIC_TESTS, FORMATION_TESTS
)
from zr_validation.config_check import missing_fields

CATEGORIES = {
    "gsfe": GSFE_TESTS,
    "neb": NEB_TESTS,
    "qsd": QSD_TESTS,
    "melting": MELT_TESTS,
    "eos": EOS_TESTS,
    "cohesive": COHESIVE_TESTS,
    "elastic": ELASTIC_TESTS,
    "formation": FORMATION_TESTS,
}


def select_tests(category, phase=None, case=None):
    selected = (set().union(*CATEGORIES.values()) if category == "all"
                else set(CATEGORIES[category]))
    if phase:
        if category not in ("elastic", "melting", "cohesive", "eos"):
            raise ValueError("--phase is valid for elastic, melting, cohesive, eos only")
        tag = {"elastic": "Elastic", "melting": "Melt", "cohesive": "Cohesive", "eos": "EOS"}[category]
        selected &= {f"{tag}-{phase}"}
        if not selected:
            raise ValueError(f"Unsupported --phase {phase} for {category}")
    if case:
        if category != "formation":
            raise ValueError("--case is valid for formation only")
        selected &= {"Formation-" + case}
    return selected


def main(argv=None):
    parser = argparse.ArgumentParser(prog="zrvalidate")
    sub = parser.add_subparsers(dest="command", required=True)
    for command in ("generate", "check"):
        p = sub.add_parser(command)
        p.add_argument("--only", choices=sorted(CATEGORIES) + ["all"], default="all")
        p.add_argument("--phase", choices=["BCC", "FCC", "HCP"])
        p.add_argument("--case", choices=["1_vac", "BO", "di_vac_inplane", "di_vac_outplane"])
        if command == "generate":
            p.add_argument("--seed", type=int, action="append")
    sub.add_parser("list")
    args = parser.parse_args(argv)
    if args.command == "list":
        for name, tests in CATEGORIES.items():
            print(f"{name:12s} {len(tests):2d}  " + ", ".join(sorted(tests)))
        return 0
    try:
        selected = select_tests(args.only, args.phase, args.case)
    except ValueError as err:
        parser.error(str(err))
    missing = missing_fields(selected)
    if missing:
        print("Please replace *** in config/paths.json:", file=sys.stderr)
        for item in missing:
            print("  " + item, file=sys.stderr)
        return 2
    if args.command == "check":
        print("Configuration placeholders resolved for selected workflows.")
        return 0
    from zr_validation.core.generator import main as generate
    generate(selected_tests=selected,
             seed_filter=set(args.seed) if args.seed else None,
             label=args.only)
    return 0
