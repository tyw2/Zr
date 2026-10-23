#!/usr/bin/env python3


def melt_submit_manager_py() -> str:
    return r'''#!/usr/bin/env python3
"""Explicit melting job submission: prepare -> afterok temperature array -> afterany post.
Usage: python3 submit_melting.py --seed 0 [--phase BCC]
       python3 submit_melting.py --all
       python3 submit_melting.py --all
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
PHASES=("Melt-BCC", "Melt-FCC", "Melt-HCP")

def submit(args):
    selected_seeds = [p for p in sorted(ROOT.glob("seed_*"), key=lambda p: int(p.name[5:]))
                      if p.is_dir() and p.name[5:].isdigit()]
    if args.seed:
        wanted=set(args.seed)
        selected_seeds=[p for p in selected_seeds if int(p.name[5:]) in wanted]
        missing=wanted-{int(p.name[5:]) for p in selected_seeds}
        if missing: raise SystemExit(f"Seed(s) not found: {sorted(missing)}")
    selected_phases=("Melt-"+args.phase,) if args.phase else PHASES
    todo=[]
    for seed in selected_seeds:
        for phase in selected_phases:
            d=seed/phase
            if (d/"simulation_config.json").exists():
                if (d/"melt_result.json").exists() and not args.rerun:
                    print(f"[skip] {d}: existing melt_result.json")
                    continue
                todo.append(d)
    if not todo:
        print("No pending phase(s) found; generator may not have selected melting.")
        return
    print(f"Selected {len(todo)} phase(s), each: 1 preparation + 201 full-node temperatures + 1 post job.")
    print("No work is submitted unless --execute is given.")
    for d in todo:
        print(f"  {d.relative_to(ROOT)}")
    if not args.execute:
        return
    with (ROOT/"melting_jobs.jsonl").open("a",encoding="utf-8") as out:
        for d in todo:
            ids={"test_dir": str(d)}
            for stage, script, depends in (
                ("prepare", "run_prepare.sh", None),
                ("scan", "run_scan_array.sh", "afterok:prepare"),
                ("post", "run_post.sh", "afterany:scan"),
            ):
                cmd=["sbatch", "--parsable"]
                if depends:
                    kind, prev=depends.split(":")
                    cmd.append(f"--dependency={kind}:{ids[prev]}")
                cmd.append(str(d/script))
                try:
                    ids[stage]=subprocess.check_output(cmd, stderr=subprocess.STDOUT).decode().strip().split(";")[0]
                except (OSError,subprocess.CalledProcessError) as exc:
                    msg=(exc.output.decode(errors="replace") if isinstance(exc,subprocess.CalledProcessError)
                         else str(exc))
                    print(f"ERROR submitting {d}/{script}: {msg}", file=sys.stderr)
                    print("Previous stage jobs may have been submitted! Record: "+str(ids),file=sys.stderr)
                    sys.exit(2)
            out.write(json.dumps(ids)+"\n")
            out.flush()
            print(f"[submitted] {d.relative_to(ROOT)} "
                  f"prep={ids['prepare']}, array={ids['scan']}, post={ids['post']}")

if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--seed", type=int, action="append", help="One or more --seed N selectors")
    g.add_argument("--all", action="store_true", help="Select all seeds (explicitly)")
    ap.add_argument("--phase", choices=("BCC","FCC","HCP"))
    ap.add_argument("--execute", action="store_true", help="Actually submit SLURM jobs; default preview")
    ap.add_argument("--rerun", action="store_true", help="Include phases with existing Tm results (scan may still skip completed temperatures)")
    submit(ap.parse_args())
'''

def melt_collect_manager_py() -> str:
    return r'''#!/usr/bin/env python3
"""Summarize available melting results across seeds and crystal prototypes."""
import csv
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
rows=[]
for fp in sorted(ROOT.glob("seed_*/Melt-*/melt_result.json")):
    result=json.loads(fp.read_text(encoding="utf-8"))
    rows.append({
        "seed": fp.parent.parent.name,
        "phase": fp.parent.name.removeprefix("Melt-"),
        "Tm_estimate_K": result["Tm_estimate_K"],
        "fit_status": result["fit_status"],
        "temperature_points": result["temperature_points"],
    })
fp=ROOT/"melting_summary.csv"
with fp.open("w", newline="",encoding="utf-8") as f:
    writer=csv.DictWriter(f,fieldnames=["seed","phase","Tm_estimate_K","fit_status","temperature_points"])
    writer.writeheader()
    writer.writerows(rows)
print(f"Summarized {len(rows)} phase results -> {fp}")
'''

