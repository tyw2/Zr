#!/bin/bash
#SBATCH --account=def-belandl1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --time=0-00:59


from pathlib import Path
import subprocess
import re
import json
import os

LMP_EXE = "/home/yiwen1/scratch/MLIP3_mlpno_lmpyes/interface-lammps-mlip-3/lmp_mpi"
structures_root = Path("/home/yiwen1/scratch/CAC/physical_properties/EOS/HCP/configuration/structures/")

PAIR_STYLE_LINES = [
    "pair_style mlip load_from=/home/yiwen1/scratch/CAC/fit_potential/delete_335/mtp_26.almtp",
    "pair_coeff * *",
]

THERMO_EVERY = 1
skip_completed = True

def make_in_lammps():
    lines = [
        "units           metal",
        "atom_style      atomic",
        "boundary        p p p",
        "read_data       data.lmp",
        "mass 1 91.224",
        "",
    ]
    lines += PAIR_STYLE_LINES
    lines += [
        "",
        "neighbor        1.0 bin",
        "neigh_modify    delay 0 every 1 check yes",
        "thermo_style    custom step pe ke etotal press lx ly lz",
        f"thermo          {THERMO_EVERY}",
        "",
        "min_style       cg",
        "min_modify      line quadratic",
        "run             0",
        "",
    ]
    return "\n".join(lines) + "\n"

def parse_energy_from_log(log_path: Path):
    if not log_path.exists():
        return None
    text = log_path.read_text(errors="ignore")

    header = None
    for line in text.splitlines():
        if re.search(r"\bStep\b", line) and re.search(r"\betotal\b|\bTotEng\b|\bPe\b|\bPotEng\b", line, re.IGNORECASE):
            header = line.split()

    data_lines = []
    if header:
        for line in text.splitlines():
            parts = line.split()
            if len(parts) == len(header):
                data_lines.append(parts)

    if header and data_lines:
        last = data_lines[-1]
        cols = {name.lower(): idx for idx, name in enumerate(header)}
        for key in ("etotal", "toteng", "pe", "poteng"):
            if key in cols:
                return float(last[cols[key]])

    m = re.findall(r"etotal\s*=\s*([-\d\.Ee+]+)", text, re.IGNORECASE)
    if m:
        return float(m[-1])
    m = re.findall(r"\bpe\s*=\s*([-\d\.Ee+]+)", text, re.IGNORECASE)
    if m:
        return float(m[-1])

    return None

def main():
    root = Path.cwd()
    if not structures_root.exists():
        raise SystemExit(f"Directory not found: {structures_root} (please run the script to generate structures first)")

    task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))
    task_cnt = int(os.environ.get("SLURM_ARRAY_TASK_COUNT", os.environ.get("SLURM_JOB_SIZE", "1")))

    all_dirs = sorted([d for d in structures_root.iterdir() if (d / "data.lmp").exists()])
    if not all_dirs:
        raise SystemExit(f"No data.lmp found under {structures_root}")

    subdirs = [all_dirs[i] for i in range(task_id, len(all_dirs), task_cnt)]
    if not subdirs:
        print(f"[INFO] Task {task_id}/{task_cnt} has no directories assigned, exiting.")
        return

    summary = root / f"energy_summary.{task_id}.txt"
    summary.write_text("a c natoms Etot(eV)\n")

    for d in subdirs:
        log = d / "log.lammps"

        if skip_completed:
            etot_prev = parse_energy_from_log(log)
            if etot_prev is not None:
                a = c = natoms = None
                meta_path = d / "meta.json"
                if meta_path.exists():
                    meta = json.loads(meta_path.read_text())
                    a = meta.get("a")
                    c = meta.get("c")
                    natoms = meta.get("natoms")
                else:
                    try:
                        a = float(d.name.split("_")[1])
                        c = float(d.name.split("_")[3])
                    except Exception:
                        pass
                summary.open("a").write(f"{a} {c} {natoms} {etot_prev}\n")
                print(f"[SKIP][{task_id}] {d.name} energy already present")
                continue

        (d / "in.lmp").write_text(make_in_lammps())

        try:
            cp = subprocess.run([LMP_EXE, "-in", "in.lmp"], cwd=d, check=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            (d / "run_stdout.txt").write_text(cp.stdout)
        except subprocess.CalledProcessError as e:
            (d / "run_stdout.txt").write_text(e.stdout or "")
            print(f"[WARN][{task_id}] Run failed: {d.name} (see {d/'run_stdout.txt'})")
            continue

        etot = parse_energy_from_log(log)
        if etot is None:
            print(f"[WARN][{task_id}] Cannot parse energy: {d.name} (check log.lammps / thermo_style)")
            continue

        a = c = natoms = None
        meta_path = d / "meta.json"
        if meta_path.exists():
            meta = json.loads(meta_path.read_text())
            a = meta.get("a")
            c = meta.get("c")
            natoms = meta.get("natoms")
        else:
            try:
                a = float(d.name.split("_")[1])
                c = float(d.name.split("_")[3])
            except Exception:
                pass
        summary.open("a").write(f"{a} {c} {natoms} {etot}\n")
        print(f"[OK][{task_id}] {d.name}: Etot={etot} eV")

    print(f"[DONE][{task_id}] Results written to {summary}")

if __name__ == "__main__":
    main()
