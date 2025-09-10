#!/usr/bin/env python3
from pathlib import Path
from decimal import Decimal, getcontext
import json

from ase.build import bulk
from ase.io import write

element  = "Zr"
c_over_a = Decimal("1.633")

a_start = Decimal("2.00")
a_end   = Decimal("5.00")
a_step  = Decimal("0.1")

repeat = (1, 1, 1)
out_root = Path("structures")

getcontext().prec = 28

def frange(start: Decimal, stop: Decimal, step: Decimal):
    x = start
    while x <= stop + Decimal("1e-15"):
        yield x
        x += step

def main():
    out_root.mkdir(exist_ok=True)
    count = 0
    for a in frange(a_start, a_end, a_step):
        c = a * c_over_a
        d = out_root / f"a_{a:.6f}_c_{c:.6f}"
        d.mkdir(exist_ok=True)

        atoms = bulk(name=element, crystalstructure="hcp", a=float(a), c=float(c)).repeat(repeat)
        write(d / "data.lmp", atoms, format="lammps-data", atom_style="atomic")

        meta = {
            "element": element,
            "a": float(a),
            "c": float(c),
            "c_over_a": float(c_over_a),
            "repeat": repeat,
            "natoms": len(atoms),
        }
        (d / "meta.json").write_text(json.dumps(meta, indent=2))
        count += 1

    print(f"Generation complete: {count} configurations stored in {out_root}/")

if __name__ == "__main__":
    main()
