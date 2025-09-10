#!/usr/bin/env python3
from pathlib import Path
from decimal import Decimal, getcontext
import json

from ase.build import bulk
from ase.io import write

element  = "Zr"              
a_start  = Decimal("2.00")   
a_end    = Decimal("5.00")    
a_step   = Decimal("0.10")    
repeat   = (1, 1, 1)          
out_root = Path("structures") 

getcontext().prec = 28

def frange(start: Decimal, stop: Decimal, step: Decimal):
    """Generate a decimal range from start to stop with step."""
    x = start
    while x <= stop + Decimal("1e-15"):
        yield x
        x += step

def main():
    out_root.mkdir(exist_ok=True)
    count = 0

    for a in frange(a_start, a_end, a_step):
        d = out_root / f"a_{a:.6f}"
        d.mkdir(exist_ok=True)

        atoms = bulk(name=element, crystalstructure="bcc", a=float(a)).repeat(repeat)

        write(d / "data.lmp", atoms, format="lammps-data", atom_style="atomic")

        meta = {
            "element": element,
            "structure": "bcc",
            "a": float(a),
            "repeat": repeat,
            "natoms": len(atoms),
        }
        (d / "meta.json").write_text(json.dumps(meta, indent=2))

        count += 1

    print(f"Generation complete: {count} bcc configurations stored in {out_root}/")

if __name__ == "__main__":
    main()
