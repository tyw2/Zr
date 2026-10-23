#!/usr/bin/env python3


import csv
import json
import math
import re
from pathlib import Path

RE = re.compile(r'FORM_RESULT\s+E=([^\s]+)\s+N=(\d+)\s+FMAX=([^\s]+)')


def collect(root):
    root = Path(root)
    conf = json.loads((root/'simulation_config.json').read_text())
    details = json.loads((root/'structure_info.json').read_text())
    data = {}
    for sub in ('bulk', 'defect') + (('monovac',) if conf['kind'] in ('basal', 'out_of_plane') else ()):
        txt = root/sub/'result.txt'
        log = root/sub/'log.lammps'
        if not txt.is_file() or not log.is_file():
            raise RuntimeError(f'Missing {sub} result or log')
        match = RE.search(txt.read_text())
        if not match or 'Total wall time:' not in log.read_text(errors='replace'):
            raise RuntimeError(f'Incomplete {sub} result/log')
        energy, count, fmax = float(match[1]), int(match[2]), float(match[3])
        if not all(math.isfinite(x) for x in (energy, fmax)):
            raise RuntimeError(f'Nonfinite {sub} energy or force')
        data[sub] = (energy, count, fmax)
    bulk_e, n, fbulk = data['bulk']
    defect_e, nd, fdefect = data['defect']
    if n != details['N_bulk'] or nd != details['N_defect']:
        raise RuntimeError(f'Atom count mismatch: {n} {nd} vs {details}')
    expected_change = 1 if conf['kind'] == 'bo' else (-1 if conf['kind'] == 'vacancy' else -2)
    if nd != n + expected_change:
        raise RuntimeError('Incorrect defect stoichiometry')
    ef = defect_e - (nd/n)*bulk_e
    row = dict(phase='HCP', defect_type=conf['kind'], case=conf['test'],
               relaxation=conf['relaxation'], N_bulk=n, N_defect=nd,
               E_bulk_eV=bulk_e, E_defect_eV=defect_e,
               formation_energy_eV=ef, bulk_fmax_eV_A=fbulk, defect_fmax_eV_A=fdefect,
               E_monovac_eV=None, monovac_fmax_eV_A=None,
               monovac_formation_energy_eV=None, divac_binding_energy_eV=None)
    if 'monovac' in data:
        e_v, nv, fv = data['monovac']
        if nv != n-1:
            raise RuntimeError('Wrong same-cell monovacancy atom count')
        efv = e_v - ((n-1)/n)*bulk_e
        row.update(E_monovac_eV=e_v, monovac_fmax_eV_A=fv,
                   monovac_formation_energy_eV=efv,
                   divac_binding_energy_eV=2*efv-ef)

    with (root/'formation_result.csv').open('w',newline='') as f:
        w = csv.DictWriter(f,fieldnames=list(row))
        w.writeheader(); w.writerow(row)
    (root/'formation_result.json').write_text(json.dumps(row,indent=2,allow_nan=False)+'\n')
    print(f'[OK] {conf["test"]}: E_form={ef:.9f} eV (N={n}, defect N={nd})')
    if row['divac_binding_energy_eV'] is not None:
        print(f'     same-cell E_binding(2V)={row["divac_binding_energy_eV"]:.9f} eV')
    return row


if __name__=='__main__':
    collect(Path(__file__).resolve().parent)
