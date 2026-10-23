#!/usr/bin/env python3


import json
import shutil
from pathlib import Path
from zr_validation.settings import (ELASTIC_HCP_A, ELASTIC_HCP_C, ELASTIC_HCP_REPEAT,
                      ELASTIC_HCP_STRAIN, ELASTIC_HCP_BOX_RELAX,
                      ELASTIC_HCP_ETOL, ELASTIC_HCP_FTOL,
                      ELASTIC_HCP_MAXITER, ELASTIC_HCP_MAXEVAL,
                      ELASTIC_HCP_DMAX, TABLE_POTENTIAL)
from zr_validation.core.common import pair_block, write_file_if_missing
from zr_validation.core.structures import write_lmp, append_cfg
from zr_validation.workflows.elastic.bcc import stress_print, execute_script, single_run_script


def hcp_cell(a, c, repeat):
    import math
    nx, ny, nz = repeat
    if min(repeat) < 1:
        raise ValueError('HCP repeat must be positive')
    lx, ly, lz = nx*a, ny*math.sqrt(3)*a/2, nz*c
    xy = ny*a/2
    cell = [(lx, 0., 0.), (xy, ly, 0.), (0., 0., lz)]
    positions=[]
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for u,v,w in ((0.,0.,0.), (1/3,1/3,1/2)):
                    positions.append(((i+u)*a+(j+v)*a/2,
                                      (j+v)*math.sqrt(3)*a/2,
                                      (k+w)*c))
    return cell, positions


def prefix(potential, infile):
    return f'''units metal
atom_style atomic
boundary p p p
box tilt large
read_data {infile}
mass 1 91.224
{pair_block(potential)}
neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
thermo 50
thermo_style custom step pe pxx pyy pzz pyz pxz pxy press vol fmax fnorm lx ly lz xy
min_style cg
min_modify line quadratic dmax {ELASTIC_HCP_DMAX}
'''


def minimize():
    return (f'minimize {ELASTIC_HCP_ETOL:.12g} {ELASTIC_HCP_FTOL:.12g} '
            f'{ELASTIC_HCP_MAXITER} {ELASTIC_HCP_MAXEVAL}')


def equil_input(pot):


    relax = ('x 0.0 y 0.0 z 0.0 couple xy scalexy yes'
             if ELASTIC_HCP_BOX_RELAX == 'coupled_xy' else 'aniso 0.0')
    return (prefix(pot, 'init.lammps')
            + ''
            + minimize() + '\n'
            + ''
            + f'fix boxeq all box/relax {relax} vmax 0.001\n'
            + minimize() + '\n'
            + 'unfix boxeq\n'
            + 'write_data equilibrated.data\n'
            + 'write_dump all custom relaxed.dump id type xu yu zu fx fy fz modify sort id\n'
            + stress_print('reference.txt','ELASTIC_REF'))


def strain_input(pot, direction, sign):
    e = sign*ELASTIC_HCP_STRAIN


    if direction == 1:
        change = (f'variable newxy equal xy*({1+e:.12f})\n'
                  f'change_box all x scale {1+e:.12f} xy final ${{newxy}} remap units box')
    elif direction in (2,3):
        dim = 'y' if direction == 2 else 'z'
        change = f'change_box all {dim} scale {1+e:.12f} remap units box'
    else:
        tilt = {4:'yz',5:'xz',6:'xy'}[direction]
        length = 'lz' if direction in (4,5) else 'ly'
        change = (f'variable dtilt equal ({e:.12g})*{length}\n'
                  f'change_box all {tilt} delta ${{dtilt}} remap units box')
    label = f'strain_{direction}_{"pos" if sign>0 else "neg"}'
    return (prefix(pot,'equilibrated.data')
            + ''
            + change+'\n'
            + ''
            + minimize()+'\n'
            + stress_print(label+'.txt','ELASTIC_STRESS'))


def setup_elastic_hcp(test_name, pot, base_dir):
    dest=base_dir / test_name
    dest.mkdir(parents=True,exist_ok=True)
    settings = dict(phase='HCP', potential=str(pot), nlh=TABLE_POTENTIAL,
                    a=ELASTIC_HCP_A, c=ELASTIC_HCP_C,
                    repeat=list(ELASTIC_HCP_REPEAT),strain=ELASTIC_HCP_STRAIN,
                    box_relax=ELASTIC_HCP_BOX_RELAX, etol=ELASTIC_HCP_ETOL,
                    ftol=ELASTIC_HCP_FTOL, maxiter=ELASTIC_HCP_MAXITER,
                    maxeval=ELASTIC_HCP_MAXEVAL,dmax=ELASTIC_HCP_DMAX)
    cfg=dest/'simulation_config.json'
    if cfg.exists() and json.loads(cfg.read_text())!=settings:
        raise RuntimeError(f'{dest}: settings changed; archive old output or choose new OUTPUT_ROOT')
    write_file_if_missing(cfg,json.dumps(settings,indent=2,ensure_ascii=False)+'\n')
    cell,positions=hcp_cell(ELASTIC_HCP_A,ELASTIC_HCP_C,ELASTIC_HCP_REPEAT)
    if not (dest/'init.lammps').exists():
        write_lmp(dest/'init.lammps',cell,positions)
    if not (dest/'init.cfg').exists():
        append_cfg(dest/'init.cfg',cell,positions,'initial_HCP_elastic')
    write_file_if_missing(dest/'in.equilibrate.lmp',equil_input(str(pot)))
    for direction in range(1,7):
        for sign in (-1,1):
            name=f'in.strain_{direction}_{"pos" if sign>0 else "neg"}.lmp'
            write_file_if_missing(dest/name,strain_input(str(pot),direction,sign))
    write_file_if_missing(dest/'execute_elastic.sh',execute_script(),executable=True)
    write_file_if_missing(dest/'run.sh',single_run_script(dest),executable=True)
    local=dest/'collect_energy.py'
    if not local.exists():
        shutil.copy2(Path(__file__).with_name('post.py'),local)
    write_file_if_missing(dest/'collect_energy.sh',
        '#!/bin/bash\nset -euo pipefail\ncd "$(dirname "$0")"\npython3 collect_energy.py\n',
        executable=True)
    print(f'  [elastic] {dest}: {len(positions)} atoms, ±{ELASTIC_HCP_STRAIN:g} strain')
