#!/usr/bin/env python3


from itertools import product
from math import sqrt
from zr_validation.core.structures import crystal, write_lmp, append_cfg


def vecmul(frac, cell):
    return tuple(sum(frac[i] * cell[i][j] for i in range(3)) for j in range(3))


def sub(a, b):
    return tuple(a[k] - b[k] for k in range(3))


def norm(v):
    return sqrt(sum(x*x for x in v))


def mic_vector(a, b, cell):

    delta = sub(b, a)
    return min((tuple(delta[k] + sum(s[j]*cell[j][k] for j in range(3))
                      for k in range(3)) for s in product((-1, 0, 1), repeat=3)),
               key=norm)


def closest_to_center(positions, cell):
    center = tuple(sum(row[j] for row in cell) * 0.5 for j in range(3))
    return min(range(len(positions)), key=lambda i: norm(sub(positions[i], center)))


def select_pair(positions, cell, mode, cutoff):

    i = closest_to_center(positions, cell)
    neighbors = []
    for j in range(len(positions)):
        if j == i:
            continue
        v = mic_vector(positions[i], positions[j], cell)
        r = norm(v)
        if r > 1e-8 and r < cutoff:
            neighbors.append((j, v, r))
    if not neighbors:
        raise ValueError(f'No neighbors within {cutoff} A; change cutoff or supercell')
    if mode == 'basal':

        j, vector, r = min(neighbors, key=lambda e: (round(abs(e[1][2]), 8), e[2], e[0]))
    elif mode == 'out_of_plane':

        j, vector, r = min(neighbors, key=lambda e: (-round(abs(e[1][2])/e[2], 8), e[2], e[0]))
    else:
        raise ValueError(mode)
    return i, j, vector, r


def make_bo(positions, cell, repeat, a, c, mode, half_separation, min_separation):


    new_positions = list(positions)
    i0 = closest_to_center(positions, cell)
    meta = {'sia_mode': mode}
    if mode == 'dumbbell':
        direction = (1., 0., 0.)
        old = positions[i0]
        new_positions[i0] = tuple(old[k] - half_separation*direction[k] for k in range(3))
        new_positions.append(tuple(old[k] + half_separation*direction[k] for k in range(3)))
        meta.update(original_index=i0, separation_A=2*half_separation)
    elif mode == 'octa':
        nx, ny, nz = repeat


        base_sites = ((2/3, 1/3, 0.), (1/3, 2/3, 0.),
                      (2/3, 1/3, .5), (1/3, 2/3, .5))
        primitive = [(a, 0., 0.), (-a/2, sqrt(3)*a/2, 0.), (0., 0., c)]
        center = vecmul((.5, .5, .5), cell)
        candidates = []
        for ix, iy, iz in product(range(nx), range(ny), range(nz)):
            for basis in base_sites:
                f = (ix+basis[0], iy+basis[1], iz+basis[2])
                p = vecmul(f, primitive)
                d_nearest = min(norm(mic_vector(p, q, cell)) for q in positions)
                if d_nearest >= min_separation:
                    candidates.append((norm(sub(p, center)), p, f, d_nearest))
        if not candidates:
            raise ValueError('BO site list has no unoccupied candidates; verify fractional convention')
        _, site, frac, nearest = min(candidates, key=lambda x: x[0])
        new_positions.append(site)
        meta.update(site_fractional_primitive=frac, initial_nearest_distance_A=nearest,
                    caution='Candidate BO site from supplied script; validate BO/O identity after CG')
    else:
        raise ValueError('FORMATION_BO_MODE must be octa or dumbbell')
    return new_positions, meta


def generate_structures(folder, cfg):

    cell, bulk = crystal('HCP', cfg['a'], c=cfg['c'], repeat=tuple(cfg['repeat']))
    n = len(bulk)
    mode = cfg['kind']
    i0 = closest_to_center(bulk, cell)
    mono = None
    details = {}
    if mode == 'vacancy':
        defect = [p for j, p in enumerate(bulk) if j != i0]
        details = {'removed_atom_indices_0based': [i0]}
        dN = -1
    elif mode == 'bo':
        defect, details = make_bo(bulk, cell, tuple(cfg['repeat']), cfg['a'], cfg['c'],
                                  cfg['bo_mode'], cfg['dumbbell_disp'], cfg['min_bo_separation'])
        dN = +1
    elif mode in ('basal', 'out_of_plane'):
        i, j, v, r = select_pair(bulk, cell, mode, cfg['nn_cut'])
        defect = [p for k, p in enumerate(bulk) if k not in (i, j)]
        mono = [p for k, p in enumerate(bulk) if k != i]
        details = {'removed_atom_indices_0based': [i, j],
                   'pair_vector_MIC_A': v, 'pair_distance_A': r,
                   'note': 'same-cell monovacancy used for E_bind; in/out selection by c projection'}
        dN = -2
    else:
        raise ValueError(f'Unknown formation kind: {mode}')
    assert len(defect) == n + dN
    if mono is not None:
        assert len(mono) == n - 1
    structs = dict(bulk=bulk, defect=defect)
    if mono is not None:
        structs['monovac'] = mono
    for name, positions in structs.items():
        write_lmp(folder/name/'data.lmp', cell, positions)
        cfg_path = folder/f'initial_{name}.cfg'
        if not cfg_path.exists():
            append_cfg(cfg_path, cell, positions,
                       f'formation_{cfg["kind"]}_{name}')
    details['N_bulk'] = n
    details['N_defect'] = len(defect)
    return details, len(structs)
