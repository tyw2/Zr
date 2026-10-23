#!/usr/bin/env python3


import csv
import json
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
COLUMNS = ['PXX', 'PYY', 'PZZ', 'PYZ', 'PXZ', 'PXY']


def load_record(path, prefix):
    if not path.exists():
        raise FileNotFoundError(f'Missing {path}')
    match = re.search(rf'{prefix}\s+([^\n]+)', path.read_text())
    if not match:
        raise ValueError(f'Incomplete result file: {path}')
    fields = dict(re.findall(r'([A-Za-z_]+)=([^\s]+)', match.group(1)))
    mandatory = ['E', 'V', 'N', 'FMAX', *COLUMNS]
    if not all(k in fields for k in mandatory):
        raise ValueError(f'Missing fields in {path}: {set(mandatory)-set(fields)}')
    record = {k: float(fields[k]) for k in mandatory}
    if not all(math.isfinite(v) for v in record.values()):
        raise ValueError(f'Nonfinite result in {path}')
    return record


def csv_write(path, rows, fieldnames=None):
    if not rows:
        return
    if fieldnames is None:
        fieldnames = list(rows[0])
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_relaxed_cfg(datafile, dest, energy, phase='BCC'):


    contents = datafile.read_text().splitlines()
    bounds={}
    tilt=[0.0,0.0,0.0]
    n=None
    start=None
    for i,line in enumerate(contents):
        bits=line.split()
        if len(bits)>=2 and bits[0].isdigit() and bits[1]=='atoms':
            n=int(bits[0])
        if len(bits)>=4 and bits[2] in ('xlo','ylo','zlo'):
            bounds[bits[2][0]]=(float(bits[0]),float(bits[1]))
        if len(bits)>=6 and bits[3:6]==['xy','xz','yz']:
            tilt=[float(v) for v in bits[:3]]
        if bits and bits[0]=='Atoms':
            start=i
    if n is None or len(bounds)!=3 or start is None:
        raise ValueError('Unsupported LAMMPS data file format for relaxed.cfg.0')
    lx,ly,lz=[bounds[k][1]-bounds[k][0] for k in 'xyz']
    xy,xz,yz=tilt
    cell=((lx,0,0),(xy,ly,0),(xz,yz,lz))
    atoms=[]
    for line in contents[start+1:]:
        bits=line.split('#')[0].split()
        if len(bits)<5:
            if atoms: break
            continue
        try:
            atid,typ=int(bits[0]),int(bits[1])
            xyz=[float(q) for q in bits[2:5]]
        except ValueError:
            if atoms: break
            continue
        if typ!=1:
            raise ValueError('Only single-type Zr data are supported')
        atoms.append((atid,xyz))
    if len(atoms)!=n or sorted(i for i,_ in atoms)!=list(range(1,n+1)):
        raise ValueError(f'Atom data mismatch: parsed {len(atoms)}, expected {n}')
    with dest.open('w') as f:
        f.write(f'BEGIN_CFG\n Size\n    {n}\n Supercell\n')
        for row in cell:
            f.write('    '+' '.join(f'{v:.12f}' for v in row)+'\n')
        f.write(' AtomData:  id type cartes_x cartes_y cartes_z\n')
        for idx,xyz in sorted(atoms):
            position=[xyz[j]-bounds['xyz'[j]][0] for j in range(3)]
            f.write(f'    {idx:8d} 0 '+' '.join(f'{v:.12f}' for v in position)+'\n')
        f.write(f' Energy\n    {energy:.12f}\n Feature   relaxed_{phase}_elastic\nEND_CFG\n')


def collect(root=ROOT):
    config=json.loads((root/'simulation_config.json').read_text())
    eps=config['strain']
    phase=config['phase'].upper()
    atoms_per_cell={'BCC':2,'FCC':4,'HCP':2}
    if phase not in atoms_per_cell:
        raise ValueError(f'Unsupported phase: {phase}')
    ref=load_record(root/'reference.txt','ELASTIC_REF')
    expected=atoms_per_cell[phase]*math.prod(config['repeat'])
    if int(ref['N'])!=expected:
        raise ValueError(f'Reference atom count {int(ref["N"])} != {expected}')
    if eps<=0:
        raise ValueError('strain must be positive')
    raw=[[0.0]*6 for _ in range(6)]
    signed=[]
    max_force=ref['FMAX']
    for j in range(6):
        stress={}
        for sign,label in [(-1,'neg'),(1,'pos')]:
            row=load_record(root/f'strain_{j+1}_{label}.txt','ELASTIC_STRESS')
            if row['N']!=ref['N']:
                raise ValueError(f'Atom count changed for strain {j+1} {label}')
            max_force=max(max_force,row['FMAX'])
            signed.append(dict(direction=j+1,sign=label,**row))
            stress[sign]=row
        for i,name in enumerate(COLUMNS):


            raw[i][j]=-(stress[1][name]-stress[-1][name])*1.0e-4/(2*eps)
    sym=[[0.5*(raw[i][j]+raw[j][i]) for j in range(6)] for i in range(6)]
    asym=max(abs(raw[i][j]-raw[j][i]) for i in range(6) for j in range(6))
    warnings=[]
    if max_force>max(1e-3,10*config['ftol']):
        warnings.append(f'maximum fmax={max_force:g} eV/A: verify minimization convergence')
    if asym>5.0:
        warnings.append(f'raw elastic tensor asymmetry={asym:g} GPa: check strain magnitude/convergence')

    result={'phase':phase,'strain':eps,'natoms':int(ref['N']),
            'a_start_A':config['a'],'equilibrated_volume_A3':ref['V'],
            'equilibrated_energy_eV':ref['E']}
    if phase in ('BCC','FCC'):
        c11=sum(sym[i][i] for i in range(3))/3
        c12=(sym[0][1]+sym[0][2]+sym[1][2])/3
        c44=sum(sym[i][i] for i in range(3,6))/3
        bulk=(c11+2*c12)/3
        g1=c44
        g2=(c11-c12)/2
        gvoigt=(c11-c12+3*c44)/5
        nu100=c12/(c11+c12) if abs(c11+c12)>1e-12 else None
        nuvoigt=((3*bulk-2*gvoigt)/(2*(3*bulk+gvoigt))
                 if abs(3*bulk+gvoigt)>1e-12 else None)
        stable=bool(c11-c12>0 and c11+2*c12>0 and c44>0)
        if not stable:
            warnings.append(f'{phase} violates cubic Born criteria')
        result.update({'C11_GPa':c11,'C12_GPa':c12,'C44_GPa':c44,
            'bulk_modulus_GPa':bulk,'shear_modulus_C44_GPa':g1,
            'shear_modulus_tetragonal_GPa':g2,'shear_modulus_Voigt_GPa':gvoigt,
            'poisson_ratio_100':nu100,'poisson_ratio_Voigt':nuvoigt,
            'stable_Born_cubic':stable})
        summary_lines=[f'C11={c11:.8f}',f'C12={c12:.8f}',f'C44={c44:.8f}',f'B={bulk:.8f}']
    else:


        c11=(sym[0][0]+sym[1][1])/2
        c12=sym[0][1]
        c13=(sym[0][2]+sym[1][2])/2
        c33=sym[2][2]
        c44=(sym[3][3]+sym[4][4])/2
        c66=sym[5][5]
        c66_theory=(c11-c12)/2
        c66_diff=c66-c66_theory
        bulk_voigt=(2*(c11+c12)+4*c13+c33)/9
        shear_voigt=(2*c11+c33-c12-2*c13+3*(2*c44+c66))/15
        stable=bool(c44>0 and c11>abs(c12) and (c11+c12)*c33>2*c13*c13)
        if not stable:
            warnings.append('HCP violates hexagonal Born mechanical stability criteria')

        ideal=[[0.]*6 for _ in range(6)]
        ideal[0][0]=ideal[1][1]=c11
        ideal[0][1]=ideal[1][0]=c12
        ideal[0][2]=ideal[2][0]=ideal[1][2]=ideal[2][1]=c13
        ideal[2][2]=c33
        ideal[3][3]=ideal[4][4]=c44
        ideal[5][5]=c66_theory
        deviation=max(abs(sym[i][j]-ideal[i][j]) for i in range(6) for j in range(6))
        if deviation>5:
            warnings.append(f'HCP tensor differs from hexagonal symmetry by {deviation:g} GPa')
        if abs(c66_diff)>5:
            warnings.append(f'C66 vs (C11-C12)/2 differs by {c66_diff:g} GPa')
        result.update({'c_start_A':config['c'],
            'C11_GPa':c11,'C12_GPa':c12,'C13_GPa':c13,'C33_GPa':c33,
            'C44_GPa':c44,'C66_GPa':c66,'C66_expected_GPa':c66_theory,
            'C66_minus_expected_GPa':c66_diff,
            'bulk_modulus_Voigt_GPa':bulk_voigt,'shear_modulus_Voigt_GPa':shear_voigt,
            'hexagonal_symmetry_deviation_GPa':deviation,
            'stable_Born_hexagonal':stable})
        summary_lines=[f'{name}={result[name+"_GPa"]:.8f}'
                       for name in ('C11','C12','C13','C33','C44','C66')]
        summary_lines += [f'C66_from_(C11-C12)/2={c66_theory:.8f}',
                          f'C66_deviation={c66_diff:.8f}',
                          f'Bulk_Voigt={bulk_voigt:.8f}']
    result.update({'max_force_eV_A':max_force,'raw_asymmetry_max_GPa':asym,
                   'warnings':warnings,'C_raw_GPa':raw,'C_symmetric_GPa':sym})
    (root/'elastic_result.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    csv_write(root/'elastic_result.csv', [{k:v for k,v in result.items()
                                        if not isinstance(v,(dict,list))}])
    csv_write(root/'elastic_tensor.csv',
              [{'component':f'C{i+1}{j+1}','raw_GPa':raw[i][j],
                'symmetric_GPa':sym[i][j]} for i in range(6) for j in range(6)])
    csv_write(root/'elastic_signed_stress.csv',signed)
    header=['Elastic tensor (symmetric GPa, Voigt xx yy zz yz xz xy)']
    header += [' '.join(f'{v:13.6f}' for v in row) for row in sym]
    header += ['', *summary_lines,f'Warnings: {warnings}']
    (root/'ElasticConstants.txt').write_text('\n'.join(header)+'\n')
    write_relaxed_cfg(root/'equilibrated.data',root/'relaxed.cfg.0',ref['E'],phase)
    print(f'[OK] {phase} elastic constants; '+', '.join(summary_lines))
    for warning in warnings:
        print('[WARN]',warning)
    return result

if __name__ == '__main__':
    collect()
