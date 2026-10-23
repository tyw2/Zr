#!/usr/bin/env python3


from zr_validation.settings import *
from zr_validation.core.common import pair_block

def melt_temps():
    if MELT_T_STEP <= 0 or (MELT_T_MAX - MELT_T_MIN) % MELT_T_STEP:
        raise ValueError("MELT_T_MIN/MAX/STEP must produce an integer number of temperatures")
    return list(range(MELT_T_MIN, MELT_T_MAX + 1, MELT_T_STEP))

def melt_interface_name(test_name: str):
    return f"interface_{MELT_CONFIG[test_name]['tag']}_Zr.dat"

def melt_prep_input(test_name: str, pot_path: str) -> str:
    c = MELT_CONFIG[test_name]
    if MELT_PREP_LIQUID_ENSEMBLE == "npt_z":


        melt_fix = ("fix f2 g_liquid npt temp ${T_LIQUID} ${T_LIQUID} 0.1 "
                    "z ${Pext} ${Pext} 1.0 dilate g_liquid")
    elif MELT_PREP_LIQUID_ENSEMBLE == "nvt":
        melt_fix = "fix f2 g_liquid nvt temp ${T_LIQUID} ${T_LIQUID} 0.1"
    else:
        raise ValueError("MELT_PREP_LIQUID_ENSEMBLE must be 'npt_z' or 'nvt'")
    return f'''
clear
units           metal
dimension       3
boundary        p p p
atom_style      atomic

variable a0 equal {c['a']}
variable Pext equal {MELT_PEXT}
variable T_SOLID equal {MELT_T_SOLID}
variable T_LIQUID equal {MELT_T_LIQUID}
variable NX equal {c['nx']}
variable NY equal {c['ny']}
variable NZ equal {c['nz']}
variable ZMID equal {c['zmid']}

lattice {c['lattice']} ${{a0}}
region box block 0 ${{NX}} 0 ${{NY}} 0 ${{NZ}} units lattice
create_box 1 box
create_atoms 1 box
mass 1 91.224

{pair_block(pot_path)}

neighbor 2.0 bin
neigh_modify delay 0 every 1 check yes
timestep {MELT_TIMESTEP}

region r_solid block INF INF INF INF 0 ${{ZMID}} units lattice
group g_solid region r_solid
group g_liquid subtract all g_solid

thermo 1000
thermo_style custom step temp pe press pxx pyy pzz lx ly lz density
dump d0 all custom 1000 interface_{c['tag']}_prep.dump id type x y z
dump_modify d0 sort id

velocity all create ${{T_SOLID}} 12345 mom yes rot yes dist gaussian
fix f1 all npt temp ${{T_SOLID}} ${{T_SOLID}} 0.1 iso ${{Pext}} ${{Pext}} 1.0
run {MELT_PREP_SOLID_STEPS}
unfix f1

velocity g_solid set 0.0 0.0 0.0
fix freeze g_solid setforce 0.0 0.0 0.0
velocity g_liquid create ${{T_LIQUID}} 23456 mom yes rot yes dist gaussian
{melt_fix}
run {MELT_PREP_LIQUID_STEPS}
unfix f2
unfix freeze
undump d0
write_data {melt_interface_name(test_name)}
write_restart interface_{c['tag']}_Zr.restart
'''

def melt_scan_input(test_name: str, pot_path: str) -> str:

    if MELT_PROD_STEPS % 100:
        raise ValueError("MELT_PROD_STEPS must be divisible by 100")
    return f'''
clear
units           metal
dimension       3
boundary        p p p
atom_style      atomic
read_data       {melt_interface_name(test_name)}

{pair_block(pot_path)}

neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes
timestep        {MELT_TIMESTEP}
thermo          2000
thermo_style    custom step temp pe press lx ly lz density

velocity        all create ${{T}} 12345 mom yes rot yes dist gaussian
fix             f1 all npt temp ${{T}} ${{T}} 0.1 aniso {MELT_PEXT} {MELT_PEXT} 1.0
run             {MELT_EQUIL_STEPS}
write_data      ${{outdir}}/equilibrated.dat
write_restart   ${{outdir}}/equilibrated.restart

comm_modify     cutoff {MELT_COMM_CUTOFF}
compute         myRDF all rdf {MELT_RDF_BINS} 1 1 cutoff {MELT_RDF_CUTOFF}
fix             rdfout all ave/time 100 {MELT_PROD_STEPS // 100} {MELT_PROD_STEPS} c_myRDF[*] file ${{outdir}}/rdf.txt mode vector
run             {MELT_PROD_STEPS}
write_data      ${{outdir}}/coexist.dat
'''

