#!/usr/bin/env python3


from zr_validation.settings import *
from zr_validation.core.common import pair_block

def prism_in_gsfe(pot_path: str) -> str:
    return f"""units          metal
atom_style     atomic
boundary       p p p

read_data input.pos${{step_index}}

{pair_block(pot_path)}

mass         1 91.224

fix 1 all setforce 0 0 NULL

thermo 1
thermo_style custom step cpu press etotal lx cella cellb cellc cellalpha cellbeta cellgamma

min_style sd
minimize 1.0e-10 1.0e-10 100000 100000

variable area equal lx*ly
print "${{area}}" append area_energy.txt
"""

def basal_in_gsfe(pot_path: str) -> str:
    return f"""units          metal
atom_style     atomic
boundary       p p p

read_data input.pos${{step_index}}

{pair_block(pot_path)}

mass         1 91.224

fix 1 all setforce 0 0 NULL

thermo 1
thermo_style custom step cpu press etotal lx cella cellb cellc cellalpha cellbeta cellgamma

min_style cg
minimize 1.0e-10 1.0e-10 100000 100000

variable area equal lx*ly
print "${{area}}" append area_energy.txt
"""

def surface_pps_in_gsfe(pot_path: str) -> str:


    return f"""units          metal
atom_style     atomic
boundary       p p s

read_data input.pos${{step_index}}

{pair_block(pot_path)}

mass         1 91.224

fix 1 all setforce 0 0 NULL

thermo 1
thermo_style custom step cpu press etotal lx cella cellb cellc cellalpha cellbeta cellgamma

min_style sd
minimize 1.0e-10 1.0e-10 100000 100000

variable area equal lx*ly
print "${{area}}" append area_energy.txt
"""

def gsfe_input_for_test(test_name: str, pot_path: str) -> str:

    input_type = GSFE_TEST_CONFIG[test_name]["input_type"]

    if input_type == "ppp_prism":
        return prism_in_gsfe(pot_path)
    if input_type == "ppp_basal":
        return basal_in_gsfe(pot_path)
    if input_type == "pps_sd":
        return surface_pps_in_gsfe(pot_path)

    raise ValueError(f"Unknown GSFE input_type for {test_name}: {input_type}")

def single_auto_gsfe_run_sh(work_dir: str, test_name: str) -> str:


    return f"""#!/bin/bash
#SBATCH --account={SLURM_ACCOUNT}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={BATCH_CORES}
#SBATCH --cpus-per-task=1
#SBATCH --time={BATCH_TIME}

module purge
module load StdEnv/2023
module load gcc/12.3
module load openmpi/4.1.5

set -e

LMP={LMP_EXE}
parent_dir={work_dir}

mapfile -t CONFIGS < <(find -L "${{parent_dir}}/configuration" -maxdepth 1 -type f \
    -name 'input.pos*' -printf '%p\\n' | sort -V)

NCONFIG=${{#CONFIGS[@]}}

if (( NCONFIG == 0 )); then
    echo "ERROR: {test_name} found no input.pos* in ${{parent_dir}}/configuration" >&2
    exit 1
fi

if (( NCONFIG < {BATCH_CORES} )); then
    MAX_PARALLEL=$NCONFIG
else
    MAX_PARALLEL={BATCH_CORES}
fi

echo "{test_name}: detected $NCONFIG configurations"
echo "{test_name}: using up to $MAX_PARALLEL concurrent cores"

count=0
for cfg in "${{CONFIGS[@]}}"; do
    fname=$(basename "$cfg")
    step_index="${{fname#input.pos}}"

    step_dir="${{parent_dir}}/Step_${{step_index}}"
    mkdir -p "$step_dir"

    cp "${{parent_dir}}/in.gsfe" "$step_dir/"
    cp "$cfg" "$step_dir/input.pos${{step_index}}"
    sed -i "s/\\${{step_index}}/$step_index/g" "$step_dir/in.gsfe"

    (
        cd "$step_dir"
        srun --exclusive --ntasks=1 --cpus-per-task=1 \
             "$LMP" -in in.gsfe > lammps.out 2>&1
    ) &

    count=$((count + 1))

    if (( count % MAX_PARALLEL == 0 )); then
        wait
    fi
done

wait
echo "{test_name}: all $NCONFIG configurations completed."
"""

def auto_gsfe_collect_sh(test_name: str) -> str:


    return f"""#!/bin/bash
set -e

raw_file="./gsfe_raw.txt"
result_file="./gsfe_result.txt"
> "$raw_file"
> "$result_file"

mapfile -t step_dirs < <(
    find . -maxdepth 1 -type d -name 'Step_*' -printf '%f\\n' | sort -V
)

NCONFIG=${{#step_dirs[@]}}

if (( NCONFIG == 0 )); then
    echo "ERROR: {test_name} found no Step_* directories." >&2
    exit 1
fi

echo "{test_name}: detected $NCONFIG Step_* directories"

for folder_name in "${{step_dirs[@]}}"; do
    step="${{folder_name#Step_}}"
    logfile="${{folder_name}}/log.lammps"

    if [ ! -f "$logfile" ]; then
        echo "ERROR: ${{logfile}} not found" >&2
        exit 1
    fi

    energy=$(awk '/Energy initial, next-to-last, final =/{{getline; print $3}}' "$logfile")

    if [ -z "$energy" ]; then
        echo "ERROR: energy not found in ${{logfile}}" >&2
        exit 1
    fi

    echo "${{step}} ${{energy}}" >> "$raw_file"
done

first_dir="${{step_dirs[0]}}"
area_file="${{first_dir}}/area_energy.txt"

if [ ! -f "$area_file" ]; then
    echo "ERROR: ${{area_file}} not found, cannot compute GSFE." >&2
    exit 1
fi

area=$(awk 'NR==1{{print $1}}' "$area_file")
if [ -z "$area" ] || [ "$area" = "0" ]; then
    echo "ERROR: invalid area in ${{area_file}}" >&2
    exit 1
fi

E0=$(awk 'NR==1{{print $2}}' "$raw_file")
if [ -z "$E0" ]; then
    echo "ERROR: reference energy not found." >&2
    exit 1
fi

echo "Reference configuration = ${{first_dir}}"
echo "Reference area          = ${{area}} Angstrom^2"
echo "Reference energy        = ${{E0}} eV"

echo "# step  GSFE(meV/Angstrom^2)" > "$result_file"

while read -r step energy; do
    gsfe=$(awk -v e="$energy" -v e0="$E0" -v a="$area" \
               'BEGIN {{ printf "%.6f", (e - e0) / a * 1000 }}')
    echo "${{step}}  ${{gsfe}}" >> "$result_file"
done < "$raw_file"

echo "Done."
echo "  Configurations: $NCONFIG"
echo "  Raw energies  : $raw_file"
echo "  GSFE results  : $result_file"
"""

