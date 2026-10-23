from pathlib import Path
import os
import json

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = Path(os.environ.get("ZRV_CONFIG", PROJECT_ROOT / "config" / "paths.json"))
if not CONFIG_PATH.is_file():
    raise RuntimeError(f"Missing configuration file: {CONFIG_PATH}")
with CONFIG_PATH.open(encoding="utf-8") as stream:
    _CONFIG = json.load(stream)


def _value(name):
    return _CONFIG[name]


def _path(name):
    return Path(str(_value(name)))


POT_ROOT = Path(os.environ.get("ZRPOT_ROOT", str(_path("potential_root"))))
OUTPUT_ROOT = Path(os.environ.get("ZROUT_ROOT", str(_path("output_root"))))
LMP_EXE = str(_value("lammps_executable"))
SLURM_ACCOUNT = str(_value("slurm_account"))
TABLE_POTENTIAL = str(_value("nlh_table"))
TABLE_NAME = str(_CONFIG.get("nlh_table_name", "NLH"))
BATCH_CORES = int(_CONFIG.get("batch_cores", 192))
BATCH_TIME = str(_CONFIG.get("batch_time", "0-05:59"))

QSD_LATTICE_MODE = "legacy_hcp"
QSD_A = 3.232
QSD_CA = 1.6
QSD_PARALLEL = 8
QSD_RELAX_OTHERS = False
QSD_D_STEP = 0.1
QSD_CONFIG = {
    "QSD-A":  {"dmax": 5.0, "drag": 51, "target": None, "legacy_txt": "qsdA_repulsive_E.txt"},
    "QSD-AC": {"dmax": 3.0, "drag": 57, "target": 60,   "legacy_txt": "qsdAC_repulsive_E.txt"},
    "QSD-C":  {"dmax": 5.0, "drag": 51, "target": None, "legacy_txt": "qsdC_repulsive_E.txt"},
}
QSD_TESTS = set(QSD_CONFIG)


EOS_TESTS = {"EOS-BCC", "EOS-HCP"}
EOS_A_START = "2.00"
EOS_A_END = "5.00"
EOS_A_STEP = "0.10"
EOS_HCP_CA = "1.633"
EOS_REPEAT = (1, 1, 1)
EOS_RELAX_ATOMS = {"BCC": True, "HCP": False}
EOS_CORES = 8


COHESIVE_CONFIG = {
    "BCC": {"a": 3.62, "repeat": (3,3,3),
            "filename": "bcc_3x3x3.data", "primitive": "bcc_primitive.data"},
    "FCC": {"a": 4.55, "repeat": (2,2,2),
            "filename": "fcc_2x2x2.data", "primitive": "fcc_primitive.data"},
    "HCP": {"a": 3.233, "c": 5.147, "repeat": (3,3,2),
            "filename": "bulk_hcp_3x3x2.data", "primitive": "bulk_hcp_primitive.data"},
}
COHESIVE_TESTS = {"Cohesive-BCC", "Cohesive-FCC", "Cohesive-HCP"}
COHESIVE_CORES = 2


ELASTIC_TESTS = {"Elastic-BCC", "Elastic-FCC", "Elastic-HCP"}
ELASTIC_BCC_A = 3.583
ELASTIC_BCC_REPEAT = (1, 1, 1)
ELASTIC_BCC_STRAIN = 1.0e-4
ELASTIC_BCC_BOX_RELAX = "iso"
ELASTIC_BCC_ETOL = 0.0
ELASTIC_BCC_FTOL = 1.0e-8
ELASTIC_BCC_MAXITER = 10000
ELASTIC_BCC_MAXEVAL = 100000
ELASTIC_BCC_DMAX = 0.01


ELASTIC_FCC_A = 4.530
ELASTIC_FCC_REPEAT = (1, 1, 1)
ELASTIC_FCC_STRAIN = 1.0e-4
ELASTIC_FCC_BOX_RELAX = "iso"
ELASTIC_FCC_ETOL = 0.0
ELASTIC_FCC_FTOL = 1.0e-10
ELASTIC_FCC_MAXITER = 1000
ELASTIC_FCC_MAXEVAL = 10000
ELASTIC_FCC_DMAX = 0.01

ELASTIC_HCP_A = 3.230
ELASTIC_HCP_C = 5.180
ELASTIC_HCP_REPEAT = (1, 1, 1)
ELASTIC_HCP_STRAIN = 1.0e-4


ELASTIC_HCP_BOX_RELAX = "coupled_xy"
ELASTIC_HCP_ETOL = 0.0
ELASTIC_HCP_FTOL = 1.0e-10
ELASTIC_HCP_MAXITER = 1000
ELASTIC_HCP_MAXEVAL = 10000
ELASTIC_HCP_DMAX = 0.01
ELASTIC_CORES = 1


FORMATION_CONFIG = {
    'Formation-1_vac': {
        'kind': 'vacancy', 'a': 3.230, 'c': 5.147, 'repeat': (3,3,1),
        'relaxation': 'aniso', 'maxiter': 100000, 'maxeval': 100000,
    },
    'Formation-BO': {
        'kind': 'bo', 'a': 3.230, 'c': 5.147, 'repeat': (3,3,2),
        'relaxation': 'fixed', 'maxiter': 5000, 'maxeval': 10000,
        'bo_mode': 'octa',
        'dumbbell_disp': 0.55, 'min_bo_separation': 1.0,
    },
    'Formation-di_vac_inplane': {
        'kind': 'basal', 'a': 3.230, 'c': 5.147, 'repeat': (3,3,2),
        'relaxation': 'aniso', 'maxiter': 100000, 'maxeval': 100000,
        'nn_cut': 3.6,
    },
    'Formation-di_vac_outplane': {
        'kind': 'out_of_plane', 'a': 3.230, 'c': 5.147, 'repeat': (2,2,2),
        'relaxation': 'aniso', 'maxiter': 100000, 'maxeval': 100000,
        'nn_cut': 3.6,
    },
}
FORMATION_TESTS = set(FORMATION_CONFIG)
FORMATION_CORES = {name:(3 if c['kind'] in ('basal','out_of_plane') else 2)
                   for name,c in FORMATION_CONFIG.items()}


MELT_CONFIG = {
    "Melt-BCC": {"lattice": "bcc", "a": 3.62, "nx": 5, "ny": 5, "nz": 40, "zmid": 20, "tag": "bcc"},
    "Melt-FCC": {"lattice": "fcc", "a": 4.55, "nx": 5, "ny": 5, "nz": 20, "zmid": 10, "tag": "fcc"},
    "Melt-HCP": {"lattice": "hcp", "a": 3.23, "nx": 5, "ny": 5, "nz": 20, "zmid": 10, "tag": "hcp"},
}
MELT_TESTS = set(MELT_CONFIG)
MELT_T_MIN = 1500
MELT_T_MAX = 2500
MELT_T_STEP = 5
MELT_T_SOLID = 1200.0
MELT_T_LIQUID = 2600.0
MELT_PEXT = 0.0
MELT_TIMESTEP = 0.001
MELT_PREP_SOLID_STEPS = 50_000
MELT_PREP_LIQUID_STEPS = 100_000
MELT_EQUIL_STEPS = 1_000_000
MELT_PROD_STEPS = 200_000
MELT_RDF_BINS = 200
MELT_RDF_CUTOFF = 8.0
MELT_COMM_CUTOFF = 11.0
MELT_PEAK_SEARCH_MAX = 4.5
MELT_SMOOTH_WINDOW = 21
MELT_LOW_T_WINDOW = (1500.0, 1650.0)
MELT_HIGH_T_WINDOW = (2350.0, 2500.0)
MELT_ARRAY_LIMIT = 4
MELT_WALLTIME = "0-23:59"


MELT_PREP_LIQUID_ENSEMBLE = "npt_z"



GSFE_SOURCE_PATHS = {k: Path(v) for k, v in _CONFIG["gsfe_configurations"].items()}
NEB_SOURCE_PATHS = {k: Path(v) for k, v in _CONFIG["neb_replicas"].items()}

PRISM_CONFIG_DIR = GSFE_SOURCE_PATHS["prism_I_1"]
PRISM_I2_CONFIG_DIR = GSFE_SOURCE_PATHS["prism_I_2"]
BASAL_CONFIG_DIR = GSFE_SOURCE_PATHS["basal_1"]
BASAL_I2_CONFIG_DIR = GSFE_SOURCE_PATHS["basal_2"]
PRISM_II1_CONFIG_DIR = GSFE_SOURCE_PATHS["prism_II_1"]
PRISM_II2_CONFIG_DIR = GSFE_SOURCE_PATHS["prism_II_2"]
PYRAM_I_1_A_CONFIG_DIR = GSFE_SOURCE_PATHS["pyram_I_1_a"]
PYRAM_I_1_AC_CONFIG_DIR = GSFE_SOURCE_PATHS["pyram_I_1_a+c"]
PYRAM_I_2_CONFIG_DIR = GSFE_SOURCE_PATHS["pyram_I_2"]
PYRAM_II_1_AC_CONFIG_DIR = GSFE_SOURCE_PATHS["pyram_II_1_a+c"]
PYRAM_II_2_CONFIG_DIR = GSFE_SOURCE_PATHS["pyram_II_2"]

SIA_BASAL_REPLICA_SRC = NEB_SOURCE_PATHS["SIA_basal"]
SIA_NONBASAL_REPLICA_SRC = NEB_SOURCE_PATHS["SIA_nonbasal"]
VAC_BASAL_REPLICA_SRC = NEB_SOURCE_PATHS["vacancy_basal"]
VAC_NONBASAL_REPLICA_SRC = NEB_SOURCE_PATHS["vacancy_nonbasal"]


NEB_NREP = 27
SIA_LMP_START_LINE = 17
SIA_LMP_END_LINE   = 161
VAC_LMP_START_LINE = 17
VAC_LMP_END_LINE   = 111


GSFE_TEST_CONFIG = {
    "prism_I_1": {
        "src": PRISM_CONFIG_DIR,
        "input_type": "ppp_prism",
    },
    "basal_1": {
        "src": BASAL_CONFIG_DIR,
        "input_type": "ppp_basal",
    },
    "prism_I_2": {
        "src": PRISM_I2_CONFIG_DIR,
        "input_type": "pps_sd",
    },
    "basal_2": {
        "src": BASAL_I2_CONFIG_DIR,
        "input_type": "pps_sd",
    },
    "prism_II_1": {
        "src": PRISM_II1_CONFIG_DIR,
        "input_type": "pps_sd",
    },
    "prism_II_2": {
        "src": PRISM_II2_CONFIG_DIR,
        "input_type": "pps_sd",
    },
    "pyram_I_1_a": {
        "src": PYRAM_I_1_A_CONFIG_DIR,
        "input_type": "pps_sd",
    },
    "pyram_I_1_a+c": {
        "src": PYRAM_I_1_AC_CONFIG_DIR,
        "input_type": "pps_sd",
    },
    "pyram_I_2": {
        "src": PYRAM_I_2_CONFIG_DIR,
        "input_type": "pps_sd",
    },
    "pyram_II_1_a+c": {
        "src": PYRAM_II_1_AC_CONFIG_DIR,
        "input_type": "pps_sd",
    },
    "pyram_II_2": {
        "src": PYRAM_II_2_CONFIG_DIR,
        "input_type": "pps_sd",
    },
}

GSFE_TESTS = set(GSFE_TEST_CONFIG.keys())
NEB_TESTS  = {"SIA_basal", "SIA_nonbasal", "vacancy_basal", "vacancy_nonbasal"}


NEB_TEST_CONFIG = {
    "SIA_basal": {
        "src":   SIA_BASAL_REPLICA_SRC,
        "input": "in_basal.neb",
        "start_line": SIA_LMP_START_LINE, "end_line": SIA_LMP_END_LINE,
        "title": "SIA basal migration",
    },
    "SIA_nonbasal": {
        "src":   SIA_NONBASAL_REPLICA_SRC,
        "input": "in_non_basal.neb",
        "start_line": SIA_LMP_START_LINE, "end_line": SIA_LMP_END_LINE,
        "title": "SIA non-basal migration",
    },
    "vacancy_basal": {
        "src":   VAC_BASAL_REPLICA_SRC,
        "input": "in.neb",
        "start_line": VAC_LMP_START_LINE, "end_line": VAC_LMP_END_LINE,
        "title": "vacancy basal migration",
    },
    "vacancy_nonbasal": {
        "src":   VAC_NONBASAL_REPLICA_SRC,
        "input": "in.neb",
        "start_line": VAC_LMP_START_LINE, "end_line": VAC_LMP_END_LINE,
        "title": "vacancy non-basal migration",
    },
}


