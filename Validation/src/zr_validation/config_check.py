from zr_validation import settings


def missing_fields(selected_tests=None):
    missing = []
    values = {
        "potential_root": str(settings.POT_ROOT),
        "output_root": str(settings.OUTPUT_ROOT),
        "lammps_executable": settings.LMP_EXE,
        "nlh_table": settings.TABLE_POTENTIAL,
        "slurm_account": settings.SLURM_ACCOUNT,
    }
    for name, value in values.items():
        if not value.strip() or value.strip() == "***":
            missing.append(name)
    if selected_tests:
        for name in sorted(selected_tests & settings.GSFE_TESTS):
            if str(settings.GSFE_SOURCE_PATHS[name]).strip() == "***":
                missing.append(f"gsfe_configurations.{name}")
        for name in sorted(selected_tests & settings.NEB_TESTS):
            if str(settings.NEB_SOURCE_PATHS[name]).strip() == "***":
                missing.append(f"neb_replicas.{name}")
    return missing


def ensure_configured(selected_tests=None):
    missing = missing_fields(selected_tests)
    if missing:
        raise RuntimeError(
            "Fill in config/paths.json (replace ***): " + ", ".join(missing)
        )
