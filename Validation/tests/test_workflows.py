import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class WorkflowGenerationSmokeTest(unittest.TestCase):
    def test_all_eight_workflow_generators(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            potentials = base / "potentials"
            potentials.mkdir()
            (potentials / "seed_0.almtp").write_text("test-only")
            config = json.loads((ROOT / "config" / "paths.json").read_text())
            config.update(potential_root=str(potentials),
                          output_root=str(base / "results"),
                          lammps_executable="test_lammps",
                          nlh_table="test_nlh.table",
                          slurm_account="test-slurm-account")
            for name in config["gsfe_configurations"]:
                folder = base / "templates" / "gsfe" / name
                folder.mkdir(parents=True)
                (folder / "input.pos0").write_text("dummy")
                config["gsfe_configurations"][name] = str(folder)
            for name in config["neb_replicas"]:
                folder = base / "templates" / "neb" / name
                folder.mkdir(parents=True)
                for idx in range(27):
                    (folder / f"replica.{idx}.coords").write_text("2\n1 0 0 0\n2 1 1 1\n")
                (folder / "replica.0.lmp").write_text("dummy")
                (folder / "replica.26.lmp").write_text("dummy")
                config["neb_replicas"][name] = str(folder)
            conf_file = base / "paths.json"
            conf_file.write_text(json.dumps(config))
            env = dict(os.environ, ZRV_CONFIG=str(conf_file))
            for category in ("gsfe", "neb", "qsd", "melting", "eos", "cohesive", "elastic", "formation"):
                with self.subTest(category=category):
                    run = subprocess.run(
                        [str(ROOT / "bin" / "zrvalidate"), "generate", "--only", category, "--seed", "0"],
                        env=env, capture_output=True, text=True, timeout=60,
                    )
                    self.assertEqual(run.returncode, 0, run.stderr + "\n" + run.stdout)
            output = base / "results"
            seed = output / "seed_0"
            for file in (seed / "prism_I_1/in.gsfe",
                         seed / "SIA_basal/in_basal.neb",
                         seed / "QSD-A/qsd.in",
                         seed / "Melt-BCC/in_scan.lmp",
                         seed / "EOS-HCP/run.sh",
                         seed / "Cohesive-FCC/run.sh",
                         seed / "Elastic-FCC/in.strain_6_pos.lmp",
                         seed / "Formation-BO/execute_formation.sh"):
                self.assertTrue(file.is_file(), str(file))
            for file in output.rglob("*.sh"):
                result = subprocess.run(["bash", "-n", str(file)], capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, str(file) + result.stderr)
