import os
import pathlib
import subprocess
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]

class TestCommandLine(unittest.TestCase):
    def call(self, *args, env=None):
        return subprocess.run([str(ROOT / "bin" / "zrvalidate"), *args],
                              text=True, capture_output=True, env=env)

    def test_categories(self):
        proc = self.call("list")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("formation", proc.stdout)
        self.assertIn("elastic", proc.stdout)

    def test_placeholder_validation(self):
        proc = self.call("check", "--only", "elastic", "--phase", "HCP")
        self.assertEqual(proc.returncode, 2)
        self.assertIn("potential_root", proc.stderr)

    def test_generate_selected_without_lammps(self):
        with tempfile.TemporaryDirectory() as td:
            tmp=pathlib.Path(td)
            pot=tmp/'potentials';pot.mkdir()
            (pot/'seed_0.almtp').write_text('fixture')
            output=tmp/'results'
            cfg_file=tmp/'paths.json'
            import json
            cfg=json.loads((ROOT/'config'/'paths.json').read_text())
            cfg.update(potential_root=str(pot),output_root=str(output),
                       lammps_executable='lammps_mock',nlh_table=str(tmp/'table'),
                       slurm_account='test-account')
            cfg_file.write_text(json.dumps(cfg))
            env=os.environ.copy();env['ZRV_CONFIG']=str(cfg_file)
            proc=self.call('generate','--only','elastic','--phase','HCP','--seed','0',env=env)
            self.assertEqual(proc.returncode,0,proc.stderr+'\n'+proc.stdout)
            dest=output/'seed_0'/'Elastic-HCP'
            self.assertTrue((dest/'in.equilibrate.lmp').is_file())
            self.assertEqual(len(list(dest.glob('in.strain_*.lmp'))),12)
            self.assertTrue((output/'submit_elastic.sh').exists())
