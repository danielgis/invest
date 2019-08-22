import unittest
import shutil
import tempfile
import os


class CLIGUITests(unittest.TestCase):
    def setUp(self):
        """Use a temporary workspace for all tests in this class."""
        self.workspace_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Remove the temporary workspace after a test run."""
        shutil.rmtree(self.workspace_dir)

    def test_run_model(self):
        """CLI-GUI: Run a model GUI through the cli."""
        from natcap.invest import cli
        from natcap.invest import delineateit
        parameter_set_path = os.path.join(
            os.path.dirname(__file__), '..', 'data', 'invest-sample-data',
            'delineateit_willamette.invs.json')

        # I tried patching the model import via mock, but GUI would hang.  I'd
        # rather have a reliable test that takes a few more seconds than a test
        # that hangs.
        cli.main([
            'quickrun',
            'delineateit',  # uses an exact modelname
            parameter_set_path,
            '--workspace', self.workspace_dir,
        ])
