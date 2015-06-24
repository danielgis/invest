import unittest
import pprint
import os

import numpy as np
from numpy import testing

import fisheries_hst as main
import fisheries_hst_io as io

pp = pprint.PrettyPrinter(indent=4)

workspace_dir = '../../test/invest-data/test/data/fisheries/'
data_dir = '../../test/invest-data/Fisheries'
inputs_dir = os.path.join(data_dir, 'input/Habitat_Scenario_Tool')
outputs_dir = os.path.join(workspace_dir, 'output')


class TestConvertSurvivalMatrix(unittest.TestCase):
    def setUp(self):
        self.args = {
            'workspace_dir': workspace_dir,
            'sexsp': 'No',
            'population_csv_uri': os.path.join(
                inputs_dir, 'pop_params.csv'),
            'habitat_dep_csv_uri': os.path.join(
                inputs_dir, 'habitat_dep_params.csv'),
            'habitat_chg_csv_uri': os.path.join(
                inputs_dir, 'habitat_chg_params.csv'),
            'gamma': 0.5,
        }
        self.check = {
            'workspace_dir': workspace_dir,
            'sexsp': 'No',
            'population_csv_uri': os.path.join(
                outputs_dir, 'pop_params_spreadsheet_mod.csv'),
            'habitat_dep_csv_uri': os.path.join(
                inputs_dir, 'habitat_dep_params.csv'),
            'habitat_chg_csv_uri': os.path.join(
                inputs_dir, 'habitat_chg_params.csv'),
            'gamma': 0.5,
        }

    def test_convert_spreadsheet(self):
        '''
        Test an example from the provided spreadsheet
        '''
        # Fetch pre and post variables
        vars_dict = io.fetch_args(self.args)
        check = io.fetch_args(self.check)

        # Run operation
        guess = main.convert_survival_matrix(vars_dict)

        # Check for correctness
        testing.assert_array_almost_equal(
            guess['Surv_nat_xsa_mod'], check['Surv_nat_xsa'])


if __name__ == '__main__':
    unittest.main()
