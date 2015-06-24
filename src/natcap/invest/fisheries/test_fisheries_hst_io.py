import unittest
import os
import pprint

import numpy as np
from numpy import testing

import fisheries_hst_io as io
from fisheries_hst_io import MissingParameter

pp = pprint.PrettyPrinter(indent=4)

workspace_dir = '../../test/invest-data/test/data/fisheries/'
data_dir = '../../test/invest-data/test/data/fisheries/'
inputs_dir = os.path.join(data_dir, 'habitat_scenario_tool/inputs')
outputs_dir = os.path.join(data_dir, 'habitat_scenario_tool/outputs')

All_Parameters = ['Classes', 'Duration', 'Exploitationfraction', 'Fecundity',
                  'Larvaldispersal', 'Maturity', 'Regions', 'Surv_nat_xsa',
                  'Weight', 'Vulnfishing']
Necessary_Params = ['Classes', 'Exploitationfraction', 'Maturity', 'Regions',
                    'Surv_nat_xsa', 'Vulnfishing']


class TestPopulationParamsIO(unittest.TestCase):
    def setUp(self):
        self.args = {
            'workspace_dir': workspace_dir,
            'sexsp': 2,
            'population_csv_uri': os.path.join(inputs_dir, 'pop_params.csv'),
            'habitat_chg_csv_uri': os.path.join(
                inputs_dir, 'habitat_chg_params.csv'),
            'habitat_dep_csv_uri': os.path.join(
                inputs_dir, 'habitat_dep_params.csv'),
            'gamma': 0.5,
        }

        self.check = {
            # User Args
            # 'workspace_dir': workspace_dir,
            # 'sexsp': 2,
            # 'population_csv_uri': os.path.join(inputs_dir, 'pop_params.csv'),
            # 'habitat_csv_uri': os.path.join(inputs_dir, 'habitat_params.csv'),
            # 'gamma': 0.5,

            # Pop CSV Args
            'Surv_nat_xsa': np.array(
                [[[1.56916669e-06, 7.25000000e-01, 7.25000000e-01,
                   7.25000000e-01, 7.25000000e-01],
                 [1.56916669e-06, 7.25000000e-01, 7.25000000e-01,
                  7.25000000e-01, 5.27500000e-01]],

                 [[2.00204026e-06, 7.25000000e-01, 7.25000000e-01,
                   7.25000000e-01, 7.25000000e-01],
                  [2.00204026e-06, 7.25000000e-01, 7.25000000e-01,
                   7.25000000e-01, 5.27500000e-01]],

                 [[1.46094830e-06, 7.25000000e-01, 7.25000000e-01,
                   7.25000000e-01, 7.25000000e-01],
                  [1.46094830e-06, 7.25000000e-01, 7.25000000e-01,
                   7.25000000e-01, 5.27500000e-01]],

                 [[1.67738509e-06, 7.25000000e-01, 7.25000000e-01,
                   7.25000000e-01, 7.25000000e-01],
                  [1.67738509e-06, 7.25000000e-01, 7.25000000e-01,
                   7.25000000e-01, 5.27500000e-01]],

                 [[2.11025866e-06, 7.25000000e-01, 7.25000000e-01,
                   7.25000000e-01, 7.25000000e-01],
                  [2.11025866e-06, 7.25000000e-01, 7.25000000e-01,
                   7.25000000e-01, 5.27500000e-01]],

                 [[2.27258625e-06, 7.25000000e-01, 7.25000000e-01,
                   7.25000000e-01, 7.25000000e-01],
                  [2.27258625e-06, 7.25000000e-01, 7.25000000e-01,
                   7.25000000e-01, 5.27500000e-01]]]),

            'Classes': np.array(['Class0', 'Class1', 'Class2', 'Class3', 'Class4']),
            'Class_vectors': {
                'Vulnfishing': np.array([[0, 0, 0, 0, 0], [0, 0, 0, 0, 1]]),
                'Maturity': np.array([[0, 0, 0, 1, 1], [0, 0, 0, 0, 0]]),
            },
            'Regions': np.array(['Region1',
                                 'Region2',
                                 'Region3',
                                 'Region4',
                                 'Region5',
                                 'Region6']),
            'Region_vectors': {
                'Exploitationfraction': np.array(
                    [0.47, 0.47, 0.47, 0.47, 0.47, 0.47]),
                'Larvaldispersal': np.array(
                    [0.09, 0.12, 0.18, 0.29, 0.17, 0.15]),
            },
        }

    def test_parse_popu_params(self):
        uri = self.args['population_csv_uri']
        sexsp = self.args['sexsp']
        guess = io._parse_population_csv(uri, sexsp)
        # pp.pprint(guess['Surv_nat_xsa'])
        # pp.pprint(guess)

        for k in self.check.keys():
            if k not in [
                # User Args not added back to dictionary yet
                'workspace_dir',
                'sexsp',
                'population_csv_uri',
                'habitat_csv_uri',
                'gamma'
            ]:
                g = guess[k]
                c = self.check[k]
                if isinstance(g, np.ndarray):
                    np.testing.assert_array_almost_equal(g, c)
                elif isinstance(g, list):
                    assert(all(map(lambda x, y: x == y, g, c)))
                elif isinstance(g, dict):
                    for m in g.keys():
                        n = g[m]
                        d = c[m]
                        if isinstance(n, np.ndarray):
                            np.testing.assert_array_almost_equal(n, d)
                        elif isinstance(n, list):
                            assert(all(map(lambda x, y: x == y, n, d)))
                        else:
                            assert(n == d)
                else:
                    assert(g == c)

    def test_verify_popu_params(self):
        # Check for exception when Survival Matrix elements not between [0, 1]
        args_mod = self.args
        args_mod['population_csv_uri'] = os.path.join(
            inputs_dir, 'pop_params_bad_S_elements.csv')
        with self.assertRaises(ValueError):
            io.read_population_csv(args_mod)

        # Check for exception when Survival Matrix does not exist
        pass


class TestHabitatDepParamsIO(unittest.TestCase):
    def setUp(self):
        self.args = {
            'workspace_dir': workspace_dir,
            'sexsp': 'Yes',
            'population_csv_uri': os.path.join(inputs_dir, 'pop_params.csv'),
            'habitat_chg_csv_uri': os.path.join(
                inputs_dir, 'habitat_chg_params.csv'),
            'habitat_dep_csv_uri': os.path.join(
                inputs_dir, 'habitat_dep_params.csv'),

            'gamma': 0.5,
        }

        self.check = {
            # User Args
            'workspace_dir': workspace_dir,
            'sexsp': 'Yes',
            'population_csv_uri': os.path.join(inputs_dir, 'pop_params.csv'),
            'habitat_dep_csv_uri': os.path.join(inputs_dir, 'habitat_dep_params.csv'),
            'habitat_chg_csv_uri': os.path.join(inputs_dir, 'habitat_chg_params.csv'),
            'gamma': 0.5,

            # Habitat CSV Args
            'Habitats': ['Habitat1', 'Habitat2', 'Habitat3'],
            'Hab_classes': ['Class0', 'Class1', 'Class2', 'Class3', 'Class4'],
            'Hab_regions': [
                'Region1',
                'Region2',
                'Region3',
                'Region4',
                'Region5',
                'Region6'],
            'Hab_chg_hx': np.array(
                [[-0.25, -0.25, -0.25, -0.25, -0.25, -0.25],
                 [-0.25, -0.25, -0.25, -0.25, -0.25, -0.25],
                 [-0.25, -0.25, -0.25, -0.25, -0.25, -0.25]]),
            'Hab_dep_ha': np.array(
                [[0.0, 0.0, 1.0, 1.0, 1.0],
                 [1.0, 1.0, 0.0, 0.0, 0.0],
                 [1.0, 1.0, 0.0, 0.0, 0.0]]),
            'Hab_class_mvmt_a': np.array([]),
            'Hab_dep_num_a': np.array([]),
        }

    def test_parse_habitat_dep_params(self):
        # Demonstrate able to parse the tables correctly
        guess = io._parse_habitat_dep_csv(self.args)
        # pp.pprint(guess)
        testing.assert_array_almost_equal(
            guess['Hab_dep_ha'], self.check['Hab_dep_ha'])
        assert(all(map(
            lambda x, y: x == y, guess['Habitats'], self.check['Habitats'])))
        assert(all(map(
            lambda x, y: x == y,
            guess['Hab_classes'],
            self.check['Hab_classes'])))

    def test_read_habitat_dep_params(self):
        # Test that it works on correct CSV
        # guess = io.read_habitat_dep_csv(self.args)
        # pp.pprint(guess)

        # Test for exceptions
        args_mod = self.args
        args_mod['habitat_dep_csv_uri'] = os.path.join(
            inputs_dir, 'habitat_dep_params_bad_elements.csv')

        # Hab_dep_ha elements must be between [0,1]
        with self.assertRaises(ValueError):
            io.read_habitat_dep_csv(args_mod)


class TestHabitatChgParamsIO(unittest.TestCase):
    def setUp(self):
        self.args = {
            'workspace_dir': workspace_dir,
            'sexsp': 'Yes',
            'population_csv_uri': os.path.join(inputs_dir, 'pop_params.csv'),
            'habitat_chg_csv_uri': os.path.join(
                inputs_dir, 'habitat_chg_params.csv'),
            'habitat_dep_csv_uri': os.path.join(
                inputs_dir, 'habitat_dep_params.csv'),
            'gamma': 0.5,
        }

        self.check = {
            # User Args
            'workspace_dir': workspace_dir,
            'sexsp': 'Yes',
            'population_csv_uri': os.path.join(inputs_dir, 'pop_params.csv'),
            'habitat_dep_csv_uri': os.path.join(inputs_dir, 'habitat_dep_params.csv'),
            'habitat_chg_csv_uri': os.path.join(inputs_dir, 'habitat_chg_params.csv'),
            'gamma': 0.5,

            # Habitat CSV Args
            'Habitats': ['Habitat1', 'Habitat2', 'Habitat3'],
            'Hab_classes': ['Class0', 'Class1', 'Class2', 'Class3', 'Class4'],
            'Hab_regions': [
                'Region1',
                'Region2',
                'Region3',
                'Region4',
                'Region5',
                'Region6'],
            'Hab_chg_hx': np.array(
                [[-0.25, -0.25, -0.25, -0.25, -0.25, -0.25],
                 [-0.25, -0.25, -0.25, -0.25, -0.25, -0.25],
                 [-0.25, -0.25, -0.25, -0.25, -0.25, -0.25]]),
            'Hab_dep_ha': np.array(
                [[0.0, 0.0, 1.0, 1.0, 1.0],
                 [1.0, 1.0, 0.0, 0.0, 0.0],
                 [1.0, 1.0, 0.0, 0.0, 0.0]]),
            'Hab_class_mvmt_a': np.array([]),
            'Hab_dep_num_a': np.array([]),
        }

    def test_parse_habitat_chg_params(self):
        # Demonstrate able to parse the tables correctly
        guess = io._parse_habitat_chg_csv(self.args)
        # pp.pprint(guess)
        testing.assert_array_almost_equal(
            guess['Hab_chg_hx'], self.check['Hab_chg_hx'])
        assert(all(map(
            lambda x, y: x == y, guess['Habitats'], self.check['Habitats'])))
        assert(all(map(
            lambda x, y: x == y,
            guess['Hab_regions'],
            self.check['Hab_regions'])))

    def test_read_habitat_chg_params(self):
        # Test that it works on correct CSV
        # guess = io.read_habitat_chg_csv(self.args)
        # pp.pprint(guess)

        # Test for exceptions
        args_mod = self.args
        args_mod['habitat_chg_csv_uri'] = os.path.join(
            inputs_dir, 'habitat_chg_params_bad_elements.csv')

        # Hab_chg_hx elements must be between [-1.0, +inf)
        with self.assertRaises(ValueError):
            io.read_habitat_chg_csv(args_mod)


class TestFetchArgs(unittest.TestCase):
    def setUp(self):
        self.args = {
            'workspace_dir': workspace_dir,
            'sexsp': 'Yes',
            'population_csv_uri': os.path.join(inputs_dir, 'pop_params.csv'),
            'habitat_chg_csv_uri': os.path.join(
                inputs_dir, 'habitat_chg_params.csv'),
            'habitat_dep_csv_uri': os.path.join(
                inputs_dir, 'habitat_dep_params.csv'),
            'gamma': 0.5,
        }

        self.check = {
            # User Args
            'workspace_dir': workspace_dir,
            'sexsp': 2,
            'population_csv_uri': os.path.join(inputs_dir, 'pop_params.csv'),
            'habitat_chg_csv_uri': os.path.join(
                inputs_dir, 'habitat_chg_params.csv'),
            'habitat_dep_csv_uri': os.path.join(
                inputs_dir, 'habitat_dep_params.csv'),
            'gamma': 0.5,

            # Pop CSV Args
            'Surv_nat_xsa': np.array(
                [[[1.56916669e-06, 7.25000000e-01, 7.25000000e-01,
                   7.25000000e-01, 7.25000000e-01],
                 [1.56916669e-06, 7.25000000e-01, 7.25000000e-01,
                  7.25000000e-01, 5.27500000e-01]],

                 [[2.00204026e-06, 7.25000000e-01, 7.25000000e-01,
                   7.25000000e-01, 7.25000000e-01],
                  [2.00204026e-06, 7.25000000e-01, 7.25000000e-01,
                   7.25000000e-01, 5.27500000e-01]],

                 [[1.46094830e-06, 7.25000000e-01, 7.25000000e-01,
                   7.25000000e-01, 7.25000000e-01],
                  [1.46094830e-06, 7.25000000e-01, 7.25000000e-01,
                   7.25000000e-01, 5.27500000e-01]],

                 [[1.67738509e-06, 7.25000000e-01, 7.25000000e-01,
                   7.25000000e-01, 7.25000000e-01],
                  [1.67738509e-06, 7.25000000e-01, 7.25000000e-01,
                   7.25000000e-01, 5.27500000e-01]],

                 [[2.11025866e-06, 7.25000000e-01, 7.25000000e-01,
                   7.25000000e-01, 7.25000000e-01],
                  [2.11025866e-06, 7.25000000e-01, 7.25000000e-01,
                   7.25000000e-01, 5.27500000e-01]],

                 [[2.27258625e-06, 7.25000000e-01, 7.25000000e-01,
                   7.25000000e-01, 7.25000000e-01],
                  [2.27258625e-06, 7.25000000e-01, 7.25000000e-01,
                   7.25000000e-01, 5.27500000e-01]]]),
            'Classes': np.array(
                ['Class0', 'Class1', 'Class2', 'Class3', 'Class4']),
            'Class_vectors': {
                'Vulnfishing': np.array([[0, 0, 0, 0, 0], [0, 0, 0, 0, 1]]),
                'Maturity': np.array([[0, 0, 0, 1, 1], [0, 0, 0, 0, 0]]),
            },
            'Regions': np.array(['Region1',
                                 'Region2',
                                 'Region3',
                                 'Region4',
                                 'Region5',
                                 'Region6']),
            'Region_vectors': {
                'Exploitationfraction': np.array(
                    [0.47, 0.47, 0.47, 0.47, 0.47, 0.47]),
                'Larvaldispersal': np.array(
                    [0.09, 0.12, 0.18, 0.29, 0.17, 0.15]),
            },

            # Habitat CSV Vars
            'Habitats': ['Habitat1', 'Habitat2', 'Habitat3'],
            'Hab_classes': ['Class0', 'Class1', 'Class2', 'Class3', 'Class4'],
            'Hab_regions': [
                'Region1',
                'Region2',
                'Region3',
                'Region4',
                'Region5',
                'Region6'],
            'Hab_chg_hx': np.array(
                [[-0.25, -0.25, -0.25, -0.25, -0.25, -0.25],
                 [-0.25, -0.25, -0.25, -0.25, -0.25, -0.25],
                 [-0.25, -0.25, -0.25, -0.25, -0.25, -0.25]]),
            'Hab_dep_ha': np.array(
                [[0.0, 0.0, 1.0, 1.0, 1.0],
                 [1.0, 1.0, 0.0, 0.0, 0.0],
                 [1.0, 1.0, 0.0, 0.0, 0.0]]),
            'Hab_class_mvmt_a': np.array([0, 0, 1, 0, 0]),
            'Hab_dep_num_a': np.array([2, 2, 1, 1, 1]),
        }

    def test_fetch_verify(self):
        guess = io.fetch_args(self.args)
        # pp.pprint(guess)
        for k in self.check.keys():
            g = guess[k]
            c = self.check[k]
            if isinstance(g, np.ndarray):
                np.testing.assert_array_almost_equal(g, c)
            elif isinstance(g, list):
                assert(all(map(lambda x, y: x == y, g, c)))
            elif isinstance(g, dict):
                for m in g.keys():
                    n = g[m]
                    d = c[m]
                    if isinstance(n, np.ndarray):
                        np.testing.assert_array_almost_equal(n, d)
                    elif isinstance(n, list):
                        assert(all(map(lambda x, y: x == y, n, d)))
                    else:
                        assert(n == d)
            else:
                assert(g == c)


class TestSavePopCSV(unittest.TestCase):
    def setUp(self):
        self.args = {
            'workspace_dir': '../../test/invest-data/test/data/fisheries/habitat_scenario_tool/',
            'sexsp': 'Yes',
            'population_csv_uri': '../../test/invest-data/test/data/fisheries/habitat_scenario_tool/output/pop_params_modified.csv',
            'habitat_chg_csv_uri': '../../test/invest-data/test/data/fisheries/habitat_scenario_tool/inputs/habitat_chg_params.csv',
            'habitat_dep_csv_uri': '../../test/invest-data/test/data/fisheries/habitat_scenario_tool/inputs/habitat_dep_params.csv',
            'gamma': 0.5,
        }
        self.vars_check = {
            'Class_vector_names': ['Maturity', 'Vulnfishing'],
            'Class_vectors': {
                'Maturity': np.array(
                  [[ 0.,  0.,  0.,  1.,  1.],
                   [ 0.,  0.,  0.,  0.,  0.]]),
                'Vulnfishing': np.array(
                    [[ 0.,  0.,  0.,  0.,  0.],
                     [ 0.,  0.,  0.,  0.,  1.]])},
            'Classes': ['Class0', 'Class1', 'Class2', 'Class3', 'Class4'],
            'Hab_chg_hx': np.array(
                [[-0.25, -0.25, -0.25, -0.25, -0.25, -0.25],
                 [-0.25, -0.25, -0.25, -0.25, -0.25, -0.25],
                 [-0.25, -0.25, -0.25, -0.25, -0.25, -0.25]]),
            'Hab_class_mvmt_a': np.array([0, 0, 1, 0, 0]),
            'Hab_classes': ['Class0', 'Class1', 'Class2', 'Class3', 'Class4'],
            'Hab_dep_ha': np.array([[ 0.,  0.,  1.,  1.,  1.],
               [ 1.,  1.,  0.,  0.,  0.],
               [ 1.,  1.,  0.,  0.,  0.]]),
            'Hab_dep_num_a': np.array([2, 2, 1, 1, 1]),
            'Hab_regions': [   'Region1',
                               'Region2',
                               'Region3',
                               'Region4',
                               'Region5',
                               'Region6'],
            'Habitats': ['Habitat1', 'Habitat2', 'Habitat3'],
            'Region_vectors': {
                'Exploitationfraction': np.array([ 0.47,  0.47,  0.47,  0.47,  0.47,  0.47]),
                'Larvaldispersal': np.array([ 0.09,  0.12,  0.18,  0.29,  0.17,  0.15])
                },
            'Region_vector_names': ['Exploitationfraction', 'Larvaldispersal'],
            'Regions': [   'Region1',
                           'Region2',
                           'Region3',
                           'Region4',
                           'Region5',
                           'Region6'],
            'Surv_nat_xsa': np.array(
                [[[  1.56916669e-06,   7.25000000e-01,   7.25000000e-01,
                     7.25000000e-01,   7.25000000e-01],
                  [  1.56916669e-06,   7.25000000e-01,   7.25000000e-01,
                     7.25000000e-01,   5.27500000e-01]],

               [[  2.00204026e-06,   7.25000000e-01,   7.25000000e-01,
                   7.25000000e-01,   7.25000000e-01],
                [  2.00204026e-06,   7.25000000e-01,   7.25000000e-01,
                   7.25000000e-01,   5.27500000e-01]],

               [[  1.46094830e-06,   7.25000000e-01,   7.25000000e-01,
                   7.25000000e-01,   7.25000000e-01],
                [  1.46094830e-06,   7.25000000e-01,   7.25000000e-01,
                   7.25000000e-01,   5.27500000e-01]],

               [[  1.67738509e-06,   7.25000000e-01,   7.25000000e-01,
                   7.25000000e-01,   7.25000000e-01],
                [  1.67738509e-06,   7.25000000e-01,   7.25000000e-01,
                   7.25000000e-01,   5.27500000e-01]],

               [[  2.11025866e-06,   7.25000000e-01,   7.25000000e-01,
                   7.25000000e-01,   7.25000000e-01],
                [  2.11025866e-06,   7.25000000e-01,   7.25000000e-01,
                   7.25000000e-01,   5.27500000e-01]],

               [[  2.27258625e-06,   7.25000000e-01,   7.25000000e-01,
                   7.25000000e-01,   7.25000000e-01],
                [  2.27258625e-06,   7.25000000e-01,   7.25000000e-01,
                   7.25000000e-01,   5.27500000e-01]]]),
            'gamma': 0.5,
            'habitat_chg_csv_uri': '../../test/invest-data/test/data/fisheries/habitat_scenario_tool/inputs/spreadsheet/habitat_chg_params.csv',
            'habitat_dep_csv_uri': '../../test/invest-data/test/data/fisheries/habitat_scenario_tool/inputs/spreadsheet/habitat_dep_params.csv',
            'output_dir': '../../test/invest-data/test/data/fisheries/habitat_scenario_tool/output',
            'population_csv_uri': '../../test/invest-data/test/data/fisheries/habitat_scenario_tool/inputs/spreadsheet/pop_params.csv',
            'sexsp': 2,
            'workspace_dir': '../../test/invest-data/test/data/fisheries/habitat_scenario_tool/',

            # Generated Variable
            'Surv_nat_xsa_mod': np.array(
                [[[  1.56916669e-06,   7.25000000e-01,   7.25000000e-01,
                   7.25000000e-01,   7.25000000e-01],
                [  1.56916669e-06,   7.25000000e-01,   7.25000000e-01,
                   7.25000000e-01,   5.27500000e-01]],

               [[  2.00204026e-06,   7.25000000e-01,   7.25000000e-01,
                   7.25000000e-01,   7.25000000e-01],
                [  2.00204026e-06,   7.25000000e-01,   7.25000000e-01,
                   7.25000000e-01,   5.27500000e-01]],

               [[  1.46094830e-06,   7.25000000e-01,   7.25000000e-01,
                   7.25000000e-01,   7.25000000e-01],
                [  1.46094830e-06,   7.25000000e-01,   7.25000000e-01,
                   7.25000000e-01,   5.27500000e-01]],

               [[  1.67738509e-06,   7.25000000e-01,   7.25000000e-01,
                   7.25000000e-01,   7.25000000e-01],
                [  1.67738509e-06,   7.25000000e-01,   7.25000000e-01,
                   7.25000000e-01,   5.27500000e-01]],

               [[  2.11025866e-06,   7.25000000e-01,   7.25000000e-01,
                   7.25000000e-01,   7.25000000e-01],
                [  2.11025866e-06,   7.25000000e-01,   7.25000000e-01,
                   7.25000000e-01,   5.27500000e-01]],

               [[  2.27258625e-06,   7.25000000e-01,   7.25000000e-01,
                   7.25000000e-01,   7.25000000e-01],
                [  2.27258625e-06,   7.25000000e-01,   7.25000000e-01,
                   7.25000000e-01,   5.27500000e-01]]]),
        }

    def test_save_population_csv(self):
        io.save_population_csv(self.vars_check)
        guess = io.fetch_args(self.args)
        # pp.pprint(self.vars_check)
        # pp.pprint(guess)

        for k in guess.keys():
            if k not in [
                # User Args not added back to dictionary yet
                'workspace_dir',
                'sexsp',
                'population_csv_uri',
                'habitat_chg_csv_uri',
                'habitat_dep_csv_uri',
                'gamma'
            ]:
                g = guess[k]
                c = self.vars_check[k]
                if isinstance(g, np.ndarray):
                    np.testing.assert_array_almost_equal(g, c)
                elif isinstance(g, list):
                    assert(all(map(lambda x, y: x == y, g, c)))
                elif isinstance(g, dict):
                    for m in g.keys():
                        n = g[m]
                        d = c[m]
                        if isinstance(n, np.ndarray):
                            np.testing.assert_array_almost_equal(n, d)
                        elif isinstance(n, list):
                            assert(all(map(lambda x, y: x == y, n, d)))
                        else:
                            assert(n == d)
                else:
                    assert(g == c)


if __name__ == '__main__':
    unittest.main()
