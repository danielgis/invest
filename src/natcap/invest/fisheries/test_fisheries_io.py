import unittest
import os
import pprint

from numpy import testing
import numpy as np

import fisheries_io

from fisheries_io import MissingParameter

data_directory = '../../test/invest-data/test/data/fisheries'
pp = pprint.PrettyPrinter(indent=4)


All_Parameters = ['Classes', 'Duration', 'Exploitationfraction', 'Fecundity',
                  'Larvaldispersal', 'Maturity', 'Regions', 'Survnaturalfrac',
                  'Weight', 'Vulnfishing']
Necessary_Params = ['Classes', 'Exploitationfraction', 'Maturity', 'Regions',
                    'Survnaturalfrac', 'Vulnfishing']


class TestPopulationParamsIO(unittest.TestCase):
    def test_parse_popu_params_blue_crab(self):
        uri = os.path.join(data_directory, 'CSVs/New_Params.csv')
        sexsp = 1
        pop_dict = fisheries_io._parse_population_csv(uri, sexsp)
        # Check that keys are correct
        Matching_Keys = [i for i in pop_dict.keys() if i in Necessary_Params]
        self.assertEqual(len(Matching_Keys), len(Necessary_Params))
        # Check that sexsp handled correctly
        self.assertEqual(len(pop_dict['Survnaturalfrac'][0]), sexsp)
        # Check that Class attribute lengths match
        self.assertEqual(
            len(pop_dict['Vulnfishing']), len(pop_dict['Maturity']))
        # Print Dictionary if debugging
        #pp.pprint(pop_dict)

    def test_parse_popu_params_sn(self):
        uri = os.path.join(data_directory, 'CSVs/TestCSV_SN_Syntax.csv')
        sexsp = 1
        pop_dict = fisheries_io._parse_population_csv(uri, sexsp)
        # Check that keys are correct
        Matching_Keys = [i for i in pop_dict.keys() if i in All_Parameters]
        self.assertEqual(len(Matching_Keys), len(All_Parameters))
        # Check that sexsp handled correctly
        self.assertEqual(len(pop_dict['Survnaturalfrac'][0]), sexsp)
        # Check that Class attribute lengths match
        self.assertEqual(
            len(pop_dict['Vulnfishing']), len(pop_dict['Maturity']))
        # Print Dictionary if debugging
        #pp.pprint(pop_dict)

    def test_parse_popu_params_ss(self):
        uri = os.path.join(data_directory, 'CSVs/TestCSV_SS_Syntax.csv')
        sexsp = 2
        pop_dict = fisheries_io._parse_population_csv(uri, sexsp)
        # Check that keys are correct
        Matching_Params = [i for i in pop_dict.keys() if i in All_Parameters]
        self.assertEqual(len(Matching_Params), len(All_Parameters))
        # Check that sexsp handled correctly
        self.assertEqual(len(pop_dict['Survnaturalfrac'][0]), sexsp)
        # Check that Class attribute lengths match
        self.assertEqual(
            len(pop_dict['Vulnfishing']), len(pop_dict['Maturity']))
        # Print Dictionary if debugging
        #pp.pprint(pop_dict)

    def test_read_popu_params(self):
        # Check that throws error when necessary information does not exist

        # Test with not all necessary params
        population_csv_uri = os.path.join(data_directory, 'CSVs/Fail/TestCSV_SN_Syntax_fail1.csv')
        args = {'population_csv_uri': population_csv_uri, 'sexsp': 1}
        with self.assertRaises(MissingParameter):
            fisheries_io.read_population_csv(args, population_csv_uri)

        # Test Stage-based without Duration vector
        population_csv_uri = os.path.join(data_directory, 'CSVs/Fail/TestCSV_SN_Syntax_fail2.csv')
        args['population_csv_uri'] = population_csv_uri
        args['recruitment_type'] = 'Beverton-Holt'
        args['population_type'] = 'Stage-Based'
        with self.assertRaises(MissingParameter):
            fisheries_io.read_population_csv(args, population_csv_uri)

        # Test B-H / Weight without Weight vector
        population_csv_uri = os.path.join(data_directory, 'CSVs/Fail/TestCSV_SN_Syntax_fail3.csv')
        args['population_csv_uri'] = population_csv_uri
        args['spawn_units'] = 'Weight'
        with self.assertRaises(MissingParameter):
            fisheries_io.read_population_csv(args, population_csv_uri)

        # Test Fecundity without Fecundity vector
        population_csv_uri = os.path.join(data_directory, 'CSVs/Fail/TestCSV_SN_Syntax_fail3.csv')
        args['population_csv_uri'] = population_csv_uri
        args['recruitment_type'] = 'Fecundity'
        args['harvest_units'] = 'Weight'
        with self.assertRaises(MissingParameter):
            fisheries_io.read_population_csv(args, population_csv_uri)

        '''
        # Check that throws error when incorrect information exists
        population_csv_uri = os.path.join(data_directory, 'CSVs/Fail/TestCSV_SN_Semantics_fail1.csv')
        args = {'population_csv_uri': population_csv_uri, 'sexsp': 1}
        self.assertRaises(
            MissingParameter, fisheries_io.read_population_csv(args))

        population_csv_uri = os.path.join(data_directory, 'CSVs/Fail/TestCSV_SN_Semantics_fail2.csv')
        args = {'population_csv_uri': population_csv_uri, 'sexsp': 1}
        self.assertRaises(
            MissingParameter, fisheries_io.read_population_csv(args))

        population_csv_uri = os.path.join(data_directory, 'CSVs/Fail/TestCSV_SN_Semantics_fail3.csv')
        args = {'population_csv_uri': population_csv_uri, 'sexsp': 1}
        self.assertRaises(
            MissingParameter, fisheries_io.read_population_csv(args))
        '''


class TestMigrationIO(unittest.TestCase):
    def test_parse_migration(self):
        uri = os.path.join(data_directory, 'migration/')
        args = {
            'migr_cont': True,
            'migration_dir': uri
        }
        class_list = ['larva', 'adult']
        mig_dict = fisheries_io._parse_migration_tables(args, class_list)
        #pp.pprint(mig_dict)
        self.assertIsInstance(mig_dict['adult'], np.matrix)
        self.assertEqual(
            mig_dict['adult'].shape[0], mig_dict['adult'].shape[1])

    def test_read_migration(self):
        uri = os.path.join(data_directory, 'migration/')
        args = {
            "migration_dir": uri,
            "migr_cont": True,
            }
        class_list = ['larva', 'other', 'other2', 'adult']
        region_list = ['Region 1', 'Region 2', '...', 'Region N']
        mig_dict = fisheries_io.read_migration_tables(
            args, class_list, region_list)
        test_matrix_dict = fisheries_io._parse_migration_tables(
            args, ['larva'])
        # pp.pprint(test_matrix_dict)
        # pp.pprint(mig_dict)
        testing.assert_array_equal(
            mig_dict['Migration'][0], test_matrix_dict['larva'])


class TestSingleParamsIO(unittest.TestCase):
    def test_verify_single_params(self):
        args = {
            'workspace_dir': '',
            'aoi_uri': None,
            'population_type': None,
            'sexsp': 1,
            'do_batch': False,
            'total_init_recruits': -1.0,
            'total_timesteps': -1,
            'recruitment_type': 'Ricker',
            'spawn_units': 'Individuals',
            'alpha': None,
            'beta': None,
            'total_recur_recruits': None,
            'migr_cont': True,
            'harvest_units': None,
            'frac_post_process': None,
            'unit_price': None,
            'val_cont': True,
            }

        # Check that path exists and user has read/write permissions along path
        with self.assertRaises(OSError):
            fisheries_io._verify_single_params(args)

        # Check timesteps positive number
        with self.assertRaises(ValueError):
            fisheries_io._verify_single_params(args, create_outputs=False)
        args['total_timesteps'] = 100

        # Check total_init_recruits for non-negative float
        with self.assertRaises(ValueError):
            fisheries_io._verify_single_params(args, create_outputs=False)
        args['total_init_recruits'] = 1.2

        # Check recruitment type's corresponding parameters exist
        with self.assertRaises(ValueError):
            fisheries_io._verify_single_params(args, create_outputs=False)
        args['alpha'] = -1.0
        args['beta'] = -1.0
        args['total_recur_recruits'] = -1.0

        # If BH or Ricker: Check alpha positive float
        with self.assertRaises(ValueError):
            fisheries_io._verify_single_params(args, create_outputs=False)
        args['alpha'] = 1.0

        # Check positive beta positive float
        with self.assertRaises(ValueError):
            fisheries_io._verify_single_params(args, create_outputs=False)
        args['beta'] = 1.0

        # Check total_recur_recruits is non-negative float
        args['recruitment_type'] = 'Fixed'
        with self.assertRaises(ValueError):
            fisheries_io._verify_single_params(args, create_outputs=False)
        args['total_recur_recruits'] = 100.0

        # If Harvest: Check frac_post_process float between [0,1]
        with self.assertRaises(ValueError):
            fisheries_io._verify_single_params(args, create_outputs=False)
        args['frac_post_process'] = 0.2

        # If Harvest: Check unit_price non-negative float
        with self.assertRaises(ValueError):
            fisheries_io._verify_single_params(args, create_outputs=False)
        args['unit_price'] = 20.2

        # Check file extension? (maybe try / except would be better)
        # Check shapefile subregions match regions in population parameters file
        args['aoi_uri'] = None


class TestFetchArgs(unittest.TestCase):
    def test_fetch_args(self):
        csv_uri = os.path.join(data_directory, 'CSVs/TestCSV_SN_Syntax.csv')
        mig_uri = os.path.join(data_directory, 'migration/')
        args = {
            'population_csv_uri': csv_uri,
            'migr_cont': True,
            'migration_dir': mig_uri,
            'workspace_dir': '',
            'aoi_uri': None,
            'population_type': "Stage-Based",
            'sexsp': 'No',
            'do_batch': False,
            'total_init_recruits': 1.2,
            'total_timesteps': 100,
            'recruitment_type': 'Ricker',
            'spawn_units': 'Individuals',
            'alpha': 1.0,
            'beta': 1.2,
            'total_recur_recruits': 100.0,
            'migr_cont': True,
            'harvest_units': "Weight",
            'frac_post_process': 0.2,
            'unit_price': 20.2,
            'val_cont': True,
        }
        vars_dict = fisheries_io.fetch_args(args, create_outputs=False)
        # pp.pprint(vars_dict)
        # with self.assertRaises():
        #    fisheries_io.fetch_args(args)

    def test_fetch_args2(self):
        csv_dir = os.path.join(data_directory, 'CSVs/Multiple_CSV_Test')
        mig_uri = os.path.join(data_directory, 'migration/')
        workspace_dir = ''
        args = {
            'population_csv_dir': csv_dir,
            'migr_cont': True,
            'migration_dir': mig_uri,
            'workspace_dir': workspace_dir,
            'aoi_uri': None,
            'population_type': "Stage-Based",
            'sexsp': 'No',
            'do_batch': True,
            'total_init_recruits': 1.2,
            'total_timesteps': 100,
            'recruitment_type': 'Ricker',
            'spawn_units': 'Individuals',
            'alpha': 1.0,
            'beta': 1.2,
            'total_recur_recruits': 100.0,
            'migr_cont': True,
            'harvest_units': "Weight",
            'frac_post_process': 0.2,
            'unit_price': 20.2,
            'val_cont': True,
        }
        # model_list = fisheries_io.fetch_args(args)
        # pp.pprint(model_list)
        # with self.assertRaises():
        #    fisheries_io.fetch_args(args)
        # os.removedirs(os.path.join(args['workspace_dir'], 'output'))


class TestCreateCSV(unittest.TestCase):
    def setUp(self):
        self.vars_dict = {
            'workspace_dir': 'path/to/workspace_dir',
            'output_dir': os.getcwd(),
            # 'aoi_uri': 'path/to/aoi_uri',
            'total_timesteps': 15,
            'population_type': 'Age-Based',
            'sexsp': 2,
            'do_batch': False,
            'spawn_units': 'Weight',
            'total_init_recruits': 100.0,
            'recruitment_type': 'Fixed',
            'alpha': 3.0,
            'beta': 4.0,
            'total_recur_recruits': 1.0,
            'migr_cont': True,
            'val_cont': True,
            'harvest_units': 'Individuals',
            'frac_post_process': 0.5,
            'unit_price': 5.0,

            # Pop Params
            # 'population_csv_uri': 'path/to/csv_uri',
            'Survnaturalfrac': np.ones([2, 2, 2]) * 0.5,  # Regions, Sexes, Classes
            'Classes': np.array(['larva', 'adult']),
            'Vulnfishing': np.array([[0.5, 0.5], [0.5, 0.5]]),
            'Maturity': np.array([[0.0, 1.0], [0.0, 1.0]]),
            'Duration': np.array([[2, 3], [2, 3]]),
            'Weight': np.array([[0.1, 1.0], [0.1, 1.0]]),
            'Fecundity': np.array([[0.1, 1.0], [0.1, 2.0]]),
            'Regions': np.array(['r1', 'r2']),
            'Exploitationfraction': np.array([0.25, 0.5]),
            'Larvaldispersal': np.array([0.5, 0.5]),

            # Mig Params
            # 'migration_dir': 'path/to/mig_dir',
            'Migration': [np.eye(2), np.eye(2)],

            # Derived Params
            'equilibrate_cycle': 10,
            'Survtotalfrac': np.array([[[0.5, 0.5], [0.5, 0.5]], [[0.5, 0.5], [0.5, 0.5]]]),  # Index Order: class, sex, region
            'G_survtotalfrac': np.ones([2, 2, 2]),  # (same)
            'P_survtotalfrac': np.ones([2, 2, 2]),  # (same)
            'N_tasx': np.ones([15, 2, 2, 2]),  # Index Order: time, class, sex, region
            'H_tx': np.ones([15, 2]),
            'V_tx': np.ones([15, 2]) * 5.0,
        }
        pass

    def test_create_csv(self):
        # fisheries_io._create_csv(self.vars_dict)
        pass


class TestCreateHTML(unittest.TestCase):
    def setUp(self):
        self.vars_dict = {
            'workspace_dir': 'path/to/workspace_dir',
            'output_dir': os.getcwd(),
            # 'aoi_uri': 'path/to/aoi_uri',
            'total_timesteps': 15,
            'population_type': 'Age-Based',
            'sexsp': 2,
            'do_batch': False,
            'spawn_units': 'Weight',
            'total_init_recruits': 100.0,
            'recruitment_type': 'Fixed',
            'alpha': 3.0,
            'beta': 4.0,
            'total_recur_recruits': 1.0,
            'migr_cont': True,
            'val_cont': True,
            'harvest_units': 'Individuals',
            'frac_post_process': 0.5,
            'unit_price': 5.0,

            # Pop Params
            # 'population_csv_uri': 'path/to/csv_uri',
            'Survnaturalfrac': np.ones([2, 2, 2]) * 0.5,  # Regions, Sexes, Classes
            'Classes': np.array(['larva', 'adult']),
            'Vulnfishing': np.array([[0.5, 0.5], [0.5, 0.5]]),
            'Maturity': np.array([[0.0, 1.0], [0.0, 1.0]]),
            'Duration': np.array([[2, 3], [2, 3]]),
            'Weight': np.array([[0.1, 1.0], [0.1, 1.0]]),
            'Fecundity': np.array([[0.1, 1.0], [0.1, 2.0]]),
            'Regions': np.array(['r1', 'r2']),
            'Exploitationfraction': np.array([0.25, 0.5]),
            'Larvaldispersal': np.array([0.5, 0.5]),

            # Mig Params
            # 'migration_dir': 'path/to/mig_dir',
            'Migration': [np.eye(2), np.eye(2)],

            # Derived Params
            'equilibrate_cycle': 10,
            'Survtotalfrac': np.array([[[0.5, 0.5], [0.5, 0.5]], [[0.5, 0.5], [0.5, 0.5]]]),  # Index Order: class, sex, region
            'G_survtotalfrac': np.ones([2, 2, 2]),  # (same)
            'P_survtotalfrac': np.ones([2, 2, 2]),  # (same)
            'N_tasx': np.ones([15, 2, 2, 2]),  # Index Order: time, class, sex, region
            'H_tx': np.ones([15, 2]),
            'V_tx': np.ones([15, 2]) * 5.0,
        }
        pass

    def test_create_html(self):
        # fisheries_io._create_html(self.vars_dict)
        pass


class TestCreateAOI(unittest.TestCase):
    def setUp(self):
        self.vars_dict = {
            'workspace_dir': 'path/to/workspace_dir',
            'output_dir': os.getcwd(),
            'aoi_uri': os.path.join(data_directory, 'Galveston_Subregion.shp'),
            'Classes': np.array(['larva']),
            'Regions': np.array(['1']),
            'N_tasx': np.ones([15, 2, 2, 2]),
            'H_tx': np.ones([15, 1]),
            'V_tx': np.ones([15, 1]) * 5.0,
        }
        pass

    def test_create_aoi(self):
        # fisheries_io._create_aoi(self.vars_dict)
        pass


if __name__ == '__main__':
    unittest.main()
