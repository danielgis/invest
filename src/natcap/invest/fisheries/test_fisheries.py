import unittest
import os
import pprint

from numpy import testing
import numpy as np
# from nose.plugins.skip import SkipTest

import fisheries

workspace_dir = '../../test/invest-data/Fisheries'
input_dir = '../../test/invest-data/Fisheries/input'
pp = pprint.PrettyPrinter(indent=4)


class TestBlueCrab(unittest.TestCase):
    def setUp(self):
        self.args = {
            'workspace_dir': workspace_dir,
            'results_suffix': '',
            'aoi_uri': os.path.join(input_dir, 'shapefile_galveston/Galveston_Subregion.shp'),
            'total_timesteps': 100,
            'population_type': 'Age-Based',
            'sexsp': 'No',
            'harvest_units': 'Individuals',
            'do_batch': False,
            'population_csv_uri': os.path.join(input_dir, 'input_blue_crab/population_params.csv'),
            'spawn_units': 'Individuals',
            'total_init_recruits': 200000.0,
            'recruitment_type': 'Ricker',
            'alpha': 6050000.0,
            'beta': 0.0000000414,
            'total_recur_recruits': 0,
            'migr_cont': False,
            'migration_dir': '',
            'val_cont': False,
            'frac_post_process': 0.0,
            'unit_price': 0.0,
        }

    def test_run(self):
        guess = fisheries.execute(self.args, create_outputs=False)
        # pp.pprint(guess[0]['N_tasx'])
        # pp.pprint(guess[0]['N_tasx'][0][0])

        # check harvest: 24,798,419
        harvest_guess = guess[0]['H_tx'][self.args['total_timesteps'] - 1].sum()
        testing.assert_approx_equal(harvest_guess, 24798419.0, significant=3)

        # check spawners: 42,644,460
        spawners_check = 42644460.0
        spawners_guess = guess[0]['Spawners_t'][self.args['total_timesteps'] - 1]
        testing.assert_approx_equal(spawners_guess, spawners_check, significant=4)


class TestDungenessCrab(unittest.TestCase):
    def setUp(self):
        self.args = {
            'workspace_dir': workspace_dir,
            'results_suffix': '',
            'aoi_uri': os.path.join(input_dir, 'shapefile_hood_canal/DC_HoodCanal_Subregions.shp'),
            'total_timesteps': 100,
            'population_type': 'Age-Based',
            'sexsp': 'Yes',
            'harvest_units': 'Individuals',
            'do_batch': False,
            'population_csv_uri': os.path.join(input_dir, 'input_dungeness_crab/population_params.csv'),
            'spawn_units': 'Individuals',
            'total_init_recruits': 2249339326901,
            'recruitment_type': 'Ricker',
            'alpha': 2000000,
            'beta': 0.000000309,
            'total_recur_recruits': 0,
            'migr_cont': False,
            'migration_dir': '',
            'val_cont': False,
            'frac_post_process': 0.0,
            'unit_price': 0.0,
        }

    def test_run(self):
        guess = fisheries.execute(self.args, create_outputs=False)
        # pp.pprint(guess[0]['Maturity'])

        # check harvest: 526,987
        harvest_check = 526987.0
        harvest_guess = guess[0]['H_tx'][self.args['total_timesteps'] - 1].sum()
        testing.assert_approx_equal(harvest_guess, harvest_check, significant=2)

        # check spawners: 4,051,538
        spawners_check = 4051538.0
        spawners_guess = guess[0]['Spawners_t'][self.args['total_timesteps'] - 1]
        testing.assert_approx_equal(spawners_guess, spawners_check, significant=3)


class TestLobster(unittest.TestCase):
    def setUp(self):
        self.args = {
            'workspace_dir': workspace_dir,
            'results_suffix': '',
            'aoi_uri': os.path.join(input_dir, 'shapefile_belize/Lob_Belize_Subregions.shp'),
            'total_timesteps': 100,
            'population_type': 'Age-Based',
            'sexsp': 'No',
            'harvest_units': 'Weight',
            'do_batch': False,
            'population_csv_uri': os.path.join(input_dir, 'input_lobster/population_params.csv'),
            'spawn_units': 'Weight',
            'total_init_recruits': 4686959.0,
            'recruitment_type': 'Beverton-Holt',
            'alpha': 5770000.0,
            'beta': 2885000.0,
            'total_recur_recruits': 0,
            'migr_cont': True,
            'migration_dir': os.path.join(input_dir, ''),
            'val_cont': True,
            'frac_post_process': 0.28633258,
            'unit_price': 29.93,
        }

    def test_run(self):
        guess = fisheries.execute(self.args, create_outputs=False)
        # pp.pprint(guess)

        # check harvest: 936,451
        harvest_guess = guess[0]['H_tx'][self.args['total_timesteps'] - 1].sum()
        testing.assert_approx_equal(harvest_guess, 963451.0, significant=4)

        # check spawners: 2,847,870
        spawners_guess = guess[0]['Spawners_t'][self.args['total_timesteps'] - 1]
        testing.assert_approx_equal(spawners_guess, 2847870.0, significant=3)


class TestShrimp(unittest.TestCase):
    def setUp(self):
        self.args = {
            'workspace_dir': workspace_dir,
            'results_suffix': '',
            'aoi_uri': os.path.join(input_dir, 'shapefile_galveston/Galveston_Subregion.shp'),
            'total_timesteps': 300,
            'population_type': 'Stage-Based',
            'sexsp': 'No',
            'harvest_units': 'Individuals',
            'do_batch': False,
            'population_csv_uri': os.path.join(input_dir, 'input_shrimp/population_params.csv'),
            'spawn_units': 'Weight',
            'total_init_recruits': 100000.0,
            'recruitment_type': 'Fixed',
            'alpha': 0,
            'beta': 0,
            'total_recur_recruits': 216000000000.0,
            'migr_cont': False,
            'migration_dir': '',
            'val_cont': False,
            'frac_post_process': 0.0,
            'unit_price': 0.0,
        }

    def test_run(self):
        guess = fisheries.execute(self.args, create_outputs=False)

        # check harvest: 456,424
        harvest_guess = guess[0]['H_tx'][self.args['total_timesteps'] - 1].sum()
        testing.assert_approx_equal(harvest_guess, 456424.0, significant=2)


class TestCustomRecruitmentFunc(unittest.TestCase):
    def setUp(self):
        Matu = np.array([0.0, 0.0, 0.0, 1.0])  # the Maturity vector in the Population Parameters File
        Weight = 1.0  # the Weight vector in the Population Parameters File
        LarvDisp = np.array([1.0])  # the LarvalDispersal vector in the Population Parameters File
        alpha = 6050000.0  # scalar value
        beta = 0.0000000414  # scalar value
        sexsp = 1   # 1 = not sex-specific, 2 = sex-specific

        def spawners(N_prev):
            return (N_prev * Matu * Weight).sum()

        def rec_func_Ricker(N_prev):
            N_0 = (LarvDisp * (alpha * spawners(N_prev) * (
                np.e ** (-beta * spawners(N_prev)))) / sexsp)
            return (N_0, spawners(N_prev))


        self.args = {
            'workspace_dir': workspace_dir,
            'results_suffix': '',
            'aoi_uri': os.path.join(input_dir, 'shapefile_galveston/Galveston_Subregion.shp'),
            'total_timesteps': 100,
            'population_type': 'Age-Based',
            'sexsp': 'No',
            'harvest_units': 'Individuals',
            'do_batch': False,
            'population_csv_uri': os.path.join(input_dir, 'input_blue_crab/population_params.csv'),
            'spawn_units': 'Individuals',
            'total_init_recruits': 200000.0,
            'recruitment_type': 'Other',
            'recruitment_func': rec_func_Ricker,
            'alpha': 6050000.0,
            'beta': 0.0000000414,
            'total_recur_recruits': 0,
            'migr_cont': False,
            'migration_dir': '',
            'val_cont': False,
            'frac_post_process': 0.0,
            'unit_price': 0.0,
        }

    def test_run(self):
        guess = fisheries.execute(self.args, create_outputs=False)
        # pp.pprint(guess[0]['N_tasx'])
        # pp.pprint(guess[0]['N_tasx'][0][0])

        # check harvest: 24,798,419
        harvest_guess = guess[0]['H_tx'][self.args['total_timesteps'] - 1].sum()
        testing.assert_approx_equal(harvest_guess, 24798419.0, significant=3)

        # check spawners: 42,644,460
        spawners_check = 42644460.0
        spawners_guess = guess[0]['Spawners_t'][self.args['total_timesteps'] - 1]
        testing.assert_approx_equal(spawners_guess, spawners_check, significant=4)


class TestCustomRecruitmentFunc2(unittest.TestCase):
    def setUp(self):
        Matu = np.array([[ 0.,  0.,  0.,  1.,  1.],
                         [ 0.,  0.,  0.,  0.,  0.]])  # the Maturity vector in the Population Parameters File
        Weight = 1.0  # the Weight vector in the Population Parameters File
        LarvDisp = np.array([ 0.09,  0.12,  0.18,  0.29,  0.17,  0.15])  # the LarvalDispersal vector in the Population Parameters File
        alpha = 2000000  # scalar value
        beta = 0.000000309  # scalar value
        sexsp = 2   # 1 = not sex-specific, 2 = sex-specific

        def spawners(N_prev):
            return (N_prev * Matu * Weight).sum()

        def rec_func_Ricker(N_prev):
            N_0 = (LarvDisp * (alpha * spawners(N_prev) * (
                np.e ** (-beta * spawners(N_prev)))) / sexsp)
            return (N_0, spawners(N_prev))

        self.args = {
            'workspace_dir': workspace_dir,
            'results_suffix': '',
            'aoi_uri': os.path.join(input_dir, 'shapefile_hood_canal/DC_HoodCanal_Subregions.shp'),
            'total_timesteps': 100,
            'population_type': 'Age-Based',
            'sexsp': 'Yes',
            'harvest_units': 'Individuals',
            'do_batch': False,
            'population_csv_uri': os.path.join(input_dir, 'input_dungeness_crab/population_params.csv'),
            'spawn_units': 'Individuals',
            'total_init_recruits': 2249339326901,
            'recruitment_type': 'Other',
            'recruitment_func': rec_func_Ricker,
            'alpha': 2000000,
            'beta': 0.000000309,
            'total_recur_recruits': 0,
            'migr_cont': False,
            'migration_dir': '',
            'val_cont': False,
            'frac_post_process': 0.0,
            'unit_price': 0.0,
        }

    def test_run(self):
        guess = fisheries.execute(self.args, create_outputs=False)
        # pp.pprint(guess)

        # check harvest: 526,987
        harvest_check = 526987.0
        harvest_guess = guess[0]['H_tx'][self.args['total_timesteps'] - 1].sum()
        testing.assert_approx_equal(harvest_guess, harvest_check, significant=2)

        # check spawners: 4,051,538
        spawners_check = 4051538.0
        spawners_guess = guess[0]['Spawners_t'][self.args['total_timesteps'] - 1]
        testing.assert_approx_equal(spawners_guess, spawners_check, significant=3)


if __name__ == '__main__':
    unittest.main()
