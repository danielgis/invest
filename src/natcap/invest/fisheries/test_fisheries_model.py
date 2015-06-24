import unittest
import pprint

from numpy import testing
import numpy as np

import fisheries_model as model

pp = pprint.PrettyPrinter(indent=4)


class TestInitializeVars(unittest.TestCase):
    def setUp(self):
        self.sample_vars = {
            #'workspace_dir': 'path/to/workspace_dir',
            #'aoi_uri': 'path/to/aoi_uri',
            'total_timesteps': 100,
            'results_suffix': '',
            'population_type': 'Stage-Based',
            'sexsp': 2,
            'spawn_units': 'Weight',
            'total_init_recruits': 100.0,
            'recruitment_type': 'Ricker',
            'alpha': 32.4,
            'beta': 54.2,
            'total_recur_recruits': 92.1,
            'migr_cont': True,
            'val_cont': True,
            'harvest_units': 'Individuals',
            'frac_post_process': 0.5,
            'unit_price': 5.0,

            # Pop Params
            #'population_csv_uri': 'path/to/csv_uri',
            'Survnaturalfrac': np.ones([2, 2, 2]),  # Regions, Sexes, Classes
            'Classes': np.array(['larva', 'adult']),
            'Vulnfishing': np.array([[0.5, 0.5], [0.5, 0.5]]),
            'Maturity': np.array([[0.0, 1.0], [0.0, 1.0]]),
            'Duration': np.array([[2, 3], [2, 3]]),
            'Weight': np.array([[0.1, 1.0], [0.1, 2.0]]),
            'Fecundity': np.array([[0.0, 1.0], [0.0, 1.0]]),
            'Regions': np.array(['r1', 'r2']),
            'Exploitationfraction': np.array([0.25, 0.5]),
            'Larvaldispersal': np.array([0.75, 0.75]),

            # Mig Params
            #'migration_dir': 'path/to/mig_dir',
            'Migration': [np.eye(2), np.eye(2)]
        }

    def test_calc_survtotalfrac(self):
        # Test very simple
        self.sample_vars['Survnaturalfrac'] = np.array([[[1, 0.5], [0.5, 1]], [[2, 1], [1, 2]]])
        self.sample_vars['Exploitationfraction'] = np.array([1.0, 2.0])
        self.sample_vars['Vulnfishing'] = np.array([[1.0, 2.0], [2.0, 1.0]])
        check = np.array([[[0, -0.5], [-0.5, 0]], [[-2, -3], [-3, -2]]])
        guess = model._calc_survtotalfrac(self.sample_vars)
        # print "Guess"
        # pp.pprint(guess)
        # print "Check"
        # pp.pprint(check)
        testing.assert_array_equal(guess, check)

    def test_p_g_survtotalfrac(self):
        # Test simple
        self.sample_vars['Survtotalfrac'] = np.array([[[1, 0.5], [0.5, 1]], [[2, 1], [1, 2]]])
        self.sample_vars['Exploitationfraction'] = np.array([1.0, 2.0])
        self.sample_vars['Vulnfishing'] = np.array([[1.0, 2.0], [2.0, 1.0]])
        self.sample_vars['Duration'] = np.array([[1.0, 2.0], [1.0, 3.0]])
        G_check = np.array([[[np.nan, 1.0/6], [0.5, np.nan]], [[2.0, np.nan], [np.nan, 8.0/7]]])
        P_check = np.array([[[np.nan, 1.0/3], [0.0, np.nan]], [[0.0, np.nan], [np.nan, 6.0/7]]])
        G_guess, P_guess = model._calc_p_g_survtotalfrac(self.sample_vars)
        # print "G_Guess"
        # pp.pprint(G_guess)
        # print "G_Check"
        # pp.pprint(G_check)
        # print "P_Guess"
        # pp.pprint(P_guess)
        # print "P_Check"
        # pp.pprint(P_check)
        testing.assert_array_equal(P_guess, P_check)

    def test_initialize_vars(self):
        # vars_dict = model.initialize_vars(self.sample_vars)
        # pp.pprint(vars_dict['Survtotalfrac'])
        # pp.pprint(vars_dict['G_survtotalfrac'])
        # pp.pprint(vars_dict['P_survtotalfrac'])
        pass


class TestSetRecruitmentFunc(unittest.TestCase):
    def setUp(self):
        self.sample_vars = {
            #'workspace_dir': 'path/to/workspace_dir',
            #'aoi_uri': 'path/to/aoi_uri',
            'total_timesteps': 100,
            'results_suffix': '',
            'population_type': 'Stage-Based',
            'sexsp': 2,
            'spawn_units': 'Weight',
            'total_init_recruits': 100.0,
            'recruitment_type': 'Ricker',
            'alpha': 3.0,
            'beta': 4.0,
            'total_recur_recruits': 92.1,
            'migr_cont': True,
            'val_cont': True,
            'harvest_units': 'Individuals',
            'frac_post_process': 0.5,
            'unit_price': 5.0,

            # Pop Params
            # 'population_csv_uri': 'path/to/csv_uri',
            'Survnaturalfrac': np.ones([2, 2, 2]),  # Regions, Sexes, Classes
            'Classes': np.array(['larva', 'adult']),
            'Vulnfishing': np.array([[0.5, 0.5], [0.5, 0.5]]),
            'Maturity': np.array([[0.0, 1.0], [0.0, 1.0]]),
            'Duration': np.array([[2, 3], [2, 3]]),
            'Weight': np.array([[0.1, 1.0], [0.1, 2.0]]),
            'Fecundity': np.array([[0.1, 1.0], [0.1, 2.0]]),
            'Regions': np.array(['r1', 'r2']),
            'Exploitationfraction': np.array([0.25, 0.5]),
            'Larvaldispersal': np.array([0.75, 0.75]),

            # Mig Params
            # 'migration_dir': 'path/to/mig_dir',
            'Migration': [np.eye(2), np.eye(2)],

            # Derived Params
            'Survtotalfrac': np.ones([2, 2, 2]),  # Index Order: region, sex, class
            'G_survtotalfrac': np.ones([2, 2, 2]),  # (same)
            'P_survtotalfrac': np.ones([2, 2, 2]),  # (same)
            'N_tasx': np.zeros([100, 2, 2, 2]),  # Index Order: time, region, sex, class
        }

    def test_spawners(self):
        def create_Spawners(Matu, Weight):
            return lambda N_prev: (N_prev * Matu * Weight).sum()

        N_prev = np.array(
            [[[5.0, 10.0], [2.0, 4.0]], [[4.0, 8.0], [1.0, 2.0]]])
        Matu = self.sample_vars['Maturity']
        Weight = self.sample_vars['Weight']
        f = create_Spawners(Matu, Weight)
        check = 30
        guess = f(N_prev)
        # print "Guess"
        # pp.pprint(guess)
        testing.assert_equal(guess, check)

    def test_set_recru_func(self):
        vars_dict = self.sample_vars
        N_prev = np.array(
            [[[5.0, 10.0], [2.0, 4.0]], [[4.0, 8.0], [1.0, 2.0]]])

        # Test B-H
        vars_dict['recruitment_type'] = 'Beverton-Holt'
        rec_func = model.set_recru_func(vars_dict)
        guess, spawners = rec_func(N_prev)
        check = np.array([(270.0 / 272.0), (270.0 / 272.0)])
        # print "Guess"
        # pp.pprint(guess)
        testing.assert_equal(guess, check)

        # Test Ricker
        vars_dict['recruitment_type'] = 'Ricker'
        rec_func = model.set_recru_func(vars_dict)
        guess, spawners = rec_func(N_prev)
        check = np.array([0.75, 0.75]) * (45.0 * np.e**-120.0)
        # print "Guess"
        # pp.pprint(guess)
        testing.assert_equal(guess, check)

        # Test Fecundity
        vars_dict['recruitment_type'] = 'Fecundity'
        rec_func = model.set_recru_func(vars_dict)
        guess, spawners = rec_func(N_prev)
        check = np.array([11.25, 11.25])
        # print "Guess"
        # pp.pprint(guess)
        testing.assert_equal(guess, check)

        # Test Fixed
        vars_dict['recruitment_type'] = 'Fixed'
        rec_func = model.set_recru_func(vars_dict)
        guess, spawners = rec_func(N_prev)
        check = np.array([0.75, 0.75]) * 92.1 / 2
        # print "Guess"
        # pp.pprint(guess)
        testing.assert_equal(guess, check)


class TestSetHarvestFunc(unittest.TestCase):
    def setUp(self):
        self.sample_vars = {
            #'workspace_dir': 'path/to/workspace_dir',
            #'aoi_uri': 'path/to/aoi_uri',
            'total_timesteps': 100,
            'results_suffix': '',
            'population_type': 'Stage-Based',
            'sexsp': 2,
            'spawn_units': 'Weight',
            'total_init_recruits': 100.0,
            'recruitment_type': 'Ricker',
            'alpha': 3.0,
            'beta': 4.0,
            'total_recur_recruits': 92.1,
            'migr_cont': True,
            'val_cont': True,
            'harvest_units': 'Individuals',
            'frac_post_process': 0.5,
            'unit_price': 5.0,

            # Pop Params
            # 'population_csv_uri': 'path/to/csv_uri',
            'Survnaturalfrac': np.ones([2, 2, 2]),  # Regions, Sexes, Classes
            'Classes': np.array(['larva', 'adult']),
            'Vulnfishing': np.array([[0.5, 0.5], [0.5, 0.5]]),
            'Maturity': np.array([[0.0, 1.0], [0.0, 1.0]]),
            'Duration': np.array([[2, 3], [2, 3]]),
            'Weight': np.array([[0.1, 1.0], [0.1, 2.0]]),
            'Fecundity': np.array([[0.1, 1.0], [0.1, 2.0]]),
            'Regions': np.array(['r1', 'r2']),
            'Exploitationfraction': np.array([0.5, 0.5]),
            'Larvaldispersal': np.array([0.75, 0.75]),

            # Mig Params
            # 'migration_dir': 'path/to/mig_dir',
            'Migration': [np.eye(2), np.eye(2)],

            # Derived Params
            'Survtotalfrac': np.ones([2, 2, 2]),  # Index Order: region, sex, class
            'G_survtotalfrac': np.ones([2, 2, 2]),  # (same)
            'P_survtotalfrac': np.ones([2, 2, 2]),  # (same)
            'N_tasx': np.ones([100, 2, 2, 2]),  # Index Order: time, region, sex, class
        }

    def test_set_harv_func(self):
        vars_dict = self.sample_vars
        harv_func = model.set_harvest_func(vars_dict)
        H_x_guess, V_x_guess = harv_func(vars_dict['N_tasx'][0])
        H_check = np.ones([2])
        # print "Harvest Guess"
        # print H_guess
        testing.assert_equal(H_x_guess, H_check)
        V_check = H_check * 2.5
        # print "Valuation Guess"
        # print V_guess
        testing.assert_equal(V_x_guess, V_check)


class TestSetInitCondFunc(unittest.TestCase):
    def setUp(self):
        self.sample_vars = {
            #'workspace_dir': 'path/to/workspace_dir',
            #'aoi_uri': 'path/to/aoi_uri',
            'total_timesteps': 100,
            'results_suffix': '',
            'population_type': 'Stage-Based',
            'sexsp': 2,
            'spawn_units': 'Weight',
            'total_init_recruits': 100.0,
            'recruitment_type': 'Ricker',
            'alpha': 3.0,
            'beta': 4.0,
            'total_recur_recruits': 10.0,
            'migr_cont': True,
            'val_cont': True,
            'harvest_units': 'Individuals',
            'frac_post_process': 0.5,
            'unit_price': 5.0,

            # Pop Params
            # 'population_csv_uri': 'path/to/csv_uri',
            'Survnaturalfrac': np.ones([2, 2, 2]),  # Regions, Sexes, Classes
            'Classes': np.array(['larva', 'adult']),
            'Vulnfishing': np.array([[0.5, 0.5], [0.5, 0.5]]),
            'Maturity': np.array([[0.0, 1.0], [0.0, 1.0]]),
            'Duration': np.array([[2, 3], [2, 3]]),
            'Weight': np.array([[0.1, 1.0], [0.1, 2.0]]),
            'Fecundity': np.array([[0.1, 1.0], [0.1, 2.0]]),
            'Regions': np.array(['r1', 'r2']),
            'Exploitationfraction': np.array([0.25, 0.5]),
            'Larvaldispersal': np.array([0.5, 0.5]),

            # Mig Params
            # 'migration_dir': 'path/to/mig_dir',
            'Migration': [np.matrix(np.eye(2)), np.matrix(np.eye(2))],

            # Derived Params
            'Survtotalfrac': np.array([[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]]]),  # Index Order: class, sex, region
            'G_survtotalfrac': np.ones([2, 2, 2]),  # (same)
            'P_survtotalfrac': np.ones([2, 2, 2]),  # (same)
            'N_tasx': np.ones([100, 2, 2, 2]),  # Index Order: time, class, sex, region
        }

    def test_stage_based(self):
        vars_dict = self.sample_vars
        init_cond_func = model.set_init_cond_func(vars_dict)
        N_0_guess = init_cond_func()
        N_0_check = np.array([[[25.0, 25.0], [25.0, 25.0]], [[1.0, 1.0], [1.0, 1.0]]])
        # print "N_0 Guess"
        # print N_0_guess
        testing.assert_equal(N_0_guess, N_0_check)

    def test_age_based(self):
        vars_dict = self.sample_vars
        vars_dict['population_type'] = 'Age-Based'
        init_cond_func = model.set_init_cond_func(vars_dict)
        N_0_guess = init_cond_func()
        N_0_check = np.array([[[25.0, 25.0], [25.0, 25.0]], [[25.0*1/-4, 25.0*2/-5], [(25.0*3/-6), (25.0*4/-7)]]])
        # print "N_0 Guess"
        # print N_0_guess
        testing.assert_equal(N_0_guess, N_0_check)

    def test_age_based_2(self):
        vars_dict = self.sample_vars
        vars_dict['Classes'] = np.array(['larva', 'middle', 'adult'])
        vars_dict['Survtotalfrac'] = np.array([[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]], [[9.0, 10.0], [11.0, 12.0]]])
        vars_dict['population_type'] = 'Age-Based'

        init_cond_func = model.set_init_cond_func(vars_dict)
        N_0_guess = init_cond_func()
        N_0_check = np.array([[[25.0, 25.0], [25.0, 25.0]],[[25.0*1, 25.0*2], [(25.0*3), (25.0*4)]],[[(25.0*1*5/-8), (25.0*2*6/-9)], [(25.0*3*7/-10), (25.0*4*8/-11)]]])
        testing.assert_equal(N_0_guess, N_0_check)


class TestSetCycleFunc(unittest.TestCase):
    def setUp(self):
        self.sample_vars = {
            # 'workspace_dir': 'path/to/workspace_dir',
            # 'aoi_uri': 'path/to/aoi_uri',
            'total_timesteps': 10,
            'results_suffix': '',
            'population_type': 'Stage-Based',
            'sexsp': 2,
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
            'Survnaturalfrac': np.ones([2, 2, 2]),  # Regions, Sexes, Classes
            'Classes': np.array(['larva', 'adult']),
            'Vulnfishing': np.array([[0.5, 0.5], [0.5, 0.5]]),
            'Maturity': np.array([[0.0, 1.0], [0.0, 1.0]]),
            'Duration': np.array([[2, 3], [2, 3]]),
            'Weight': np.array([[0.0, 1.0], [0.0, 1.0]]),
            'Fecundity': np.array([[0.1, 1.0], [0.1, 2.0]]),
            'Regions': np.array(['r1', 'r2']),
            'Exploitationfraction': np.array([0.25, 0.5]),
            'Larvaldispersal': np.array([0.5, 0.5]),

            # Mig Params
            # 'migration_dir': 'path/to/mig_dir',
            'Migration': [np.matrix(np.eye(2)), np.matrix(np.eye(2))],

            # Derived Params
            'Survtotalfrac': np.array([[[1.0, 2.0], [3.0, 4.0]], [[5.0, 6.0], [7.0, 8.0]]]),  # Index Order: class, sex, region
            'G_survtotalfrac': np.ones([2, 2, 2]),  # (same)
            'P_survtotalfrac': np.ones([2, 2, 2]),  # (same)
            'N_tasx': np.ones([10, 2, 2, 2]),  # Index Order: time, class, sex, region
        }

    def test_stage_based(self):
        vars_dict = self.sample_vars
        rec_func = model.set_recru_func(vars_dict)
        cycle_func = model.set_cycle_func(vars_dict, rec_func)

        N_prev = np.ones([2, 2, 2])

        N_cur_guess, spawners = cycle_func(N_prev)
        # N_cur_check = np.array([])
        # testing.assert_equal(N_cur_guess, N_cur_check)

    def test_age_based(self):
        pass


class TestRunPopulationModel(unittest.TestCase):
    def setUp(self):
        self.sample_vars = {
            # 'workspace_dir': 'path/to/workspace_dir',
            # 'aoi_uri': 'path/to/aoi_uri',
            'total_timesteps': 100,
            'results_suffix': '',
            'population_type': 'Age-Based',
            'sexsp': 2,
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
            'Migration': [np.matrix(np.eye(2)), np.matrix(np.eye(2))],

            # Derived Params
            'Survtotalfrac': np.array([[[0.5, 0.5], [0.5, 0.5]], [[0.5, 0.5], [0.5, 0.5]]]),  # Index Order: class, sex, region
            'G_survtotalfrac': np.ones([2, 2, 2]),  # (same)
            'P_survtotalfrac': np.ones([2, 2, 2]),  # (same)
            'N_tasx': np.ones([100, 2, 2, 2]),  # Index Order: time, class, sex, region
            'H_tx': np.ones([100, 2]),
            'V_tx': np.ones([100, 2]) * 5.0,
            'Spawners_t': np.zeros([100]),
        }

    def test_run_population_model(self):
        ## UNTESTED
        vars_dict = self.sample_vars
        recru_func = model.set_recru_func(vars_dict)
        init_cond_func = model.set_init_cond_func(vars_dict)
        cycle_func = model.set_cycle_func(vars_dict, recru_func)
        harvest_func = model.set_harvest_func(vars_dict)

        # Run Model
        vars_dict = model.run_population_model(
            vars_dict, init_cond_func, cycle_func, harvest_func)

        # pp.pprint(vars_dict['N_tasx'])


if __name__ == '__main__':
    unittest.main()
