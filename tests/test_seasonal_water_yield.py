"""InVEST Seasonal water yield model tests that use the InVEST sample data"""

import unittest
import tempfile
import shutil
import os

import numpy
from osgeo import ogr
from pygeoprocessing.testing import scm

SAMPLE_DATA = os.path.join(
    os.path.dirname(__file__), '..', 'data', 'invest-data', 'seasonal_water_yield')
REGRESSION_DATA = os.path.join(
    os.path.dirname(__file__), 'data', 'seasonal_water_yield')


class SeasonalWaterYieldRegressionTests(unittest.TestCase):
    """Regression tests for the InVEST Seasonal Water Yield model"""

    def setUp(self):
        # this lets us delete the workspace after its done no matter the
        # the rest result
        self.workspace_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.workspace_dir)

    @staticmethod
    def generate_base_args(workspace_dir):
        """Generate an args list that is consistent across all three regression
        tests"""
        args = {
            'alpha_m': '1/12',
            'aoi_path': os.path.join(SAMPLE_DATA, 'watershed.shp'),
            'beta_i': '1.0',
            'biophysical_table_path': os.path.join(
                SAMPLE_DATA, 'biophysical_table.csv'),
            'dem_raster_path': os.path.join(SAMPLE_DATA, 'dem.tif'),
            'et0_dir': os.path.join(SAMPLE_DATA, 'eto_dir'),
            'gamma': '1.0',
            'lulc_raster_path': os.path.join(SAMPLE_DATA, 'lulc.tif'),
            'precip_dir': os.path.join(SAMPLE_DATA, 'precip_dir'),
            'rain_events_table_path': os.path.join(
                SAMPLE_DATA, 'rain_events_table.csv'),
            'results_suffix': '',
            'soil_group_path': os.path.join(SAMPLE_DATA, 'soil_group.tif'),
            'threshold_flow_accumulation': '1000',
            'workspace_dir': workspace_dir,
        }
        return args

    @scm.skip_if_data_missing(SAMPLE_DATA)
    @scm.skip_if_data_missing(REGRESSION_DATA)
    def test_base_regression(self):
        """SWY base regression test

        Executes SWY in default mode and checks that the output files are
        generated and that the aggregate shapefile fields are the same as the
        regression case."""
        from natcap.invest.seasonal_water_yield import seasonal_water_yield

        # use predefined directory so test can clean up files during teardown
        args = SeasonalWaterYieldRegressionTests.generate_base_args(
            self.workspace_dir)
        # make args explicit that this is a base run of SWY
        args['user_defined_climate_zones'] = False
        args['user_defined_local_recharge'] = False
        args['results_suffix'] = ''

        seasonal_water_yield.execute(args)

        SeasonalWaterYieldRegressionTests._assert_regression_results_equal(
            args['workspace_dir'],
            os.path.join(REGRESSION_DATA, 'file_list_base.txt'),
            os.path.join(args['workspace_dir'], 'aggregated_results.shp'),
            os.path.join(REGRESSION_DATA, 'agg_results_base.csv'))

    def test_climate_zones_regression(self):
        """SWY climate zone regression test

        Executes SWY in climate zones mode and checks that the output files are
        generated and that the aggregate shapefile fields are the same as the
        regression case."""
        from natcap.invest.seasonal_water_yield import seasonal_water_yield

        # use predefined directory so test can clean up files during teardown
        args = SeasonalWaterYieldRegressionTests.generate_base_args(
            self.workspace_dir)
        # modify args to account for climate zones defined
        args['climate_zone_raster_path'] = os.path.join(
            SAMPLE_DATA, 'climate_zones.tif')
        args['climate_zone_table_path'] = os.path.join(
            SAMPLE_DATA, 'climate_zone_events.csv')
        args['user_defined_climate_zones'] = True
        args['user_defined_local_recharge'] = False
        args['results_suffix'] = 'cz'

        seasonal_water_yield.execute(args)

        SeasonalWaterYieldRegressionTests._assert_regression_results_equal(
            args['workspace_dir'],
            os.path.join(REGRESSION_DATA, 'file_list_cz.txt'),
            os.path.join(
                args['workspace_dir'], 'aggregated_results_cz.shp'),
            os.path.join(REGRESSION_DATA, 'agg_results_cz.csv'))

    def test_user_recharge(self):
        """SWY user recharge regression test.

        Executes SWY in user defined local recharge mode and checks that the
        output files are generated and that the aggregate shapefile fields
        are the same as the regression case."""
        from natcap.invest.seasonal_water_yield import seasonal_water_yield

        # use predefined directory so test can clean up files during teardown
        args = SeasonalWaterYieldRegressionTests.generate_base_args(
            self.workspace_dir)
        # modify args to account for user recharge
        args['user_defined_climate_zones'] = False
        args['user_defined_local_recharge'] = True
        args['results_suffix'] = ''
        args['l_path'] = os.path.join(REGRESSION_DATA, 'L.tif')

        seasonal_water_yield.execute(args)

        SeasonalWaterYieldRegressionTests._assert_regression_results_equal(
            args['workspace_dir'],
            os.path.join(REGRESSION_DATA, 'file_list_user_recharge.txt'),
            os.path.join(args['workspace_dir'], 'aggregated_results.shp'),
            os.path.join(REGRESSION_DATA, 'agg_results_base.csv'))

    @staticmethod
    def _assert_regression_results_equal(
            workspace_dir, file_list_path, result_vector_path,
            agg_results_path):
        """Test the state of the workspace against the expected list of files
        and aggregated results.

        Parameters:
            workspace_dir (string): path to the completed model workspace
            file_list_path (string): path to a file that has a list of all
                the expected files relative to the workspace base
            result_vector_path (string): path to the summary shapefile
                produced by the SWY model.
            agg_results_path (string): path to a csv file that has the
                expected aggregated_results.shp table in the form of
                fid,vri_sum,qb_val per line

        Returns:
            None

        Raises:
            AssertionError if any files are missing or results are out of
            range by `tolerance_places`
        """

         # Test that the workspace has the same files as we expect
        SeasonalWaterYieldRegressionTests._test_same_files(
            file_list_path, workspace_dir)

        # we expect a file called 'aggregated_results.shp'
        result_vector = ogr.Open(result_vector_path)
        result_layer = result_vector.GetLayer()

        # The tolerance of 3 digits after the decimal was determined by
        # experimentation on the application with the given range of numbers.
        # This is an apparently reasonable approach as described by ChrisF:
        # http://stackoverflow.com/a/3281371/42897
        # and even more reading about picking numerical tolerance (it's hard):
        # https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
        tolerance_places = 3

        with open(agg_results_path, 'rb') as agg_result_file:
            for line in agg_result_file:
                fid, vri_sum, qb_val = [float(x) for x in line.split(',')]
                feature = result_layer.GetFeature(int(fid))
                for field, value in [('vri_sum', vri_sum), ('qb', qb_val)]:
                    numpy.testing.assert_almost_equal(
                        feature.GetField(field), value,
                        decimal=tolerance_places)

    @staticmethod
    def _test_same_files(base_list_path, directory_path):
        """Assert that the files listed in `base_list_file` are also in the
        directory pointed to by `path`.

        Parameters:
            base_list_file (string): a path to a file that has one relative file
                path per line.
            directory_path (string): a path to a directory whose contents will be
                checked against the files listed in `base_list_file`

        Returns:
            None

        Raises:
            AssertionError when there are files listed in `base_list_file` that
                don't exist in the directory indicated by `path`"""

        missing_files = []
        with open(base_list_path, 'r') as file_list:
            for file_path in file_list:
                full_path = os.path.join(directory_path, file_path.rstrip())
                if full_path == '':
                    #skip blank lines
                    continue
                if not os.path.isfile(full_path):
                    missing_files.append(full_path)
        if len(missing_files) > 0:
            raise AssertionError(
                "The following files were expected but not found: " +
                '\n'.join(missing_files))
