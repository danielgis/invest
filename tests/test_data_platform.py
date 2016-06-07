"""Module for Testing Data Platform modules."""
import tempfile
import shutil
import unittest
import os
import collections

from pygeoprocessing.testing import scm

TEST_DATA = os.path.join(
    os.path.dirname(__file__), '..', 'data', 'invest-test-data',
    'data_platform')


class DataPlatformTests(unittest.TestCase):
    """Regression Tests for the Timber Model."""

    def setUp(self):
        """Overriding setUp function to create temp workspace directory."""
        # this lets us delete the workspace after its done no matter the
        # the rest result
        self.workspace_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Overriding tearDown function to remove temporary directory."""
        shutil.rmtree(self.workspace_dir)

    @scm.skip_if_data_missing(TEST_DATA)
    def test_data_platform_server(self):
        """Data Platform: testing the Data Platform Server."""
        from natcap.invest.data_platform import gdap_server

        database_filepath = os.path.join(self.workspace_dir, 'database.db')
        dataserver = gdap_server.DataServer(database_filepath)
        self.assertTrue(dataserver)
        dataserver.add_search_directory(
            [os.path.join(TEST_DATA, 'server_data')])

    def test_string_passed_as_list(self):
        """Data Platform: ensure ValueError on bad search_directory args."""
        from natcap.invest.data_platform import gdap_server

        database_filepath = os.path.join(self.workspace_dir, 'database.db')
        dataserver = gdap_server.DataServer(database_filepath)
        self.assertTrue(dataserver)
        self.assertRaises(
            ValueError, dataserver.add_search_directory, TEST_DATA)

    def test_coverage_polygon(self):
        """Data Platform: test database coverage given a polygon."""
        from natcap.invest.data_platform import gdap_server

        database_filepath = os.path.join(self.workspace_dir, 'database.db')
        dataserver = gdap_server.DataServer(database_filepath)
        self.assertTrue(dataserver)
        dataserver.add_search_directory(
            [os.path.join(TEST_DATA, 'server_data')])

        aoi_path = os.path.join(TEST_DATA, 'single_aoi')

        aoi_binary_zipstring = gdap_server.binaryzip_path(aoi_path)

        coverage_list = dataserver.get_data_coverage_polygon(
            aoi_binary_zipstring, [])

        # our sample data has two types and only 1 raster of each
        expected_coverage_types = {
            'dem': 1,
            'pawc': 1,
        }
        actual_coverage_types = collections.Counter(
            [_[1] for _ in coverage_list])
        self.assertEqual(actual_coverage_types, expected_coverage_types)
