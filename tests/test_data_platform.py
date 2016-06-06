"""Module for Testing Data Platform modules."""
import tempfile
import shutil
import unittest
import os

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
