'''
python -m unittest test_crop_production_model
'''

import unittest
import os
import pprint
import tempfile

import gdal
from numpy import testing
import numpy as np

import crop_production_model as model

workspace_dir = '../../test/invest-data/test/data/crop_production/'
input_dir = '../../test/invest-data/test/data/crop_production/input/'
pp = pprint.PrettyPrinter(indent=4)


class TestModelCreateObservedYieldMaps(unittest.TestCase):
    def setUp(self):
        self.vars_dict = {
            'lulc_map_uri': os.path.join(input_dir, 'lulc_map.tif'),
            'crop_lookup_dict': {
                1: 'corn',
            },
            'observed_yields_maps_dict': {
                'corn': os.path.join(  # global raster of 5's
                    input_dir,
                    'spatial_dataset/observed_yield/corn_yield_map.tif'),
            },
            'tmp_observed_dir': tempfile.mkdtemp(),
        }

    def test_run(self):
        guess = model.create_observed_yield_maps(self.vars_dict)
        corn_yield_map_uri = os.path.join(
            guess['tmp_observed_dir'],
            'yield/corn_yield_map.tif')
        a = read_raster(corn_yield_map_uri)
        pp.pprint(a)


class TestModelCreatePercentileYieldMaps(unittest.TestCase):
    def setUp(self):
        self.vars_dict = {

        }

    def test_run1(self):
        guess = model.create_percentile_yield_maps(self.vars_dict)
        pass

    def test_run2(self):
        guess = model.create_percentile_yield_maps(self.vars_dict)
        pass


class TestModelCreateRegressionYieldMaps(unittest.TestCase):
    def setUp(self):
        self.vars_dict = {

        }

    def test_run1(self):
        guess = model.create_regression_yield_maps(self.vars_dict)
        pass

    def test_run2(self):
        guess = model.create_regression_yield_maps(self.vars_dict)
        pass


class TestModelCreateProductionMaps(unittest.TestCase):
    def setUp(self):
        self.vars_dict = {

        }

    def test_run1(self):
        guess = model.create_production_maps(self.vars_dict)
        pass

    def test_run2(self):
        guess = model.create_production_maps(self.vars_dict)
        pass


class TestModelCreateEconomicReturnsMap(unittest.TestCase):
    def setUp(self):
        self.vars_dict = {

        }

    def test_run1(self):
        guess = model.create_economic_returns_map(self.vars_dict)
        pass

    def test_run2(self):
        guess = model.create_economic_returns_map(self.vars_dict)
        pass


class TestModelCreateResultsTable(unittest.TestCase):
    def setUp(self):
        self.vars_dict = {

        }

    def test_run1(self):
        guess = model.create_results_table(self.vars_dict)
        pass

    def test_run2(self):
        guess = model.create_results_table(self.vars_dict)
        pass


# class TestModelSumCellsInRaster(unittest.TestCase):
#     def setUp(self):
#         self.vars_dict = {

#         }

#     def test_run1(self):
#         guess = model.sum_cells_in_raster(self.vars_dict)
#         pass

#     def test_run2(self):
#         guess = model.sum_cells_in_raster(self.vars_dict)
#         pass

def read_raster(raster_uri):
    dataset = gdal.Open(raster_uri)
    band = dataset.GetRasterBand(1)
    a = band.ReadAsArray()
    dataset = None
    band = None
    return a

if __name__ == '__main__':
    unittest.main()
