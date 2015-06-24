'''
python -m unittest test_crop_production
'''

import unittest
import os
import pprint

from numpy import testing
import numpy as np

import crop_production

workspace_dir = '../../test/invest-data/Crop_Production'
input_dir = '../../test/invest-data/Crop_Production/input'
pp = pprint.PrettyPrinter(indent=4)


class TestOverallModel1(unittest.TestCase):
    def setUp(self):
        self.args = {

        }

    def test_run(self):
        guess = crop_production.execute(self.args, create_outputs=False)
        pass


class TestOverallModel2(unittest.TestCase):
    def setUp(self):
        self.args = {

        }

    def test_run(self):
        guess = crop_production.execute(self.args, create_outputs=False)
        pass


if __name__ == '__main__':
    unittest.main()
