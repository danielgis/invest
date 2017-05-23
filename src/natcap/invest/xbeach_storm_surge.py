"""InVEST XBeach Storm Surge model."""
import os
import logging

from osgeo import gdal
from osgeo import ogr
import numpy

from . import utils

logging.basicConfig(format='%(asctime)s %(name)-20s %(levelname)-8s \
%(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

LOGGER = logging.getLogger('natcap.invest.xbeach_storm_surge')

def execute(args):
    """XBeach Storm Surge Model."""
    LOGGER.info("Nothing implemented %s.", args)
