import logging
from osgeo import gdal, ogr, osr
from natcap.invest import raster_utils
import sys

logging.basicConfig(format='%(asctime)s %(name)-20s %(levelname)-8s \
%(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

LOGGER = logging.getLogger('crop_production_modeled')

def execute(args):
    pass
