"""InVEST Bathymetry Profile Generator."""
import os
import logging

from osgeo import gdal
from osgeo import ogr
import numpy

import pygeoprocessing
from . import utils

logging.basicConfig(format='%(asctime)s %(name)-20s %(levelname)-8s \
%(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

LOGGER = logging.getLogger('natcap.invest.profile_generator')
_OUTPUT_BASE_FILES = {
    }

_INTERMEDIATE_BASE_FILES = {
    'land_mask': 'land_mask.tif',
    'shore_mask': 'shore_mask.tif',
    }

_TMP_BASE_FILES = {
    }

_MASK_NODATA = -1

def execute(args):
    """Profile generator.

    Parameters:
        args['workspace_dir'] (string): output directory for intermediate,
            temporary, and final files
        args['results_suffix'] (string): (optional) string to append to any
            output file names
        args['bathymetry_path'] (string): path to a single band bathymetry
            raster that is projected in linear units and values represent
            elevations.
        args['shore_height'] (float): value in bathymetry raster that
            represents the shoreline elevation.  In most cases this would be
            0.
        args['sample_point_vector_path'] (string): Path to a point vector file
            that contains points from which to sample bathymetry.
        args['step_size'] (float): the number of linear units per step to
            sample the profile
        args['habitat_vector_path_list'] (list): List of paths to vector files
            that represent habitat layers.  The presence of overlap/no overlap
            will be included in the profile results.

    Returns:
        None.
    """
    file_suffix = utils.make_suffix_string(args, 'results_suffix')
    intermediate_output_dir = os.path.join(
        args['workspace_dir'], 'intermediate_outputs')
    output_dir = os.path.join(args['workspace_dir'])
    for dir_path in [output_dir, intermediate_output_dir]:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
    f_reg = utils.build_file_registry(
        [(_OUTPUT_BASE_FILES, output_dir),
         (_INTERMEDIATE_BASE_FILES, intermediate_output_dir),
         (_TMP_BASE_FILES, output_dir)], file_suffix)

    dem_nodata = pygeoprocessing.get_nodata_from_uri(args['bathymetry_path'])

    def _land_mask_op(dem):
        result = numpy.empty(dem.shape, dtype=numpy.int16)
        result[:] = _MASK_NODATA
        valid_mask = dem != dem_nodata
        result[valid_mask] = numpy.where(
            dem[valid_mask] >= args['shore_height'], 1, 0)
        return result

    dem_pixel_size = pygeoprocessing.get_cell_size_from_uri(
        args['bathymetry_path'])
    pygeoprocessing.vectorize_datasets(
        [args['bathymetry_path']], _land_mask_op, f_reg['land_mask'],
        gdal.GDT_Int16, _MASK_NODATA, dem_pixel_size, 'intersection',
        vectorize_op=False, datasets_are_pre_aligned=True)

    pygeoprocessing.convolve_2d_uri(
        f_reg['land_mask'], kernel_path, output_path)

    # GENERATE SHORELINE PIXELS
    # GENERATE SHORELINE SHAPE
    # SMOOTH SHORLINE SHAPE?
    # FOR EACH POINT:
    #   SNAP POINT TO NEAREST POINT ON LINE
    #   POINT IN DOWNWARD DIRECTION?
    #   WALK ALONG LINE FOR EACH STEP:
    #       CALCULATE COORDINATE
    #       SAMPLE RASTER UNDERNEATH
    #       SAMPLE HABITAT LAYER UNDERNEATH
