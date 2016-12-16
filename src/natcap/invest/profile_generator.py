"""InVEST Bathymetry Profile Generator."""
import os
import logging

import scipy.interpolate
from osgeo import gdal
from osgeo import osr
from osgeo import ogr
import numpy
import rtree
import pygeoprocessing

from . import utils

LOGGER = logging.getLogger('natcap.invest.profile_generator')
_OUTPUT_BASE_FILES = {
    }

_INTERMEDIATE_BASE_FILES = {
    'land_mask': 'land_mask.tif',
    'shore_mask': 'shore_mask.tif',
    'shore_convolution': 'shore_convolution.tif',
    'shore_points': 'shore_points.shp',
    'profile_lines': 'profile_lines.shp',
    'sample_points': 'sample_points.shp'
    }

_TMP_BASE_FILES = {
    'shore_kernel': 'shore_kernel.tif',
    'clipped_bathymetry': 'clipped_bathymetry.tif'
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

    bathymetry_pixel_size = pygeoprocessing.get_cell_size_from_uri(
        args['bathymetry_path'])
    pygeoprocessing.vectorize_datasets(
        [args['bathymetry_path']], _land_mask_op, f_reg['land_mask'],
        gdal.GDT_Int16, _MASK_NODATA, bathymetry_pixel_size, 'intersection',
        vectorize_op=False, datasets_are_pre_aligned=True)

    _make_shore_kernel(f_reg['shore_kernel'])

    pygeoprocessing.convolve_2d_uri(
        f_reg['land_mask'], f_reg['shore_kernel'], f_reg['shore_convolution'])

    shore_convolution_nodata = pygeoprocessing.get_nodata_from_uri(
        f_reg['shore_convolution'])

    def _shore_mask(shore_convolution):
        result = numpy.empty(shore_convolution.shape, dtype=numpy.int16)
        result[:] = _MASK_NODATA
        valid_mask = shore_convolution != shore_convolution_nodata
        result[valid_mask] = numpy.where(
            (shore_convolution[valid_mask] >= 9) &
            (shore_convolution[valid_mask] < 17), 1, _MASK_NODATA)
        return result
    pygeoprocessing.vectorize_datasets(
        [f_reg['shore_convolution']], _shore_mask, f_reg['shore_mask'],
        gdal.GDT_Int16, _MASK_NODATA, bathymetry_pixel_size, 'intersection',
        vectorize_op=False, datasets_are_pre_aligned=True)

    shore_raster = gdal.Open(f_reg['shore_mask'])
    shore_geotransform = shore_raster.GetGeoTransform()

    esri_driver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(f_reg['shore_points']):
        os.remove(f_reg['shore_points'])
    shore_points = esri_driver.CreateDataSource(f_reg['shore_points'])

    target_sr = osr.SpatialReference(shore_raster.GetProjection())
    shore_point_layer = shore_points.CreateLayer(
        'shore_points', srs=target_sr, geom_type=ogr.wkbPoint)
    shore_point_layer_defn = shore_point_layer.GetLayerDefn()

    # PUT SHORELINE PIXELS INTO R-TREE
    shore_point_index = rtree.index.Index()
    LOGGER.info('Building spatial index for shore points')
    for offset_info, data_block in pygeoprocessing.iterblocks(
            f_reg['shore_mask']):
        row_indexes, col_indexes = numpy.mgrid[
            offset_info['yoff']:offset_info['yoff']+offset_info['win_ysize'],
            offset_info['xoff']:offset_info['xoff']+offset_info['win_xsize']]
        valid_mask = data_block == 1
        x_coordinates = (
            shore_geotransform[0] +
            shore_geotransform[1] * (col_indexes[valid_mask] + 0.5) +
            shore_geotransform[2] * (row_indexes[valid_mask] + 0.5))
        y_coordinates = (
            shore_geotransform[3] +
            shore_geotransform[4] * (col_indexes[valid_mask] + 0.5) +
            shore_geotransform[5] * (row_indexes[valid_mask] + 0.5))
        for x_coord, y_coord in zip(x_coordinates, y_coordinates):
            point_feature = ogr.Feature(shore_point_layer_defn)
            shore_point_index.insert(
                0, (x_coord, y_coord), obj=(x_coord, y_coord))
            point_geometry = ogr.Geometry(ogr.wkbPoint)
            point_geometry.AddPoint(x_coord, y_coord)
            point_feature.SetGeometry(point_geometry)
            shore_point_layer.CreateFeature(point_feature)

    LOGGER.info("Constructing offshore profiles")
    sample_point_vector = ogr.Open(args['sample_point_vector_path'])
    for sample_point_layer in sample_point_vector:
        for sample_point in sample_point_layer:
            LOGGER.debug(sample_point)
            sample_point_geometry = sample_point.GetGeometryRef()
            sample_point = sample_point_geometry.GetPoint(0)
            closest_point = shore_point_index.nearest(
                (sample_point[0], sample_point[1]),
                objects=True).next().object
            LOGGER.debug("Closest point: %s", closest_point)

            if os.path.exists(f_reg['sample_points']):
                os.remove(f_reg['sample_points'])
            sample_points_vector = esri_driver.CreateDataSource(
                f_reg['sample_points'])
            target_sr = osr.SpatialReference(shore_raster.GetProjection())
            sample_points_layer = sample_points_vector.CreateLayer(
                'sample_points', srs=target_sr, geom_type=ogr.wkbPoint)
            sample_points_layer_defn = sample_points_layer.GetLayerDefn()

            vector_length = (
                (sample_point[0]-closest_point[0]) ** 2 +
                (sample_point[1]-closest_point[1]) ** 2) ** 0.5
            x_size = (sample_point[0]-closest_point[0]) / vector_length
            y_size = (sample_point[1]-closest_point[1]) / vector_length

            step_size = 10
            length = 2000

            for step in numpy.arange(0, length, step_size):
                point_feature = ogr.Feature(sample_points_layer_defn)
                sample_point_geometry = ogr.Geometry(ogr.wkbPoint)
                sample_point_geometry.AddPoint(
                    closest_point[0] + x_size * step,
                    closest_point[1] + y_size * step)
                point_feature.SetGeometry(sample_point_geometry)
                sample_points_layer.CreateFeature(point_feature)
            sample_points_layer.SyncToDisk()
            sample_points_layer = None
            sample_points_vector = None

            if os.path.exists(f_reg['profile_lines']):
                os.remove(f_reg['profile_lines'])
            profile_lines_vector = esri_driver.CreateDataSource(
                f_reg['profile_lines'])
            target_sr = osr.SpatialReference(shore_raster.GetProjection())
            profile_lines_layer = profile_lines_vector.CreateLayer(
                'profile_lines', srs=target_sr, geom_type=ogr.wkbLineString)
            profile_lines_layer_defn = profile_lines_layer.GetLayerDefn()
            line_feature = ogr.Feature(profile_lines_layer_defn)
            profile_line_geometry = ogr.Geometry(ogr.wkbLineString)
            profile_line_geometry.AddPoint(
                sample_point[0], sample_point[1])
            profile_line_geometry.AddPoint(
                closest_point[0], closest_point[1])
            line_feature.SetGeometry(profile_line_geometry)
            profile_lines_layer.CreateFeature(line_feature)
            profile_lines_layer.SyncToDisk()
            profile_lines_layer = None
            profile_lines_vector = None

            pygeoprocessing.align_dataset_list(
                [args['bathymetry_path']], [f_reg['clipped_bathymetry']],
                ['bilinear'], bathymetry_pixel_size, 'intersection', -1,
                aoi_uri=f_reg['profile_lines'], all_touched=True)

            shore_geotransform = shore_raster.GetGeoTransform()

            clipped_bathymetry_raster = gdal.Open(f_reg['clipped_bathymetry'])
            clipped_bathymetry_band = (
                clipped_bathymetry_raster.GetRasterBand(1))
            clipped_bathymetry_array = clipped_bathymetry_band.ReadAsArray()

            clipped_bathymetry_gt = pygeoprocessing.get_geotransform_uri(
                f_reg['clipped_bathymetry'])

            row_indexes, col_indexes = numpy.mgrid[
                0:clipped_bathymetry_array.shape[0],
                0:clipped_bathymetry_array.shape[1]]
            x_coordinates = (
                clipped_bathymetry_gt[0] +
                clipped_bathymetry_gt[1] * (col_indexes + 0.5) +
                clipped_bathymetry_gt[2] * (row_indexes + 0.5))
            y_coordinates = (
                clipped_bathymetry_gt[3] +
                clipped_bathymetry_gt[4] * (col_indexes + 0.5) +
                clipped_bathymetry_gt[5] * (row_indexes + 0.5))
            LOGGER.debug(x_coordinates.shape)
            LOGGER.debug(y_coordinates.shape)
            LOGGER.debug(clipped_bathymetry_array.shape)
            interp_fn = scipy.interpolate.interp2d(
                x_coordinates, y_coordinates,
                clipped_bathymetry_array)
            interp_fn(clipped_bathymetry_gt[0], clipped_bathymetry_gt[3])
    # GENERATE SHORELINE PIXELS
    # FOR EACH POINT:
    #   FIND NEAREST SHORELINE POINT
    #   CALCUALTE DIRECTION AS SHORE POINT TO SAMPLE POINT
    #   CREATE A LINE ROOTED AT SHORE PIXEL AS LONG AS REQUEST
    #   WALK ALONG LINE FOR EACH STEP:
    #       CALCULATE COORDINATE
    #       SAMPLE RASTER UNDERNEATH
    #       SAMPLE HABITAT LAYER UNDERNEATH


def _make_shore_kernel(kernel_path):
    """Makes a 3x3 raster with a 9 in the middle and 1s on the outside."""
    driver = gdal.GetDriverByName('GTiff')
    kernel_raster = driver.Create(
        kernel_path.encode('utf-8'), 3, 3, 1,
        gdal.GDT_Byte)

    # Make some kind of geotransform, it doesn't matter what but
    # will make GIS libraries behave better if it's all defined
    kernel_raster.SetGeoTransform([444720, 30, 0, 3751320, 0, -30])
    srs = osr.SpatialReference()
    srs.SetUTM(11, 1)
    srs.SetWellKnownGeogCS('NAD27')
    kernel_raster.SetProjection(srs.ExportToWkt())

    kernel_band = kernel_raster.GetRasterBand(1)
    kernel_band.SetNoDataValue(255)
    kernel_band.WriteArray(numpy.array([[1, 1, 1], [1, 9, 1], [1, 1, 1]]))
