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
        args['representative_point_vector_path'] (string): Path to a point
            vector file that contains points from which to sample bathymetry.
        args['step_size'] (float): the number of linear units per step to
            sample the profile.
        args['profile_length'] (float): the length of the profile from off
            the shore in linear units.
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

    bathymetry_nodata = pygeoprocessing.get_nodata_from_uri(
        args['bathymetry_path'])

    def _land_mask_op(bathymetry):
        """Mask values >= shore height."""
        result = numpy.empty(bathymetry.shape, dtype=numpy.int16)
        result[:] = _MASK_NODATA
        valid_mask = bathymetry != bathymetry_nodata
        result[valid_mask] = numpy.where(
            bathymetry[valid_mask] >= args['shore_height'], 1, 0)
        return result

    bathymetry_pixel_size = pygeoprocessing.get_cell_size_from_uri(
        args['bathymetry_path'])
    pygeoprocessing.vectorize_datasets(
        [args['bathymetry_path']], _land_mask_op, f_reg['land_mask'],
        gdal.GDT_Int16, _MASK_NODATA, bathymetry_pixel_size, 'intersection',
        vectorize_op=False, datasets_are_pre_aligned=True)

    bathymetry_raster = gdal.Open(args['bathymetry_path'])
    bathymetry_band = bathymetry_raster.GetRasterBand(1)
    _make_shore_kernel(f_reg['shore_kernel'])

    pygeoprocessing.convolve_2d_uri(
        f_reg['land_mask'], f_reg['shore_kernel'], f_reg['shore_convolution'])

    shore_convolution_nodata = pygeoprocessing.get_nodata_from_uri(
        f_reg['shore_convolution'])

    def _shore_mask(shore_convolution):
        """Mask values on land that border water."""
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

    if os.path.exists(f_reg['sample_points']):
        os.remove(f_reg['sample_points'])
    sample_points_vector = esri_driver.CreateDataSource(
        f_reg['sample_points'])
    target_sr = osr.SpatialReference(shore_raster.GetProjection())
    sample_points_layer = sample_points_vector.CreateLayer(
        'sample_points', srs=target_sr, geom_type=ogr.wkbPoint)
    sample_points_layer_defn = sample_points_layer.GetLayerDefn()

    LOGGER.info("Constructing offshore profiles")
    representative_point_vector = ogr.Open(
        args['representative_point_vector_path'])
    for representative_point_layer in representative_point_vector:
        for representative_point in representative_point_layer:
            representative_point_geometry = (
                representative_point.GetGeometryRef())
            representative_point = representative_point_geometry.GetPoint(0)
            closest_point = shore_point_index.nearest(
                (representative_point[0], representative_point[1]),
                objects=True).next().object

            vector_length = (
                (representative_point[0]-closest_point[0]) ** 2 +
                (representative_point[1]-closest_point[1]) ** 2) ** 0.5
            x_size = (
                (representative_point[0]-closest_point[0]) / vector_length)
            y_size = (
                (representative_point[1]-closest_point[1]) / vector_length)
            LOGGER.debug("%s, %s", x_size, y_size)
            sample_point_list = []
            for step in numpy.arange(
                    0, args['profile_length'], args['step_size']):
                point_feature = ogr.Feature(sample_points_layer_defn)
                sample_point_geometry = ogr.Geometry(ogr.wkbPoint)
                sample_point_list.append(
                    (closest_point[0] + x_size * step,
                     closest_point[1] + y_size * step))
                sample_point_geometry.AddPoint(
                    sample_point_list[-1][0], sample_point_list[-1][1])
                point_feature.SetGeometry(sample_point_geometry)
                sample_points_layer.CreateFeature(point_feature)

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
                closest_point[0], closest_point[1])
            profile_line_geometry.AddPoint(
                sample_point_list[-1][0], sample_point_list[-1][1])
            line_feature.SetGeometry(profile_line_geometry)
            profile_lines_layer.CreateFeature(line_feature)
            profile_lines_layer.SyncToDisk()
            extent = profile_lines_layer.GetExtent()
            # convert extent to bathymetry index extent
            bathymetry_gt = pygeoprocessing.get_geotransform_uri(
                args['bathymetry_path'])

            # reverse last two because y coord moves up while pixels move down
            extent_in_pixel_coords = (
                int((extent[0] - bathymetry_gt[0]) / bathymetry_gt[1]),
                int((extent[1] - bathymetry_gt[0]) / bathymetry_gt[1]),
                int(round(0.5+(extent[3] - bathymetry_gt[3]) /
                          bathymetry_gt[5])),
                int(round(0.5+(extent[2] - bathymetry_gt[3]) /
                          bathymetry_gt[5])))
            offset_dict = {
                'xoff': extent_in_pixel_coords[0],
                'yoff': extent_in_pixel_coords[2],
                'win_xsize': (
                    extent_in_pixel_coords[1]-extent_in_pixel_coords[0]),
                'win_ysize': (
                    extent_in_pixel_coords[3]-extent_in_pixel_coords[2]),
            }
            LOGGER.debug(bathymetry_gt)
            LOGGER.debug(extent)
            LOGGER.debug(extent_in_pixel_coords)
            LOGGER.debug(offset_dict)
            profile_lines_layer = None
            profile_lines_vector = None

            clipped_bathymetry_array = bathymetry_band.ReadAsArray(
                **offset_dict)

            x_coordinates = numpy.arange(clipped_bathymetry_array.shape[1])
            x_coordinates = (
                bathymetry_gt[0] +
                (x_coordinates + offset_dict['xoff']) * bathymetry_gt[1])

            y_coordinates = numpy.arange(clipped_bathymetry_array.shape[0])
            y_coordinates = (
                bathymetry_gt[3] +
                (y_coordinates + offset_dict['yoff']) * bathymetry_gt[5])[::-1]

            LOGGER.debug(
                "%s, %s, %s", x_coordinates.shape, y_coordinates.shape,
                clipped_bathymetry_array.shape)
            interp_fn = scipy.interpolate.RectBivariateSpline(
                y_coordinates, x_coordinates, clipped_bathymetry_array,
                kx=1, ky=1)
            for x_coord, y_coord in sample_point_list:
                value = interp_fn(y_coord, x_coord)
                LOGGER.debug("%s, %s, %s", y_coord, x_coord, value)
    # GENERATE SHORELINE PIXELS
    # FOR EACH POINT:
    #   FIND NEAREST SHORELINE POINT
    #   CALCUALTE DIRECTION AS SHORE POINT TO SAMPLE POINT
    #   CREATE A LINE ROOTED AT SHORE PIXEL AS LONG AS REQUEST
    #   WALK ALONG LINE FOR EACH STEP:
    #       CALCULATE COORDINATE
    #       SAMPLE RASTER UNDERNEATH
    #       SAMPLE HABITAT LAYER UNDERNEATH
    sample_points_layer.SyncToDisk()
    sample_points_layer = None
    sample_points_vector = None


def _make_shore_kernel(kernel_path):
    """Make a 3x3 raster with a 9 in the middle and 1s on the outside."""
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
