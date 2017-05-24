"""InVEST XBeach Storm Surge model."""
import os
import logging

import scipy.signal
import scipy.ndimage.filters
import scipy.interpolate
from osgeo import gdal
from osgeo import osr
from osgeo import ogr
import numpy
import rtree
import pygeoprocessing
import shapely.wkb
import shapely.prepared

from . import utils

logging.basicConfig(format='%(asctime)s %(name)-20s %(levelname)-8s \
%(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

LOGGER = logging.getLogger('natcap.invest.xbeach_storm_surge')

_PROFILE_WORK_DIRECTORY = 'shore_profile_work_directory'
_TEMPORARY_FILE_DIRECTORY = 'tmp'

# All files with %s must be replaced by file suffix
_LAND_MASK_FILE_PATTERN = os.path.join(
    _PROFILE_WORK_DIRECTORY, 'land_mask%s.tif')
_SHORE_KERNEL_FILE_PATTERN = os.path.join(
    _TEMPORARY_FILE_DIRECTORY, 'shore_kernel%s.tif')
_SHORE_CONVOLUTION_FILE_PATTERN = os.path.join(
    _PROFILE_WORK_DIRECTORY, 'shore_convolution%s.tif')
_SHORE_MASK_FILE_PATTERN = os.path.join(
    _PROFILE_WORK_DIRECTORY, 'shore_mask%s.tif')
_SHORE_POINTS_FILE_PATTERN = os.path.join(
    _PROFILE_WORK_DIRECTORY, 'shore_points%s.shp')

# Masks are 0 or 1, so 127 is a good value for a signed or unsigned byte
_MASK_NODATA = 127


def execute(args):
    """InVEST XBeach Storm Surge Model.

    Parameters:
        args['workspace_dir'] (string): output directory for intermediate,
            temporary, and final files.
        args['results_suffix'] (string): (optional) string to append to any
            file names produced by this model.
        args['bathymetry_path'] (string): path to a single band bathymetry
            raster that is projected in linear units and values represent
            elevations.
        args['sea_bed_depth'] (float): value in bathymetry raster that
            represents the depth of the sea bed.  In most cases this would be
            0.
        args['nearshore_step_size'] (float): the number of linear units per
            step to sample the profile when in nearshore mode.
        args['offshore_step_size'] (float): the number of linear units per
            step to sample the profile when in offshore mode.
        args['nearshore_depth'] (float): a value in the bathemetry raster
            that indicates where the nearshore 'mode' starts when calculating
            stepsizes; depths less than this have a stepsize of
            `args['nearshore_step_size']`.
        args['offshore_depth'] (float): a value in the bathymetry raster
            indicating where the offshore 'mode' starts; depths greater than
            this have a stepsize value of `args['offshore_step_size']`.
        args['onshore_depth_threshold'] (float): the depth (or elevation) at
            which to cutoff the profile extractor onshore.
        args['offshore_depth_threshold'] (float): the depth at which to
            cutoffthe profile extractor sample length.
        args['max_profile_length'] (float): the maximum cutoff threshold for
            a profile which overrides the `args['*_depth_threshold']` values.
        args['representative_point_vector_path'] (string): Path to a point
            vector file that contains points from which to sample bathymetry.
        args['habitat_vector_directory_path'] (list): Path to a directory
            that contains habitat layers.  The presence of overlap/no overlap
            will be included in the profile results.

    Returns:
        None.
    """
    # Make initial directory structure
    file_suffix = utils.make_suffix_string(args, 'results_suffix')
    profile_work_dir = os.path.join(
        args['workspace_dir'], _PROFILE_WORK_DIRECTORY)
    for dir_path in [args['workspace_dir'], profile_work_dir]:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
    f_reg = {}

    bathymetry_info = pygeoprocessing.get_raster_info(args['bathymetry_path'])
    if bathymetry_info['n_bands'] > 1:
        raise ValueError(
            "The bathymetry raster at '%s' has more than one band, expected "
            "a single band raster.", args['bathymetry_path'])

    sea_bed_depth = float(args['sea_bed_depth'])

    def _land_mask_op(bathymetry):
        """Mask values >= shore height."""
        result = numpy.empty(bathymetry.shape, dtype=numpy.int16)
        result[:] = _MASK_NODATA
        valid_mask = bathymetry != bathymetry_info['nodata'][0]
        result[valid_mask] = numpy.where(
            bathymetry[valid_mask] >= sea_bed_depth, 1, 0)
        return result

    if (
        abs(bathymetry_info['pixel_size'][0]) !=
            abs(bathymetry_info['pixel_size'][1])):
        LOGGER.warn(
            "Bathymetry pixels are not square as %s, this may incorrectly "
            "affect profile direction sampling",
            bathymetry_info['pixel_size'])

    LOGGER.info("Calculating land mask.")
    f_reg['land_mask_path'] = os.path.join(
        args['workspace_dir'], _LAND_MASK_FILE_PATTERN % file_suffix)
    pygeoprocessing.raster_calculator(
        [(args['bathymetry_path'], 1)], _land_mask_op, f_reg['land_mask_path'],
        gdal.GDT_Byte, _MASK_NODATA, calc_raster_stats=False)

    LOGGER.info("Detecting shore pixels.")
    f_reg['shore_kernel_path'] = os.path.join(
        args['workspace_dir'], _SHORE_KERNEL_FILE_PATTERN % file_suffix)
    _make_shore_kernel(f_reg['shore_kernel_path'])
    f_reg['shore_convolution_path'] = os.path.join(
        args['workspace_dir'], _SHORE_CONVOLUTION_FILE_PATTERN % file_suffix)
    pygeoprocessing.convolve_2d(
        (f_reg['land_mask_path'], 1), (f_reg['shore_kernel_path'], 1),
        f_reg['shore_convolution_path'], target_datatype=gdal.GDT_Byte)
    shore_convolution_info = pygeoprocessing.get_raster_info(
        f_reg['shore_convolution_path'])

    def _shore_mask_op(shore_convolution):
        """Mask values on land that border water."""
        result = numpy.empty(shore_convolution.shape, dtype=numpy.int16)
        result[:] = _MASK_NODATA
        valid_mask = shore_convolution != shore_convolution_info['nodata'][0]
        # If a pixel is on land, it gets at least a 9, but if it's all on
        # land it gets an 17 (8 neighboring pixels), so we search between 9
        # and 17 to determine a shore pixel
        result[valid_mask] = numpy.where(
            (shore_convolution[valid_mask] >= 9) &
            (shore_convolution[valid_mask] < 17), 1, _MASK_NODATA)
        return result

    f_reg['shore_mask_path'] = os.path.join(
        args['workspace_dir'], _SHORE_MASK_FILE_PATTERN % file_suffix)
    pygeoprocessing.raster_calculator(
        [(f_reg['shore_convolution_path'], 1)], _shore_mask_op, f_reg['shore_mask_path'],
        gdal.GDT_Byte, _MASK_NODATA, calc_raster_stats=False)

    LOGGER.info("Calculate shore points as a vector.")
    f_reg['shore_points_path'] = os.path.join(
        args['workspace_dir'], _SHORE_POINTS_FILE_PATTERN % file_suffix)
    esri_driver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(f_reg['shore_points_path']):
        os.remove(f_reg['shore_points_path'])
    shore_points_vector = esri_driver.CreateDataSource(f_reg['shore_points_path'])
    shore_points_srs = osr.SpatialReference(
        pygeoprocessing.get_raster_info(
            f_reg['shore_mask_path'])['projection'])
    shore_points_layer = shore_points_vector.CreateLayer(
        'shore_points_path', srs=shore_points_srs, geom_type=ogr.wkbPoint)
    shore_points_layer_defn = shore_points_layer.GetLayerDefn()

    sys.exit()

    # PUT SHORELINE PIXELS INTO R-TREE
    shore_point_index = rtree.index.Index()
    LOGGER.info('Building spatial index for shore points')
    for offset_info, data_block in pygeoprocessing.iterblocks(
            f_reg['shore_mask_path']):
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
            point_feature = ogr.Feature(shore_points_layer_defn)
            shore_point_index.insert(
                0, (x_coord, y_coord), obj=(x_coord, y_coord))
            point_geometry = ogr.Geometry(ogr.wkbPoint)
            point_geometry.AddPoint(x_coord, y_coord)
            point_feature.SetGeometry(point_geometry)
            shore_points_layer.CreateFeature(point_feature)

    LOGGER.info("Constructing offshore profiles")
    representative_point_vector = ogr.Open(
        args['representative_point_vector_path'])
    bathy_rows, bathy_cols = pygeoprocessing.get_row_col_from_uri(
        args['bathymetry_path'])
    for representative_point_layer in representative_point_vector:
        for representative_point in representative_point_layer:
            point_name = representative_point.GetField('name')

            if 'sample_points' not in f_reg:
                f_reg['sample_points'] = {}

            f_reg['sample_points'][point_name] = os.path.join(
                args['workspace_dir'], 'sample_points_%s_%s.shp' % (
                    point_name, args['results_suffix']))
            if os.path.exists(f_reg['sample_points'][point_name]):
                os.remove(f_reg['sample_points'][point_name])
            sample_points_vector = esri_driver.CreateDataSource(
                f_reg['sample_points'][point_name])
            target_sr = osr.SpatialReference(shore_raster.GetProjection())
            sample_points_layer = sample_points_vector.CreateLayer(
                'sample_points_%s' % point_name, srs=target_sr,
                geom_type=ogr.wkbPoint)
            sample_points_layer.CreateField(
                ogr.FieldDefn("s_depth", ogr.OFTReal))
            sample_points_layer.CreateField(
                ogr.FieldDefn("i_depth", ogr.OFTReal))
            sample_points_layer.CreateField(
                ogr.FieldDefn("s_dist", ogr.OFTReal))
            sample_points_layer_defn = sample_points_layer.GetLayerDefn()

            representative_point_geometry = (
                representative_point.GetGeometryRef())
            representative_point = representative_point_geometry.GetPoint(0)
            shore_point = shore_point_index.nearest(
                (representative_point[0], representative_point[1]),
                objects=True).next().object

            vector_length = (
                (representative_point[0]-shore_point[0]) ** 2 +
                (representative_point[1]-shore_point[1]) ** 2) ** 0.5
            x_size = (
                (representative_point[0]-shore_point[0]) / vector_length)
            y_size = (
                (representative_point[1]-shore_point[1]) / vector_length)
            sample_point_list = []

            # test if represenative point is on land or water
            bathymetry_gt = pygeoprocessing.get_geotransform_uri(
                args['bathymetry_path'])
            representative_point_coords = (
                int((representative_point[0] - bathymetry_gt[0]) /
                    bathymetry_gt[1]),
                int((representative_point[1] - bathymetry_gt[3]) /
                    bathymetry_gt[5]))

            # check if representative point is on land, if so; flip direction
            representative_point_elevation = bathymetry_band.ReadAsArray(
                xoff=representative_point_coords[0],
                yoff=representative_point_coords[1],
                win_xsize=1, win_ysize=1)
            if representative_point_elevation >= args['shore_height']:
                x_size *= -1
                y_size *= -1

            # construct step distances away from shore
            offshore_steps = []
            onshore_steps = [0.0]
            current_distance = 0.0
            near_dist = args['step_size'][0][1]
            far_dist = args['step_size'][1][1]
            near_stepsize = args['step_size'][0][0]
            far_stepsize = args['step_size'][1][0]
            current_step = near_stepsize
            while True:
                out_of_range = True
                current_distance += current_step
                if current_distance <= args['offshore_profile_length']:
                    offshore_steps.append(current_distance)
                    out_of_range = False
                if current_distance <= args['onshore_profile_length']:
                    onshore_steps.append(-current_distance)
                    out_of_range = False
                if out_of_range:
                    break
                if current_distance <= near_dist:
                    current_step = near_stepsize
                elif current_distance >= far_dist:
                    current_step = far_stepsize
                else:
                    # bargain basement linear interpolation
                    t_pos = (current_distance - far_dist) / (
                        near_dist - far_dist)
                    current_step = (
                        t_pos * near_stepsize + (1 - t_pos) * far_stepsize)
            # work from offshore inward
            onshore_steps.reverse()

            for step in onshore_steps + offshore_steps:
                point_feature = ogr.Feature(sample_points_layer_defn)
                sample_point_geometry = ogr.Geometry(ogr.wkbPoint)
                sample_point_list.append(
                    (shore_point[0] + x_size * step,
                     shore_point[1] + y_size * step))
                sample_point_geometry.AddPoint(
                    sample_point_list[-1][0], sample_point_list[-1][1])
                point_feature.SetGeometry(sample_point_geometry)
                # make sure we account for step size going negative when we're
                # on land
                point_feature.SetField('s_dist', step)
                sample_points_layer.CreateFeature(point_feature)

            sample_points_layer.SyncToDisk()
            if 'profile_lines' not in f_reg:
                f_reg['profile_lines'] = {}
            f_reg['profile_lines'][point_name] = os.path.join(
                profile_work_dir, 'profile_line_%s_%s.shp' % (
                    args['results_suffix'], point_name))
            if os.path.exists(f_reg['profile_lines'][point_name]):
                os.remove(f_reg['profile_lines'][point_name])
            profile_lines_vector = esri_driver.CreateDataSource(
                f_reg['profile_lines'][point_name])
            target_sr = osr.SpatialReference(shore_raster.GetProjection())
            profile_lines_layer = profile_lines_vector.CreateLayer(
                'profile_lines_%s' % point_name, srs=target_sr,
                geom_type=ogr.wkbLineString)
            profile_lines_layer_defn = profile_lines_layer.GetLayerDefn()
            line_feature = ogr.Feature(profile_lines_layer_defn)
            profile_line_geometry = ogr.Geometry(ogr.wkbLineString)
            profile_line_geometry.AddPoint(
                sample_point_list[0][0], sample_point_list[0][1])
            profile_line_geometry.AddPoint(
                sample_point_list[-1][0], sample_point_list[-1][1])
            line_feature.SetGeometry(profile_line_geometry)
            profile_lines_layer.CreateFeature(line_feature)
            profile_lines_layer.SyncToDisk()
            # extent is xmin, xmax, ymin, ymax
            extent = profile_lines_layer.GetExtent()

            # always want to round up on the max and round down on the min
            # reverse last two because y coord moves up while pixels move down
            extent_in_pixel_coords = (
                int((extent[0] - bathymetry_gt[0]) / bathymetry_gt[1]) - 1,
                int(round(0.5 + (extent[1] - bathymetry_gt[0]) /
                          bathymetry_gt[1])) + 1,
                int(((extent[3] - bathymetry_gt[3]) /
                     bathymetry_gt[5])) - 1,
                int(round(0.5+(extent[2] - bathymetry_gt[3]) /
                          bathymetry_gt[5])) + 1)

            bounding_box_point = ogr.Geometry(ogr.wkbLineString)
            bounding_box_point.AddPoint(
                extent[0], extent[2])
            bounding_box_point.AddPoint(
                extent[0], extent[3])
            bounding_box_point.AddPoint(
                extent[1], extent[3])
            bounding_box_point.AddPoint(
                extent[1], extent[2])
            line_feature.SetGeometry(bounding_box_point)
            profile_lines_layer.CreateFeature(line_feature)

            offset_dict = {
                'xoff': extent_in_pixel_coords[0],
                'yoff': extent_in_pixel_coords[2],
                'win_xsize': (
                    extent_in_pixel_coords[1]-extent_in_pixel_coords[0]),
                'win_ysize': (
                    extent_in_pixel_coords[3]-extent_in_pixel_coords[2]),
            }
            for offset, winsize in [
                    ('xoff', 'win_xsize'), ('yoff', 'win_ysize')]:
                if offset_dict[offset] < 0:
                    # this works because offset will be negative
                    offset_dict[winsize] += offset_dict[offset]
                    LOGGER.error(
                        "Profile sample outside of bathymetry raster, results "
                        "will be clipped to edge.")
                    offset_dict[offset] = 0
            # now check if window is too large
            for offset, winsize, rastersize in [
                    ('xoff', 'win_xsize', bathy_cols),
                    ('yoff', 'win_ysize', bathy_rows)]:
                if offset_dict[offset] + offset_dict[winsize] >= rastersize:
                    # this works because offset will be negative
                    offset_dict[winsize] = rastersize - offset_dict[offset]
                    LOGGER.error(
                        "Profile sample outside of bathymetry raster, results "
                        "will be clipped to edge.")

            profile_lines_layer = None
            profile_lines_vector = None

            clipped_bathymetry_array = bathymetry_band.ReadAsArray(
                **offset_dict)
            x_coordinates = numpy.arange(clipped_bathymetry_array.shape[1])
            x_coordinates = (
                bathymetry_gt[0] + (
                    x_coordinates + offset_dict['xoff'] + 0.5) * bathymetry_gt[1])
            y_coordinates = numpy.arange(clipped_bathymetry_array.shape[0])
            # reverse the y coordinates so they are increasing
            y_coordinates = numpy.flipud(
                bathymetry_gt[3] + (
                    y_coordinates + offset_dict['yoff'] + 0.5) * bathymetry_gt[5])
            # reverse the rows in the array so they match y coordinates
            clipped_bathymetry_array = numpy.flipud(clipped_bathymetry_array)

            interp_fn = scipy.interpolate.RectBivariateSpline(
                y_coordinates, x_coordinates, clipped_bathymetry_array,
                kx=1, ky=1)

            sampled_depth_array = []
            distance_array = []
            for sample_point in sample_points_layer:
                sample_point_geometry = sample_point.GetGeometryRef()
                x_coord = sample_point_geometry.GetX()
                y_coord = sample_point_geometry.GetY()
                step_distance = sample_point.GetField('s_dist')
                distance_array.append(step_distance)
                sampled_depth = interp_fn(y_coord, x_coord)
                sampled_depth_array.append(sampled_depth[0][0])
                sample_point.SetField(
                    's_depth', float(sampled_depth_array[-1]))
                sample_points_layer.SetFeature(sample_point)
            sample_points_layer.ResetReading()
            #smoothed_depth_array = numpy.array(sampled_depth_array)
            sampled_depth_array = numpy.array(sampled_depth_array)
            smoothed_depth_array = scipy.signal.savgol_filter(
                sampled_depth_array,
                (sampled_depth_array.size/4) * 2 - 1, 2)

            if 'profile_table' not in f_reg:
                f_reg['profile_table'] = {}
            f_reg['profile_table'][point_name] = os.path.join(
                args['workspace_dir'], 'profile_table_%s_%s.csv' % (
                    point_name, args['results_suffix']))

            with open(f_reg['profile_table'][point_name], 'w') as profile_table:
                profile_table.write('distance (m),depth (m),smoothed depth (m)')
                habitat_name_list = []
                habitat_geometry_name_list = []
                for path_name_pair in args['habitat_vector_path_list']:
                    habitat_vector = ogr.Open(path_name_pair[0])
                    LOGGER.info("Parsing habitat layer %s", path_name_pair[0])
                    n_layers = habitat_vector.GetLayerCount()
                    for layer_index, habitat_layer in enumerate(
                            habitat_vector):
                        LOGGER.info(
                            "Working on habitat layer %d of %d in %s",
                            layer_index+1, n_layers, path_name_pair[0])
                        n_features = habitat_layer.GetFeatureCount()
                        for feature_index, habitat_feature in enumerate(
                                habitat_layer):
                            LOGGER.info(
                                "Analyzing feature %d of %d", feature_index+1,
                                n_features)
                            habitat_feature_geom = (
                                habitat_feature.GetGeometryRef())
                            # sometimes there's a feature, but no geometry
                            # it's a valid shapefile still.
                            if habitat_feature_geom is None:
                                continue
                            habitat_name = habitat_feature.GetField(
                                path_name_pair[1])
                            if habitat_name not in habitat_name_list:
                                habitat_name_list.append(habitat_name)
                                profile_table.write(',%s' % habitat_name)

                            shapely_geom = shapely.wkb.loads(
                                habitat_feature_geom.ExportToWkb())
                            preped_shapely_geom = (
                                shapely.prepared.prep(shapely_geom))

                            habitat_geometry_name_list.append(
                                (preped_shapely_geom, habitat_name))
                profile_table.write('\n')
                for sample_point, samp_depth, smooth_depth, dist in zip(
                        sample_points_layer, sampled_depth_array,
                        smoothed_depth_array, distance_array):

                    habitat_crossing = [0] * len(habitat_name_list)
                    sample_point_geometry = sample_point.GetGeometryRef()
                    for habitat_geometry, habitat_name in habitat_geometry_name_list:
                        shapely_point = shapely.wkb.loads(
                            sample_point_geometry.ExportToWkb())
                        if habitat_geometry.contains(shapely_point):
                            habitat_crossing[habitat_name_list.index(habitat_name)] = 1

                    x_coord = sample_point_geometry.GetX()
                    y_coord = sample_point_geometry.GetY()
                    profile_table.write(
                        '%f,%f,%f' % (dist, samp_depth, smooth_depth))
                    for crossing_value in habitat_crossing:
                        profile_table.write(',%d' % crossing_value)
                    profile_table.write('\n')
                    sample_point.SetField('i_depth', float(smooth_depth))
                    sample_points_layer.SetFeature(sample_point)
            sample_points_layer.SyncToDisk()
            sample_points_layer = None
            sample_points_vector = None
    # GENERATE SHORELINE PIXELS
    # FOR EACH POINT:
    #   FIND NEAREST SHORELINE POINT
    #   CALCUALTE DIRECTION AS SHORE POINT TO SAMPLE POINT
    #   CREATE A LINE ROOTED AT SHORE PIXEL AS LONG AS REQUEST
    #       OR PRINT AN ERROR IF THE PIXELS GO OUT OF BOUNDS ON THE BATHYMETRY
    #   WALK ALONG LINE FOR EACH STEP:
    #       CALCULATE COORDINATE
    #       SAMPLE RASTER UNDERNEATH
    #       SAMPLE HABITAT LAYER UNDERNEATH


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
