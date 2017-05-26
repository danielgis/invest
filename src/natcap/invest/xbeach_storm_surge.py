"""InVEST XBeach Storm Surge model."""
import shutil
import os
import logging
from datetime import datetime

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

# Replacement pattern for the following is for file_suffix
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

# Replacement pattern for the following is (point_name, file_suffix)
_CLIPPED_BATHYMETRY_FILE_PATTERN = os.path.join(
    _TEMPORARY_FILE_DIRECTORY, 'clipped_bathymetry_%s%s.tif')
_SAMPLE_POINTS_FILE_PATTERN = os.path.join(
    _PROFILE_WORK_DIRECTORY, 'sample_points_%s%s.shp')
_PROFILE_TABLE_FILE_PATTERN = 'profile_table_%s%s.csv'
_X_GRID_FILE_PATTERN = 'x_%s%s.grd'
_Y_GRID_FILE_PATTERN = 'y_%s%s.grd'
_BED_GRID_FILE_PATTERN = 'bed_%s%s.grd'
_VEGGIE_GRID_FILE_PATTERN = 'veggie_%s%s.grd'
_VEGGIEFILE_FILE_PATTERN = 'veggiefile_%s%s.txt'
_XBEACH_WORKSPACE_DIRPATTERN = 'xbeach_workspace%s%s'

_PARAMETERFILE_FILE_PATTERN = 'params.txt'

_REPRESENTATIVE_POINT_ID_FIELDNAME = 'id'
_HABITAT_IDENTIFIER_FIELDNAME = 'hab_type'

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
            cutoff the profile extractor sample length.
        args['max_profile_length'] (float): the maximum length of a profile
            ray that will otherwise override the `args['*_depth_threshold']`
            values.
        args['representative_point_vector_path'] (string): Path to a point
            vector file that contains points from which to sample bathymetry.
        args['habitat_vector_path'] (list): Path to a polygon vector that
            that contains habitat layers.  The presence of overlap/no overlap
            will be included in the profile results.  Must contain the field
            name _HABITAT_IDENTIFIER_FIELDNAME to identify the potential
            habitat IDs for reporting.
        args['habitat_parameter_dir'] (string): path to a directory that
            containts [HABITAT_NAME].txt files corresponding to the habitat
            names provided in args['habitat_vector_path'].  The directory
            contains .txt files matching the values in the 'hab_type' fields
            of the habitat layer that represent biophysical properties of the
            habitat types.  Each file contains attributes for its
            corresponding habitat's morphology including number of vertical
            sections (N), canopy density (Cd), stem diameter (ah), stem
            height, and drag coefficient (bv).
        args['storm_parameter_path'] (string): path to storm parameter file
            for XBeach.  Example contents below:
                Hm0        =     7.0000
                Tp         =     10.000
                mainang    =   270.0000
                gammajsp   =     3.3000
                s          =    20.0000
                fnyq       =     1.0000

    Returns:
        None.
    """
    max_profile_length = float(args['max_profile_length'])
    # TODO: ensure that representative points have an _REPRESENTATIVE_POINT_ID_FIELDNAME field

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
        [(f_reg['shore_convolution_path'], 1)], _shore_mask_op,
        f_reg['shore_mask_path'], gdal.GDT_Byte, _MASK_NODATA,
        calc_raster_stats=False)

    LOGGER.info("Calculate shore points as a vector.")
    f_reg['shore_points_path'] = os.path.join(
        args['workspace_dir'], _SHORE_POINTS_FILE_PATTERN % file_suffix)
    esri_driver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(f_reg['shore_points_path']):
        os.remove(f_reg['shore_points_path'])
    shore_points_vector = esri_driver.CreateDataSource(
        f_reg['shore_points_path'])
    shore_points_srs = osr.SpatialReference(
        pygeoprocessing.get_raster_info(
            f_reg['shore_mask_path'])['projection'])
    shore_points_layer = shore_points_vector.CreateLayer(
        'shore_points_path', srs=shore_points_srs, geom_type=ogr.wkbPoint)
    shore_points_layer_defn = shore_points_layer.GetLayerDefn()

    shore_geotransform = pygeoprocessing.get_raster_info(
        f_reg['shore_mask_path'])['geotransform']
    LOGGER.info('Building spatial index for shore points')
    shore_point_index = rtree.index.Index()
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
    bathy_cols, bathy_rows = bathymetry_info['raster_size']
    representative_point_vector = ogr.Open(
        args['representative_point_vector_path'])
    bathymetry_raster = gdal.Open(args['bathymetry_path'])
    bathymetry_band = bathymetry_raster.GetRasterBand(1)
    for representative_point_layer in representative_point_vector:
        for representative_point in representative_point_layer:
            point_name = representative_point.GetField(
                _REPRESENTATIVE_POINT_ID_FIELDNAME)
            if 'sample_points' not in f_reg:
                f_reg['sample_points'] = {}

            f_reg['sample_points'][point_name] = os.path.join(
                args['workspace_dir'], _SAMPLE_POINTS_FILE_PATTERN % (
                    point_name, file_suffix))
            if os.path.exists(f_reg['sample_points'][point_name]):
                os.remove(f_reg['sample_points'][point_name])
            sample_points_vector = esri_driver.CreateDataSource(
                f_reg['sample_points'][point_name])
            sample_points_layer = sample_points_vector.CreateLayer(
                'sample_points_%s' % point_name, srs=shore_points_srs,
                geom_type=ogr.wkbPoint)
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
            direction_x = (
                (representative_point[0]-shore_point[0]) / vector_length)
            direction_y = (
                (representative_point[1]-shore_point[1]) / vector_length)

            # create a bounding box that will capture the maximum length of
            # the sample vector with form [minx, miny, maxx, maxy]
            profile_bounding_box = [
                representative_point[0] -
                abs(direction_x) * max_profile_length,
                representative_point[1] -
                abs(direction_y) * max_profile_length,
                representative_point[0] +
                abs(direction_x) * max_profile_length,
                representative_point[1] +
                abs(direction_y) * max_profile_length]

            clipped_bathymetry_path = os.path.join(
                args['workspace_dir'], _CLIPPED_BATHYMETRY_FILE_PATTERN % (
                    point_name, file_suffix))
            pygeoprocessing.warp_raster(
                args['bathymetry_path'], bathymetry_info['pixel_size'],
                clipped_bathymetry_path, 'nearest',
                target_bb=profile_bounding_box,
                gtiff_creation_options=['TILED=YES', 'BIGTIFF=IF_SAFER'])
            clipped_bathymetry_info = pygeoprocessing.get_raster_info(
                clipped_bathymetry_path)
            clipped_bathymetry_geotransform = (
                clipped_bathymetry_info['geotransform'])
            clipped_bathymetry_raster = gdal.Open(clipped_bathymetry_path)
            clipped_bathymetry_band = (
                clipped_bathymetry_raster.GetRasterBand(1))

            # test if representative point is on land or water
            representative_point_coords = (
                int((representative_point[0] - shore_geotransform[0]) /
                    shore_geotransform[1]),
                int((representative_point[1] - shore_geotransform[3]) /
                    shore_geotransform[5]))

            # check if representative point is on land by reading one pixel
            # from the raster.  And if so; flip direction of tracing vector
            representative_point_elevation = bathymetry_band.ReadAsArray(
                xoff=representative_point_coords[0],
                yoff=representative_point_coords[1],
                win_xsize=1, win_ysize=1)
            if representative_point_elevation >= sea_bed_depth:
                direction_x *= -1
                direction_y *= -1

            # Start at shore point, walk landward and seaward

            # bounding box coordinates are this [minx, miny, maxx, maxy]
            # initialize to crazy bounding box so test for bounds fails
            clipped_bathymetry_block = None
            sample_points = []
            first_point_added = False
            for direction_scale in [-1, 1]:
                window_bounding_box = [1, 1, -1, -1]
                current_point = shore_point
                distance_travelled = 0.0
                while True:
                    # test if point outside clipped bathymetry bounding box
                    # and if so load a
                    # new bounding box
                    if (current_point[0] < window_bounding_box[0] or
                            current_point[0] > window_bounding_box[2] or
                            current_point[1] < window_bounding_box[1] or
                            current_point[1] > window_bounding_box[3]):
                        # convert coordinates to raster indexes
                        current_point_coords = (
                            int((current_point[0] -
                                 clipped_bathymetry_geotransform[0]) /
                                clipped_bathymetry_geotransform[1]),
                            int((current_point[1] -
                                 clipped_bathymetry_geotransform[3]) /
                                clipped_bathymetry_geotransform[5]))
                        block_index = [
                            current_point_coords[0] /
                            clipped_bathymetry_info['block_size'][0],
                            current_point_coords[1] /
                            clipped_bathymetry_info['block_size'][1]]
                        # make raster pixel bounds to get an efficient
                        # ReadAsArray call, format: [minx, miny, maxx, maxy]
                        raster_window_bounds = [
                            block_index[0] *
                            clipped_bathymetry_info['block_size'][0],
                            block_index[1] *
                            clipped_bathymetry_info['block_size'][1],
                            (block_index[0] + 1) *
                            clipped_bathymetry_info['block_size'][0],
                            (block_index[1] + 1) *
                            clipped_bathymetry_info['block_size'][1],
                        ]
                        # update bounds by 1 pixel if possible
                        if raster_window_bounds[0] > 0:
                            raster_window_bounds[0] -= 1
                        if raster_window_bounds[1] > 0:
                            raster_window_bounds[1] -= 1
                        if (raster_window_bounds[2] <
                                clipped_bathymetry_info['raster_size'][0] - 1):
                            raster_window_bounds[2] += 1
                        elif (raster_window_bounds[2] >
                              clipped_bathymetry_info['raster_size'][0]):
                            # fix a case where a boundary might be bigger than
                            # the raster
                            raster_window_bounds[2] = (
                                clipped_bathymetry_info['raster_size'][0])
                        if (raster_window_bounds[3] <
                                clipped_bathymetry_info['raster_size'][1] - 1):
                            raster_window_bounds[3] += 1
                        elif (raster_window_bounds[3] >
                              clipped_bathymetry_info['raster_size'][1]):
                            # fix a case where a boundary might be bigger than
                            # the raster
                            raster_window_bounds[3] = (
                                clipped_bathymetry_info['raster_size'][1])

                        # update the bounding box in projected coordinates
                        # format: [minx, miny, maxx, maxy]

                        window_bounding_box = [
                            clipped_bathymetry_geotransform[0] +
                            raster_window_bounds[0] *
                            clipped_bathymetry_geotransform[1],
                            clipped_bathymetry_geotransform[3] +
                            raster_window_bounds[1] *
                            clipped_bathymetry_geotransform[5],
                            clipped_bathymetry_geotransform[0] +
                            raster_window_bounds[2] *
                            clipped_bathymetry_geotransform[1],
                            clipped_bathymetry_geotransform[3] +
                            raster_window_bounds[3] *
                            clipped_bathymetry_geotransform[5],
                            ]

                        clipped_bathymetry_block = (
                            clipped_bathymetry_band.ReadAsArray(
                                xoff=raster_window_bounds[0],
                                yoff=raster_window_bounds[1],
                                win_xsize=(
                                    raster_window_bounds[2] -
                                    raster_window_bounds[0]),
                                win_ysize=(
                                    raster_window_bounds[3] -
                                    raster_window_bounds[1])))

                        x_coordinates = numpy.arange(
                            clipped_bathymetry_block.shape[1])
                        x_coordinates = (
                            clipped_bathymetry_geotransform[0] + (
                                x_coordinates + raster_window_bounds[0] + 0.5) *
                            clipped_bathymetry_geotransform[1])
                        y_coordinates = numpy.arange(
                            clipped_bathymetry_block.shape[0])
                        # reverse the y coordinates so they are increasing
                        y_coordinates = numpy.flipud(
                            clipped_bathymetry_geotransform[3] + (
                                y_coordinates + raster_window_bounds[1] + 0.5) *
                            clipped_bathymetry_geotransform[5])
                        # reverse rows in the array so they match y coordinates
                        clipped_bathymetry_block = numpy.flipud(
                            clipped_bathymetry_block)

                        # create a linear interpolator for the current block
                        interp_fn = scipy.interpolate.RectBivariateSpline(
                            y_coordinates, x_coordinates,
                            clipped_bathymetry_block, kx=1, ky=1)

                        raster_window_bounds = None

                    # append the current point and its interpolated depth
                    interpolated_depth = interp_fn(
                        current_point[1], current_point[0])[0, 0]

                    if distance_travelled > 0 or not first_point_added:
                        sample_points.append(
                            (distance_travelled * direction_scale,
                             current_point, interpolated_depth))
                        first_point_added = True
                    if (
                            interpolated_depth >=
                            float(args['onshore_depth_threshold']) or
                            interpolated_depth <=
                            float(args['offshore_depth_threshold'])):
                        # we're too high or too low, okay to quit
                        break
                    step_size = None
                    if interpolated_depth >= float(args['nearshore_depth']):
                        step_size = float(args['nearshore_step_size'])
                    elif interpolated_depth <= float(args['offshore_depth']):
                        step_size = float(args['offshore_step_size'])
                    else:
                        interp_parameter = (
                            interpolated_depth -
                            float(args['nearshore_depth'])) / (
                            float(args['offshore_depth']) -
                            float(args['nearshore_depth']))
                        step_size = (
                            float(args['nearshore_step_size']) * (
                                1 - interp_parameter) +
                            float(args['offshore_step_size']) * (
                                interp_parameter))

                    # udpate to next sample step
                    current_point = (
                        current_point[0] +
                        direction_x * step_size * direction_scale,
                        current_point[1] +
                        direction_y * step_size * direction_scale)
                    distance_travelled += step_size
                    if distance_travelled > float(args['max_profile_length']):
                        # we're done if we travelled too far
                        break

            clipped_bathymetry_block = None
            clipped_bathymetry_band = None
            clipped_bathymetry_raster = None

            for distance, (point_x, point_y), depth in sample_points:
                point_feature = ogr.Feature(sample_points_layer_defn)
                point_geometry = ogr.Geometry(ogr.wkbPoint)
                point_geometry.AddPoint(point_x, point_y)
                point_feature.SetGeometry(point_geometry)
                point_feature.SetField('i_depth', float(depth))
                point_feature.SetField('s_dist', float(distance))
                sample_points_layer.CreateFeature(point_feature)

            sample_points_vector.SyncToDisk()
            sample_points_layer = None
            sample_points_vector = None

            if 'profile_table' not in f_reg:
                f_reg['profile_table'] = {}
            f_reg['profile_table'][point_name] = os.path.join(
                args['workspace_dir'], _PROFILE_TABLE_FILE_PATTERN % (
                    point_name, file_suffix))

            with open(f_reg['profile_table'][point_name],
                      'w') as profile_table:
                profile_table.write(
                    'distance (m),depth (m)')
                habitat_name_list = []
                habitat_geometry_name_list = []
                habitat_vector = ogr.Open(args['habitat_vector_path'])
                LOGGER.info(
                    "Parsing habitat layer %s", args['habitat_vector_path'])
                n_layers = habitat_vector.GetLayerCount()
                for layer_index, habitat_layer in enumerate(
                        habitat_vector):
                    LOGGER.info(
                        "Working on habitat layer %d of %d in %s",
                        layer_index+1, n_layers, args['habitat_vector_path'])
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
                            _HABITAT_IDENTIFIER_FIELDNAME)
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

                for dist, (point_x, point_y), samp_depth in sorted(
                        sample_points):
                    habitat_crossing = [0] * len(habitat_name_list)
                    point_geometry = ogr.Geometry(ogr.wkbPoint)
                    point_geometry.AddPoint(point_x, point_y)
                    for habitat_geometry, habitat_name in \
                            habitat_geometry_name_list:
                        shapely_point = shapely.wkb.loads(
                            point_geometry.ExportToWkb())
                        if habitat_geometry.contains(shapely_point):
                            habitat_crossing[
                                habitat_name_list.index(habitat_name)] = 1

                    profile_table.write(
                        '%f,%f' % (dist, samp_depth))
                    for crossing_value in habitat_crossing:
                        profile_table.write(',%d' % crossing_value)
                    profile_table.write('\n')

            LOGGER.info(
                "Create XBeach parameter run files for %s", point_name)
            xbeach_workspace_path = os.path.join(
                args['workspace_dir'], _XBEACH_WORKSPACE_DIRPATTERN % (
                    point_name, file_suffix))
            if not os.path.exists(xbeach_workspace_path):
                os.makedirs(xbeach_workspace_path)

            # write out the file that lists the habitats in the order they're
            # reported in.
            veggiefile_path = os.path.join(
                xbeach_workspace_path, _VEGGIEFILE_FILE_PATTERN % (
                    point_name, file_suffix))
            veggiefile_file = open(veggiefile_path, 'w')
            for habitat in habitat_name_list:
                base_habitat_path = os.path.join(
                    args['habitat_parameter_dir'], '%s.txt' % habitat)
                target_habitat_path = os.path.join(
                    xbeach_workspace_path,
                    os.path.basename(base_habitat_path))
                LOGGER.debug("%s, %s" % (base_habitat_path, target_habitat_path))
                shutil.copyfile(base_habitat_path, target_habitat_path)
                veggiefile_file.write('%s.txt\n' % habitat)

            with open(f_reg['profile_table'][point_name],
                      'r') as profile_table:
                profile_table.readline()  # toss the first line
                x_grid_path = os.path.join(
                    xbeach_workspace_path, _X_GRID_FILE_PATTERN % (
                        point_name, file_suffix))
                y_grid_path = os.path.join(
                    xbeach_workspace_path, _Y_GRID_FILE_PATTERN % (
                        point_name, file_suffix))
                bed_grid_path = os.path.join(
                    xbeach_workspace_path, _BED_GRID_FILE_PATTERN % (
                        point_name, file_suffix))
                veggiemap_path = os.path.join(
                    xbeach_workspace_path, _VEGGIE_GRID_FILE_PATTERN % (
                        point_name, file_suffix))
                x_grid_file = open(x_grid_path, 'w')
                y_grid_file = open(y_grid_path, 'w')
                bed_grid_file = open(bed_grid_path, 'w')
                veggiemap_file = open(veggiemap_path, 'w')
                for line in profile_table:
                    line_values = [
                        float(x) for x in line.split(',')]
                    # format is x, depth, *habitat_overlap
                    # index 0 is the distance along the ray
                    x_grid_file.write('%f,' % line_values[0])
                    y_grid_file.write('0,')
                    # index 1 is the height of the bathymetry
                    bed_grid_file.write('%f,' % line_values[1])
                    try:
                        # write the first index we see, we can only write one
                        veggiemap_file.write('%d,' % (
                            line_values[2:].index(1) + 1))
                    except ValueError:
                        # no habitat in this sample
                        veggiemap_file.write('0')
                for file in [
                        x_grid_file, y_grid_file, bed_grid_file,
                        veggiemap_file]:
                    file.write('NaN')
                    file.close()
            parameter_file_path = os.path.join(
                xbeach_workspace_path, _PARAMETERFILE_FILE_PATTERN)
            n_points = len(sample_points)
            xbeach_storm_parameter_path = os.path.join(
                xbeach_workspace_path, os.path.basename(
                    args['storm_parameter_path']))
            shutil.copyfile(
                args['storm_parameter_path'], xbeach_storm_parameter_path)
            _write_xbeach_parameter_file(
                parameter_file_path, n_points,
                os.path.basename(bed_grid_path),
                os.path.basename(x_grid_path),
                os.path.basename(y_grid_path),
                os.path.basename(xbeach_storm_parameter_path),
                os.path.basename(veggiefile_path),
                os.path.basename(veggiemap_path))


def _write_xbeach_parameter_file(
        parameter_file_path, n_points, bed_grid_path, x_grid_path,
        y_grid_path, storm_parameter_path, veggiefile_path, veggie_grid_path):
    """Write the XBeach parameter file as described in the design doc.

    This is the example given in the design doc: https://drive.google.com/file/d/0B--GW7O9bzP9TWlmV3BxTGVHalU/view

    Parameters:
        parameter_file_path (string): path to desired output parameter file.

    Returns:
        None
    """
    param_file = open(parameter_file_path, 'w')
    param_file.write("%% XBeach parameter settings input file\n")
    param_file.write("%% date: %s\n" % str(datetime.now()))
    param_file.write("\n")
    param_file.write("%% Grid parameters \n")
    param_file.write("depfile   = %s\n" % bed_grid_path)
    param_file.write("posdwn    = 0\n")
    param_file.write("nx        = %d\n" % n_points)
    param_file.write("ny        = 0\n")
    param_file.write("alfa      = 0\n")
    param_file.write("vardx     = 1\n")
    param_file.write("xfile     = %s\n" % x_grid_path)
    param_file.write("yfile     = %s\n" % y_grid_path)
    param_file.write("xori      = 0\n")
    param_file.write("yori      = 0\n")
    param_file.write("thetamin   225\n")
    param_file.write("thetamax   315\n")
    param_file.write("dtheta     90\n")
    param_file.write("thetanaut  1\n")
    param_file.write("\n")
    param_file.write("%% Initial conditions \n")
    param_file.write("zs0       = 1.0\n")
    param_file.write("\n")
    param_file.write("%% Model time \n")
    param_file.write("tstop     = 21600\n")
    param_file.write("tintg     = 10\n")
    param_file.write("\n")
    param_file.write("%% Wave boundary condition parameters \n")
    param_file.write("instat    = jons\n")
    param_file.write("\n")
    param_file.write("%% Wave-spectrum boundary condition parameters \n")
    param_file.write("bcfile    = %s\n" % storm_parameter_path)
    param_file.write("rt        = 1800\n")
    param_file.write("dtbc      = 1\n")
    param_file.write("\n")
    param_file.write("%% Output variables \n")
    param_file.write("\n")
    param_file.write("%% Vegetation (added by SMV)\n")
    param_file.write("\n")
    param_file.write("vegetation = 1\n")
    param_file.write("veggiefile = %s\n" % veggiefile_path)
    param_file.write("veggiemapfile = %s\n" % veggie_grid_path)
    param_file.write("morphology = 0\n")
    param_file.write("nonh = 1\n")
    param_file.write("swave = 0\n")
    param_file.write("wind = 1\n")
    param_file.write("windv = 30.0\n")
    param_file.write("windth = 270\n")
    param_file.write("nglobalvar = 3\n")
    param_file.write("zs\n")
    param_file.write("zb\n")
    param_file.write("hh\n")
    param_file.close()


def _make_shore_kernel(kernel_path):
    """Make a 3x3 raster with a 9 in the middle and 1s on the outside."""
    driver = gdal.GetDriverByName('GTiff')
    kernel_raster = driver.Create(
        kernel_path.encode('utf-8'), 3, 3, 1,
        gdal.GDT_Byte)

    # Make some kind of geotransform, it doesn't matter what but
    # will make GIS libraries behave better if it's all defined
    kernel_raster.SetGeoTransform([0, 1, 0, 0, 0, -1])
    srs = osr.SpatialReference()
    srs.SetWellKnownGeogCS('WGS84')
    kernel_raster.SetProjection(srs.ExportToWkt())

    kernel_band = kernel_raster.GetRasterBand(1)
    kernel_band.SetNoDataValue(127)
    kernel_band.WriteArray(numpy.array([[1, 1, 1], [1, 9, 1], [1, 1, 1]]))
