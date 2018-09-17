"""InVEST Carbon Edge Effect Model.

An implementation of the model described in 'Degradation in carbon stocks
near tropical forest edges', by Chaplin-Kramer et. al (in review).
"""
from __future__ import absolute_import
import os
import logging
import tempfile
import time
import uuid

from . import utils
import numpy
from osgeo import gdal
from osgeo import ogr
import pygeoprocessing
import scipy.spatial

from . import validation

LOGGER = logging.getLogger('natcap.invest.carbon_edge_effect')

# grid cells are 100km.  becky says 500km is a good upper bound to search
DISTANCE_UPPER_BOUND = 500e3

# helpful to have a global nodata defined for all the carbon map rasters
CARBON_MAP_NODATA = -9999


def execute(args):
    """Forest Carbon Edge Effect.

    InVEST Carbon Edge Model calculates the carbon due to edge effects in
    tropical forest pixels.

    Parameters:
        args['workspace_dir'] (string): a uri to the directory that will write
            output and other temporary files during calculation. (required)
        args['results_suffix'] (string): a string to append to any output file
            name (optional)
        args['n_nearest_model_points'] (int): number of nearest neighbor model
            points to search for
        args['aoi_vector_path'] (string): (optional) if present, a path to a
            shapefile that will be used to aggregate carbon stock results at
            the end of the run.
        args['biophysical_table_uri'] (string): a path to a CSV table that has
            at least the fields 'lucode' and 'c_above'. If
            ``args['compute_forest_edge_effects'] == True``, table must
            also contain an 'is_tropical_forest' field.  If
            ``args['pools_to_calculate'] == 'all'``, this table must contain
            the fields 'c_below', 'c_dead', and 'c_soil'.

                * ``lucode``: an integer that corresponds to landcover codes in
                  the raster ``args['lulc_raster_path']``
                * ``is_tropical_forest``: either 0 or 1 indicating whether the
                  landcover type is forest (1) or not (0).  If 1, the value
                  in ``c_above`` is ignored and instead calculated from the
                  edge regression model.
                * ``c_above``: floating point number indicating tons of above
                  ground carbon per hectare for that landcover type
                * ``{'c_below', 'c_dead', 'c_soil'}``: three other optional
                  carbon pools that will statically map landcover types to the
                  carbon densities in the table.

                Example::

                    lucode,is_tropical_forest,c_above,c_soil,c_dead,c_below
                    0,0,32.8,5,5.2,2.1
                    1,1,n/a,2.5,0.0,0.0
                    2,1,n/a,1.8,1.0,0.0
                    16,0,28.1,4.3,0.0,2.0

                Note the "n/a" in ``c_above`` are optional since that field
                is ignored when ``is_tropical_forest==1``.
        args['lulc_raster_path'] (string): path to a integer landcover code raster
        args['pools_to_calculate'] (string): if "all" then all carbon pools
            will be calculted.  If any other value only above ground carbon
            pools will be calculated and expect only a 'c_above' header in
            the biophysical table. If "all" model expects 'c_above',
            'c_below', 'c_dead', 'c_soil' in header of biophysical_table and
            will make a translated carbon map for each based off the landcover
            map.
        args['compute_forest_edge_effects'] (boolean): if True, requires
            biophysical table to have 'is_tropical_forest' forest field, and
            any landcover codes that have a 1 in this column calculate carbon
            stocks using the Chaplin-Kramer et. al method and ignore 'c_above'.
        args['tropical_forest_edge_carbon_model_shape_uri'] (string): path to a
            shapefile that defines the regions for the local carbon edge
            models.  Has at least the fields 'method', 'theta1', 'theta2',
            'theta3'.  Where 'method' is an int between 1..3 describing the
            biomass regression model, and the thetas are floating point numbers
            that have different meanings depending on the 'method' parameter.
            Specifically,

                * method 1 (asymptotic model)::

                    biomass = theta1 - theta2 * exp(-theta3 * edge_dist_km)

                * method 2 (logarithmic model)::

                    # NOTE: theta3 is ignored for this method
                    biomass = theta1 + theta2 * numpy.log(edge_dist_km)

                * method 3 (linear regression)::

                    biomass = theta1 + theta2 * edge_dist_km

        args['biomass_to_carbon_conversion_factor'] (string/float): Number by
            which to multiply forest biomass to convert to carbon in the edge
            effect calculation.

    Returns:
        None
    """
    # just check that the AOI exists since it wouldn't crash until the end of
    # the whole model run if it didn't.
    if 'aoi_vector_path' in args and args['aoi_vector_path'] != '':
        aoi_vector = gdal.OpenEx(args['aoi_vector_path'], gdal.OF_VECTOR)
        if not aoi_vector:
            raise ValueError("Unable to open aoi at: %s" % args['aoi_vector_path'])
        aoi_vector = None

    output_dir = args['workspace_dir']
    intermediate_dir = os.path.join(
        args['workspace_dir'], 'intermediate_outputs')
    utils.make_directories(
        [output_dir, intermediate_dir])
    file_suffix = utils.make_suffix_string(args, 'results_suffix')

    # used to keep track of files generated by this module
    output_file_registry = {
        'c_above_map': os.path.join(
            intermediate_dir, 'c_above_carbon_stocks%s.tif' % file_suffix),
        'carbon_map': os.path.join(
            output_dir, 'carbon_map%s.tif' % file_suffix),
        'aggregated_result_vector': os.path.join(
            output_dir, 'aggregated_carbon_stocks%s.shp' % file_suffix)
    }

    if args['pools_to_calculate'] == 'all':
        output_file_registry['c_below_map'] = os.path.join(
            intermediate_dir, 'c_below_carbon_stocks%s.tif' % file_suffix)
        output_file_registry['c_soil_map'] = os.path.join(
            intermediate_dir, 'c_soil_carbon_stocks%s.tif' % file_suffix)
        output_file_registry['c_dead_map'] = os.path.join(
            intermediate_dir, 'c_dead_carbon_stocks%s.tif' % file_suffix)

    if args['compute_forest_edge_effects']:
        output_file_registry['edge_distance'] = os.path.join(
            intermediate_dir, 'edge_distance%s.tif' % file_suffix)
        output_file_registry['tropical_forest_edge_carbon_map'] = os.path.join(
            intermediate_dir, 'tropical_forest_edge_carbon_stocks%s.tif' %
            file_suffix)

    # Map non-forest landcover codes to carbon biomasses
    LOGGER.info('calculating direct mapped carbon stocks')
    carbon_maps = []
    biophysical_table = utils.build_lookup_from_csv(
        args['biophysical_table_uri'], 'lucode', to_lower=False)
    biophysical_keys = [
        x.lower() for x in biophysical_table.itervalues().next().keys()]
    pool_list = [('c_above', True)]
    if args['pools_to_calculate'] == 'all':
        pool_list.extend([
            ('c_below', False), ('c_soil', False), ('c_dead', False)])
    for carbon_pool_type, ignore_tropical_type in pool_list:
        if carbon_pool_type in biophysical_keys:
            carbon_maps.append(
                output_file_registry[carbon_pool_type+'_map'])
            _calculate_lulc_carbon_map(
                args['lulc_raster_path'], args['biophysical_table_uri'],
                carbon_pool_type, ignore_tropical_type,
                args['compute_forest_edge_effects'], carbon_maps[-1])

    if args['compute_forest_edge_effects']:
        # generate a map of pixel distance to forest edge from the landcover map
        LOGGER.info('calculating distance from forest edge')
        _map_distance_from_tropical_forest_edge(
            args['lulc_raster_path'], args['biophysical_table_uri'],
            output_file_registry['edge_distance'])

        # Build spatial index for gridded global model for closest 3 points
        LOGGER.info('Building spatial index for forest edge models.')
        kd_tree, theta_model_parameters, method_model_parameter = (
            _build_spatial_index(
                args['lulc_raster_path'], intermediate_dir,
                args['tropical_forest_edge_carbon_model_shape_uri']))

        # calculate the edge carbon effect on forests
        LOGGER.info('calculating forest edge carbon')
        _calculate_tropical_forest_edge_carbon_map(
            output_file_registry['edge_distance'], kd_tree,
            theta_model_parameters, method_model_parameter,
            int(args['n_nearest_model_points']),
            float(args['biomass_to_carbon_conversion_factor']),
            output_file_registry['tropical_forest_edge_carbon_map'])

        # This is also a carbon stock
        carbon_maps.append(
            output_file_registry['tropical_forest_edge_carbon_map'])

    # combine maps into a single output
    LOGGER.info('combining carbon maps into single raster')

    def combine_carbon_maps(*carbon_maps):
        """This combines the carbon maps into one and leaves nodata where all
        the inputs maps were nodata"""
        result = numpy.zeros(carbon_maps[0].shape)
        nodata_mask = numpy.empty(carbon_maps[0].shape, dtype=numpy.bool)
        nodata_mask[:] = True
        for carbon_map in carbon_maps:
            valid_mask = carbon_map != CARBON_MAP_NODATA
            nodata_mask &= ~valid_mask
            result[valid_mask] += carbon_map[valid_mask]
        result[nodata_mask] = CARBON_MAP_NODATA
        return result

    # aligned_raster_list = [
    #     os.path.join(intermediate_dir, os.path.basename(path).replace(
    #         '.tif', '_aligned.tif')) for path in carbon_maps]

    # pygeoprocessing.align_and_resize_raster_stack(
    #     carbon_maps, aligned_raster_list,
    #     ['near']*len(carbon_maps), cell_size_in_meters, 'intersection')

    # carbon_maps = [
    #     os.path.join(intermediate_dir, os.path.basename(path).replace(
    #         '.tif', '_aligned.tif')) for path in carbon_maps]

    carbon_maps_band_list = [(path, 1) for path in carbon_maps]

    pygeoprocessing.raster_calculator(
        carbon_maps_band_list, combine_carbon_maps,
        output_file_registry['carbon_map'], gdal.GDT_Float32, CARBON_MAP_NODATA
        )

    # generate report (optional) by aoi if they exist
    if 'aoi_vector_path' in args and args['aoi_vector_path'] != '':
        LOGGER.info('aggregating carbon map by aoi')
        _aggregate_carbon_map(
            args['aoi_vector_path'], output_file_registry['carbon_map'],
            output_file_registry['aggregated_result_vector'])


def _aggregate_carbon_map(
        aoi_vector_path, carbon_map_path, target_aggregate_vector_path):
    """Helper function to aggregate carbon values for the given serviceshed.
    Generates a new shapefile that's a copy of 'aoi_vector_path' in
    'workspace_dir' with mean and sum values from the raster at
    'carbon_map_path'

    Parameters:
        aoi_vector_path (string): path to shapefile that will be used to
            aggregate raster at'carbon_map_path'.
        workspace_dir (string): path to a directory that function can copy
            the shapefile at aoi_vector_path into.
        carbon_map_path (string): path to raster that will be aggregated by
            the given serviceshed polygons
        target_aggregate_vector_path (string): path to an ESRI shapefile that
            will be created by this function as the aggregating output.

    Returns:
        None
    """
    aoi_vector = gdal.OpenEx(aoi_vector_path, gdal.OF_VECTOR)
    driver = gdal.GetDriverByName('ESRI Shapefile')

    if os.path.exists(target_aggregate_vector_path):
        os.remove(target_aggregate_vector_path)
    target_aggregate_vector = driver.CreateCopy(
        target_aggregate_vector_path, aoi_vector)
    aoi_vector = None
    target_aggregate_layer = target_aggregate_vector.GetLayer()

    # make an identifying id per polygon that can be used for aggregation
    while True:
        serviceshed_defn = target_aggregate_layer.GetLayerDefn()
        poly_id_field = str(uuid.uuid4())[-8:]
        if serviceshed_defn.GetFieldIndex(poly_id_field) == -1:
            break
    layer_id_field = ogr.FieldDefn(poly_id_field, ogr.OFTInteger)
    target_aggregate_layer.CreateField(layer_id_field)
    target_aggregate_layer.StartTransaction()
    for poly_index, poly_feat in enumerate(target_aggregate_layer):
        poly_feat.SetField(poly_id_field, poly_index)
        target_aggregate_layer.SetFeature(poly_feat)
    target_aggregate_layer.CommitTransaction()
    target_aggregate_layer.SyncToDisk()

    # aggregate carbon stocks by the new ID field
    try:
        serviceshed_stats = pygeoprocessing.zonal_statistics(
            (carbon_map_path, 1), target_aggregate_vector_path, poly_id_field)
    except ValueError:
        raise ValueError("There is no intersection between the land cover map "
                         "and service areas of interest.")

    # don't need a random poly id anymore
    target_aggregate_layer.DeleteField(
        serviceshed_defn.GetFieldIndex(poly_id_field))

    carbon_sum_field = ogr.FieldDefn('c_sum', ogr.OFTReal)
    carbon_sum_field.SetWidth(24)
    carbon_sum_field.SetPrecision(11)
    carbon_mean_field = ogr.FieldDefn('c_ha_mean', ogr.OFTReal)
    carbon_mean_field.SetWidth(24)
    carbon_mean_field.SetPrecision(11)

    target_aggregate_layer.CreateField(carbon_sum_field)
    target_aggregate_layer.CreateField(carbon_mean_field)

    target_aggregate_layer.ResetReading()
    target_aggregate_layer.StartTransaction()

    for poly_index, poly_feat in enumerate(target_aggregate_layer):
        poly_feat.SetField(
            'c_sum', serviceshed_stats[poly_index]['sum'])
        # calculates mean pixel value per ha in for each feature in AOI
        poly_geom = poly_feat.GetGeometryRef()
        poly_area_ha = poly_geom.GetArea() / 1e4  # converts m^2 to hectare
        poly_geom = None
        poly_feat.SetField(
            'c_ha_mean', serviceshed_stats[poly_index]['sum']/poly_area_ha)

        target_aggregate_layer.SetFeature(poly_feat)
    target_aggregate_layer.CommitTransaction()


def _calculate_lulc_carbon_map(
        lulc_raster_path, biophysical_table_uri, carbon_pool_type,
        ignore_tropical_type, compute_forest_edge_effects, carbon_map_path):
    """Calculates the carbon on the map based on non-forest landcover types
    only.

    Parameters:
        lulc_raster_path (string): a filepath to the landcover map that contains
            integer landcover codes
        biophysical_table_uri (string): a filepath to a csv table that indexes
            landcover codes to surface carbon, contains at least the fields
            'lucode' (landcover integer code), 'is_tropical_forest' (0 or 1
            depending on landcover code type), and 'c_above' (carbon density in
            terms of Mg/Ha)
        carbon_pool_type (string): a carbon mapping field in
            biophysical_table_uri.  ex. 'c_above', 'c_below', ...
        ignore_tropical_type (boolean): if true, any landcover type whose
            'is_tropical_forest' field == 1 will be ignored for mapping the
            carbon pool type.
        compute_forest_edge_effects (boolean): if true the 'is_tropical_forest'
            header will be considered, if not, it is ignored
        carbon_map_path (string): a filepath to the output raster
            that will contain total mapped carbon per cell.

    Returns:
        None"""

    # classify forest pixels from lulc
    biophysical_table = utils.build_lookup_from_csv(
        biophysical_table_uri, 'lucode', to_lower=False)

    lucode_to_per_pixel_carbon = {}
    cell_area_ha = pygeoprocessing.get_raster_info(
        lulc_raster_path)['mean_pixel_size'] ** 2 / 10000.0

    # Build a lookup table
    for lucode in biophysical_table:
        if compute_forest_edge_effects:
            is_tropical_forest = (
                int(biophysical_table[int(lucode)]['is_tropical_forest']))
        else:
            is_tropical_forest = 0
        if ignore_tropical_type and is_tropical_forest == 1:
            # if forest, lookup table is nodata
            lucode_to_per_pixel_carbon[int(lucode)] = CARBON_MAP_NODATA
        else:
            try:
                lucode_to_per_pixel_carbon[int(lucode)] = float(
                    biophysical_table[lucode][carbon_pool_type]) * cell_area_ha
            except ValueError:
                raise ValueError(
                    "Could not interpret carbon pool value as a number. "
                    "lucode: %s, pool_type: %s, value: %s" %
                    (lucode, carbon_pool_type,
                     biophysical_table[lucode][carbon_pool_type]))

    # map aboveground carbon from table to lulc that is not forest
    pygeoprocessing.reclassify_raster(
        (lulc_raster_path, 1), lucode_to_per_pixel_carbon,
        carbon_map_path, gdal.GDT_Float32, CARBON_MAP_NODATA)


def _map_distance_from_tropical_forest_edge(
        lulc_raster_path, biophysical_table_uri, edge_distance_path):
    """Generates a raster of forest edge distances where each pixel is the
    distance to the edge of the forest in meters.

    Parameters:
        lulc_raster_path (string): path to the landcover raster that contains integer
            landcover codes
        biophysical_table_uri (string): a path to a csv table that indexes
            landcover codes to forest type, contains at least the fields
            'lucode' (landcover integer code) and 'is_tropical_forest' (0 or 1
            depending on landcover code type)
        edge_distance_path (string): path to output raster where each pixel
            contains the euclidian pixel distance to nearest forest edges on
            all non-nodata values of lulc_raster_path

    Returns:
        None"""

    # Build a list of forest lucodes
    biophysical_table = utils.build_lookup_from_csv(
        biophysical_table_uri, 'lucode', to_lower=False)
    forest_codes = [
        lucode for (lucode, ludata) in biophysical_table.iteritems()
        if int(ludata['is_tropical_forest']) == 1]

    # Make a raster where 1 is non-forest landcover types and 0 is forest
    forest_mask_nodata = 255
    lulc_nodata = pygeoprocessing.get_raster_info(lulc_raster_path)['nodata']

    def mask_non_forest_op(lulc_array):
        """converts forest lulc codes to 1"""
        non_forest_mask = ~numpy.in1d(
            lulc_array.flatten(), forest_codes).reshape(lulc_array.shape)
        nodata_mask = lulc_array == lulc_nodata
        return numpy.where(nodata_mask, forest_mask_nodata, non_forest_mask)
    non_forest_mask_path_handle, non_forest_mask_path = tempfile.mkstemp()
    os.close(non_forest_mask_path_handle)
    pygeoprocessing.raster_calculator(
        [(lulc_raster_path, 1)], mask_non_forest_op, non_forest_mask_path,
        gdal.GDT_Byte, forest_mask_nodata)

    # Do the distance transform on non-forest pixels
    pygeoprocessing.distance_transform_edt(
        (non_forest_mask_path, 1), edge_distance_path)

    # good practice to delete temporary files when we're done with them
    os.remove(non_forest_mask_path)


def _build_spatial_index(
        base_raster_path, local_model_dir,
        tropical_forest_edge_carbon_model_shapefile_path):
    """Build a kd-tree index of the locally projected globally georeferenced
    carbon edge model parameters.

    Parameters:
        base_raster_path (string): path to a raster that is used to define the
            bounding box and projection of the local model.
        local_model_dir (string): path to a directory where we can write a
            shapefile of the locally projected global data model grid.
            Function will create a file called 'local_carbon_shape.shp' in
            that location and overwrite one if it exists.
        tropical_forest_edge_carbon_model_shapefile_path (string): a path to an
            OGR shapefile that has the parameters for the global carbon edge
            model. Each georeferenced feature should have fields 'theta1',
            'theta2', 'theta3', and 'method'

    Returns:
        a tuple of:
            scipy.spatial.cKDTree (georeferenced locally projected model
                points),
            theta_model_parameters (parallel Nx3 array of theta parameters),
            method_model_parameter (parallel N array of model numbers (1..3))
    """

    # Reproject the global model into local coordinate system
    carbon_model_reproject_path = os.path.join(
        local_model_dir, 'local_carbon_shape.shp')
    lulc_projection_wkt = pygeoprocessing.get_raster_info(
        base_raster_path)['projection']
    pygeoprocessing.reproject_vector(
        tropical_forest_edge_carbon_model_shapefile_path, lulc_projection_wkt,
        carbon_model_reproject_path)

    model_shape_ds = gdal.OpenEx(carbon_model_reproject_path)
    model_shape_layer = model_shape_ds.GetLayer()

    kd_points = []
    theta_model_parameters = []
    method_model_parameter = []

    # put all the polygons in the kd_tree because it's fast and simple
    for poly_feature in model_shape_layer:
        poly_geom = poly_feature.GetGeometryRef()
        poly_centroid = poly_geom.Centroid()
        # put in row/col order since rasters are row/col indexed
        kd_points.append([poly_centroid.GetY(), poly_centroid.GetX()])

        theta_model_parameters.append([
            poly_feature.GetField(feature_id) for feature_id in
            ['theta1', 'theta2', 'theta3']])
        method_model_parameter.append(poly_feature.GetField('method'))

    method_model_parameter = numpy.array(
        method_model_parameter, dtype=numpy.int32)
    theta_model_parameters = numpy.array(
        theta_model_parameters, dtype=numpy.float32)

    LOGGER.info('building kd_tree')
    kd_tree = scipy.spatial.cKDTree(kd_points)
    LOGGER.info('done building kd_tree with %d points', len(kd_points))
    return kd_tree, theta_model_parameters, method_model_parameter


def _calculate_tropical_forest_edge_carbon_map(
        edge_distance_path, kd_tree, theta_model_parameters,
        method_model_parameter, n_nearest_model_points,
        biomass_to_carbon_conversion_factor,
        tropical_forest_edge_carbon_map_path):
    """Calculates the carbon on the forest pixels accounting for their global
    position with respect to precalculated edge carbon models.

    Parameters:
        edge_distance_path (string): path to the a raster where each pixel
            contains the pixel distance to forest edge.
        kd_tree (scipy.spatial.cKDTree): a kd-tree that has indexed the valid
            model parameter points for fast nearest neighbor calculations.
        theta_model_parameters (numpy.array Nx3): parallel array of model
            theta parameters consistent with the order in which points were
            inserted into 'kd_tree'
        method_model_parameter (numpy.array N): parallel array of method
            numbers (1..3) consistent with the order in which points were
            inserted into 'kd_tree'.
        n_nearest_model_points (int): number of nearest model points to search
            for.
        biomass_to_carbon_conversion_factor (float): number by which to
            multiply the biomass by to get carbon.
        tropical_forest_edge_carbon_map_path (string): a filepath to the output
            raster which will contain total carbon stocks per cell of forest
            type.

    Returns:
        None"""

    # create output raster and open band for writing
    # fill nodata, in case we skip entire memory blocks that are non-forest
    pygeoprocessing.new_raster_from_base(
        edge_distance_path, tropical_forest_edge_carbon_map_path,
        gdal.GDT_Float32, [CARBON_MAP_NODATA],
        fill_value_list=[CARBON_MAP_NODATA])
    edge_carbon_dataset = gdal.OpenEx(
        tropical_forest_edge_carbon_map_path, gdal.GA_Update)
    edge_carbon_band = edge_carbon_dataset.GetRasterBand(1)
    edge_carbon_geotransform = edge_carbon_dataset.GetGeoTransform()

    # create edge distance band for memory block reading
    n_rows = edge_carbon_dataset.RasterYSize
    n_cols = edge_carbon_dataset.RasterXSize
    n_cells = n_rows * n_cols
    n_cells_processed = 0
    # timer to give updates per call
    last_time = time.time()

    cell_area_ha = pygeoprocessing.get_raster_info(
        edge_distance_path)['mean_pixel_size'] ** 2 / 10000.0
    cell_size_km = pygeoprocessing.get_raster_info(
        edge_distance_path)['mean_pixel_size'] / 1000.0

    # Loop memory block by memory block, calculating the forest edge carbon
    # for every forest pixel.
    for edge_distance_data, edge_distance_block in pygeoprocessing.iterblocks(
            edge_distance_path, largest_block=2**12):
        current_time = time.time()
        if current_time - last_time > 5.0:
            LOGGER.info(
                'carbon edge calculation approx. %.2f%% complete',
                (n_cells_processed / float(n_cells) * 100.0))
            last_time = current_time
        n_cells_processed += (
            edge_distance_data['win_xsize'] * edge_distance_data['win_ysize'])
        valid_edge_distance_mask = (edge_distance_block > 0)

        # if no valid forest pixels to calculate, skip to the next block
        if not valid_edge_distance_mask.any():
            continue

        # calculate local coordinates for each pixel so we can test for
        # distance to the nearest carbon model points
        col_range = numpy.linspace(
            edge_carbon_geotransform[0] +
            edge_carbon_geotransform[1] * edge_distance_data['xoff'],
            edge_carbon_geotransform[0] +
            edge_carbon_geotransform[1] * (
                edge_distance_data['xoff'] + edge_distance_data['win_xsize']),
            num=edge_distance_data['win_xsize'], endpoint=False)
        row_range = numpy.linspace(
            edge_carbon_geotransform[3] +
            edge_carbon_geotransform[5] * edge_distance_data['yoff'],
            edge_carbon_geotransform[3] +
            edge_carbon_geotransform[5] * (
                edge_distance_data['yoff'] + edge_distance_data['win_ysize']),
            num=edge_distance_data['win_ysize'], endpoint=False)
        col_coords, row_coords = numpy.meshgrid(col_range, row_range)

        # query nearest points for every point in the grid
        # n_jobs=-1 means use all available cpus
        coord_points = zip(
            row_coords[valid_edge_distance_mask].ravel(),
            col_coords[valid_edge_distance_mask].ravel())
        # note, the 'n_jobs' parameter was introduced in SciPy 0.16.0
        distances, indexes = kd_tree.query(
            coord_points, k=n_nearest_model_points,
            distance_upper_bound=DISTANCE_UPPER_BOUND, n_jobs=-1)

        if n_nearest_model_points == 1:
            distances = distances.reshape(distances.shape[0], 1)
            indexes = indexes.reshape(indexes.shape[0], 1)

        # the 3 is for the 3 thetas in the carbon model
        thetas = numpy.zeros((indexes.shape[0], indexes.shape[1], 3))
        valid_index_mask = (indexes != kd_tree.n)
        thetas[valid_index_mask] = theta_model_parameters[
            indexes[valid_index_mask]]

        # the 3 is for the 3 models (asym, exp, linear)
        biomass_model = numpy.zeros(
            (indexes.shape[0], indexes.shape[1], 3))
        # reshape to an N,nearest_points so we can multiply by thetas
        valid_edge_distances_km = numpy.repeat(
            edge_distance_block[valid_edge_distance_mask] * cell_size_km,
            n_nearest_model_points).reshape(-1, n_nearest_model_points)

        # asymptotic model
        # biomass_1 = t1 - t2 * exp(-t3 * edge_dist_km)
        biomass_model[:, :, 0] = (
            thetas[:, :, 0] - thetas[:, :, 1] * numpy.exp(
                -thetas[:, :, 2] * valid_edge_distances_km)
            ) * cell_area_ha

        # logarithmic model
        # biomass_2 = t1 + t2 * numpy.log(edge_dist_km)
        biomass_model[:, :, 1] = (
            thetas[:, :, 0] + thetas[:, :, 1] * numpy.log(
                valid_edge_distances_km)) * cell_area_ha

        # linear regression
        # biomass_3 = t1 + t2 * edge_dist_km
        biomass_model[:, :, 2] = (
            (thetas[:, :, 0] + thetas[:, :, 1] * valid_edge_distances_km) *
            cell_area_ha)

        # Collapse the biomass down to the valid models
        model_index = numpy.zeros(indexes.shape, dtype=numpy.int8)
        model_index[valid_index_mask] = (
            method_model_parameter[indexes[valid_index_mask]] - 1)

        # reduce the axis=1 dimensionality of the model by selecting the
        # appropriate value via the model_index array
        # Got this trick from http://stackoverflow.com/questions/18702746/reduce-a-dimension-of-numpy-array-by-selecting
        biomass_y, biomass_x = numpy.meshgrid(
            numpy.arange(biomass_model.shape[1]),
            numpy.arange(biomass_model.shape[0]))
        biomass = biomass_model[biomass_x, biomass_y, model_index]

        # reshape the array so that each set of points is in a separate
        # dimension, here distances are distances to each valid model
        # point, not distance to edge of forest
        weights = numpy.zeros(distances.shape)
        valid_distance_mask = (distances > 0) & (distances < numpy.inf)
        weights[valid_distance_mask] = (
            n_nearest_model_points / distances[valid_distance_mask])

        # Denominator is the sum of the weights per nearest point (axis 1)
        denom = numpy.sum(weights, axis=1)
        # To avoid a divide by 0
        valid_denom = denom != 0
        average_biomass = numpy.zeros(distances.shape[0])
        average_biomass[valid_denom] = (
            numpy.sum(weights[valid_denom] *
                      biomass[valid_denom], axis=1) / denom[valid_denom])

        # Ensure the result has nodata everywhere the distance was invalid
        result = numpy.empty(
            edge_distance_block.shape, dtype=numpy.float32)
        result[:] = CARBON_MAP_NODATA
        # convert biomass to carbon in this stage
        result[valid_edge_distance_mask] = (
            average_biomass * biomass_to_carbon_conversion_factor)
        edge_carbon_band.WriteArray(
            result, xoff=edge_distance_data['xoff'],
            yoff=edge_distance_data['yoff'])
    LOGGER.info('carbon edge calculation 100.0% complete')


@validation.invest_validator
def validate(args, limit_to=None):
    """Validate args to ensure they conform to `execute`'s contract.

    Parameters:
        args (dict): dictionary of key(str)/value pairs where keys and
            values are specified in `execute` docstring.
        limit_to (str): (optional) if not None indicates that validation
            should only occur on the args[limit_to] value. The intent that
            individual key validation could be significantly less expensive
            than validating the entire `args` dictionary.

    Returns:
        list of ([invalid key_a, invalid_keyb, ...], 'warning/error message')
            tuples. Where an entry indicates that the invalid keys caused
            the error message in the second part of the tuple. This should
            be an empty list if validation succeeds.
    """
    missing_key_list = []
    no_value_list = []
    validation_error_list = []

    for key in [
            'workspace_dir',
            'lulc_raster_path',
            'biophysical_table_uri',
            'pools_to_calculate',
            'compute_forest_edge_effects',
            'tropical_forest_edge_carbon_model_shape_uri',
            'n_nearest_model_points',
            'biomass_to_carbon_conversion_factor',
            ]:
        if limit_to is None or limit_to == key:
            if key not in args:
                missing_key_list.append(key)
            elif args[key] in ['', None]:
                no_value_list.append(key)

    if len(missing_key_list) > 0:
        # if there are missing keys, we have raise KeyError to stop hard
        raise KeyError(
            "The following keys were expected in `args` but were missing" +
            ', '.join(missing_key_list))

    if len(no_value_list) > 0:
        validation_error_list.append(
            (no_value_list, 'parameter has no value'))

    # check required files exist
    for key in [
            'lulc_raster_path',
            'biophysical_table_uri']:
        if (limit_to is None or limit_to == key) and (
                not os.path.exists(args[key])):
            validation_error_list.append(
                ([key], 'not found on disk'))

    optional_file_type_list = [('lulc_raster_path', 'raster', True)]
    if args['compute_forest_edge_effects']:
        optional_file_type_list.extend(
            [('tropical_forest_edge_carbon_model_shape_uri', 'vector', True),
             ('aoi_vector_path', 'vector', False)])

    # check that existing/optional files are the correct types
    with utils.capture_gdal_logging():
        for key, key_type, required in optional_file_type_list:
            if (limit_to is None or limit_to == key) and key in args:
                if len(args[key]) == 0 and not required:
                    continue

                if not os.path.exists(args[key]):
                    validation_error_list.append(
                        ([key], 'not found on disk'))
                    continue
                if key_type == 'raster':
                    raster = gdal.OpenEx(args[key])
                    if raster is None:
                        validation_error_list.append(
                            ([key], 'not a raster'))
                    del raster
                elif key_type == 'vector':
                    vector = gdal.OpenEx(args[key])
                    if vector is None:
                        validation_error_list.append(
                            ([key], 'not a vector'))
                    del vector

    return validation_error_list
