"""InVEST Seasonal Water Yield Model."""

import itertools
import heapq
import struct
import collections
import os
import logging
import re
import uuid
import warnings

import numpy
from osgeo import gdal
from osgeo import ogr
import pygeoprocessing
import pygeoprocessing.routing
import pygeoprocessing.routing.routing_core
import utils

logging.basicConfig(format='%(asctime)s %(name)-20s %(levelname)-8s \
%(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

LOGGER = logging.getLogger(
    'natcap.invest.seasonal_water_yield.seasonal_water_yield')

_N_MONTHS = 12
MONTH_ID_TO_LABEL = [
    'jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct',
    'nov', 'dec']

_OUTPUT_BASE_FILES = {
    'l_path_list': ['l_%d.tif' % _ for _ in xrange(_N_MONTHS)],
    'l_annual_path': 'l_annual.tif',
    }

_INTERMEDIATE_BASE_FILES = {
    'aetm_path_list': ['aetm_%d.tif' % (_+1) for _ in xrange(_N_MONTHS)],
    'flow_dir_path': 'flow_dir.tif',
    'ti_path': 'ti.tif',
    'outflow_direction_path': 'outflow_direction.tif',
    'outflow_weights_path': 'outflow_weights.tif',
    'lulc_aligned_path': 'lulc_aligned.tif',
    'dem_aligned_path': 'dem_aligned.tif',
    'lulc_valid_path': 'lulc_valid.tif',
    'dem_valid_path': 'dem_valid.tif',
    'slope_path': 'slope.tif',
    'soil_depth_aligned_path': 'soil_depth.tif',
    'flow_accum_path': 'flow_accum.tif',
    'root_depth_path': 'root_depth.tif',
    'precip_path_aligned_list': [
        'prcp_a%d.tif' % _ for _ in xrange(_N_MONTHS)],
    'valid_precip_path_list':  [
        'valid_prcp_a%d.tif' % _ for _ in xrange(_N_MONTHS)],
    'n_events_path_list': ['n_events_%d.tif' % _ for _ in xrange(_N_MONTHS)],
    'et0_path_aligned_list': ['et0_a_%d.tif' % _ for _ in xrange(_N_MONTHS)],
    'pet_path_aligned_list': ['pet_%d.tif' % _ for _ in xrange(_N_MONTHS)],
    'pet_annual_path': 'pet_annual.tif',
    'kc_path_list': ['kc_%d.tif' % _ for _ in xrange(_N_MONTHS)],
    'z_rm_path_list': ['z_rm_%d.tif' % _ for _ in xrange(_N_MONTHS)],
    'pawc_aligned_path': 'pawc_aligned.tif',
    'wm_path_list': ['wm_%d.tif' % _ for _ in xrange(_N_MONTHS)],
    'l1_path_list': ['l1_%d.tif' % _ for _ in xrange(_N_MONTHS)],
    'l1_annual_path': 'l1_annual.tif',
    'l2_path_list': ['l2_%d.tif' % _ for _ in xrange(_N_MONTHS)],
    'l2_annual_path': 'l2_annual.tif',
    'cz_aligned_raster_path': 'cz_aligned.tif',
    'subsidized_path_list': [
        'subsidized_%d.tif' % _ for _ in xrange(_N_MONTHS)],
    'subsidized_annual_path': 'subsidized_annual.tif',
    'temporary_subwatershed_path': 'temporary_subwatershed.shp',
    'l1_upstream_path_list': [
        'l1_upstream_sum_%d.tif' % _ for _ in xrange(_N_MONTHS)],
    'l1_annual_upstream_path': 'l1_annual_upstream.tif',
}

_TMP_BASE_FILES = {
    }

ROOT_DEPTH_NODATA = -1.0
N_EVENTS_NODATA = -1
QF_NODATA = -1.0
KC_NODATA = -1.0
SI_NODATA = -1.0
CN_NODATA = -1.0
AET_NODATA = -1.0
PET_NODATA = -1.0
L1_NODATA = -1.0
L2_NODATA = -1.0
L_NODATA = -1.0
TI_NODATA = -1.0
SUBSIDIZED_NODATA = -1.0


def execute(args):
    """InVEST seasonal water yield model.

    This function invokes the InVEST seasonal water yield model described in
    "Spatial attribution of baseflow generation at the parcel level for
    ecosystem-service valuation", Guswa, et. al (under review in "Water
    Resources Research")

    Parameters:
        args['workspace_dir'] (string): output directory for intermediate,
        temporary, and final files
        args['results_suffix'] (string): (optional) string to append to any
            output files
        args['et0_dir'] (string): Path to a directory that contains rasters of
            monthly reference evapotranspiration; units in mm.
        args['precip_dir'] (string): Path to a directory that contains rasters
            of monthly precipitation; units in mm.
        args['dem_raster_path'] (string): a path to a digital elevation
            raster.
        args['pawc_raster_path'] (string): a path to the plant available water
            content raster (PAWC) where each pixel indicates the plant
            available water content as a fraction between 0 and 1.
        args['lulc_raster_path'] (string): a path to a land cover raster used
            to classify biophysical properties of pixels.
        args['soil_depth_raster_path'] (string): path to raster that indicates
            soil depth per pixel in mm.
        args['aoi_path'] (string): path to a vector that indicates the area
            over which the model should be run, as well as the area in which
            to aggregate over when calculating the output Qb.
        args['biophysical_table_path'] (string): path to a CSV table that maps
            landcover codes to monthly Kc values. Headers must include
            'lucode', 'Kc_1', 'Kc_2', 'Kc_3', 'Kc_4', 'Kc_5', 'Kc_6', 'Kc_7',
            'Kc_8', 'Kc_9', 'Kc_10', 'Kc_11', 'Kc_12'.
        args['qb_value'] (string/float): value for Qb for the watershed.
    """
    # This upgrades warnings to exceptions across this model.
    # I found this useful to catch all kinds of weird inputs to the model
    # during debugging and think it makes sense to have in production of this
    # model too.
    try:
        warnings.filterwarnings('error')
        _execute(args)
    finally:
        warnings.resetwarnings()


def _execute(args):
    """Execute the seasonal water yield model.

    Parameters:
        See the parameters for
        `natcap.invest.seasonal_water_yield.seasonal_wateryield.execute`.

    Returns:
        None
    """
    LOGGER.info('prepare and test inputs for common errors')

    # fail early if the AOI is not a single polygon
    aoi_vector = ogr.Open(args['aoi_path'])
    aoi_n_layers = aoi_vector.GetLayerCount()
    if aoi_n_layers != 1:
        raise ValueError(
            "Expected one AOI layer but got %d" % aoi_n_layers)
    aoi_layer = aoi_vector.GetLayer()
    aoi_geom_type = aoi_layer.GetGeomType()
    if aoi_geom_type not in [ogr.wkbMultiPolygon, ogr.wkbPolygon]:
        raise ValueError(
            "Expected an AOI polygon type, but got %d" % aoi_geom_type)
    n_features = aoi_layer.GetFeatureCount()
    if n_features != 1:
        raise ValueError(
            "Expected 1 polygon in AOI layer but got %d" % n_features)
    aoi_layer = None
    aoi_vector = None

    # fail early if the climate zone table has more than 1 row but no raster
    cz_rain_events_lookup = (
        pygeoprocessing.get_lookup_from_table(
            args['climate_zone_table_path'], 'cz_id'))
    if ('climate_zone_raster_path' not in args and
            len(cz_rain_events_lookup) != 1):
        raise ValueError(
            "No climate zone raster was defined, but %d climate zones "
            "were defined in the table." % len(cz_rain_events_lookup))

    biophysical_table = pygeoprocessing.get_lookup_from_table(
        args['biophysical_table_path'], 'lucode')

    qb_value = float(args['qb_value'])
    pixel_size = pygeoprocessing.get_cell_size_from_uri(
        args['dem_raster_path'])
    file_suffix = utils.make_suffix_string(args, 'results_suffix')
    intermediate_output_dir = os.path.join(
        args['workspace_dir'], 'intermediate_outputs')
    output_dir = args['workspace_dir']
    pygeoprocessing.create_directories(
        [intermediate_output_dir, output_dir])

    LOGGER.info('Building file registry')
    file_registry = utils.build_file_registry(
        [(_OUTPUT_BASE_FILES, output_dir),
         (_INTERMEDIATE_BASE_FILES, intermediate_output_dir),
         (_TMP_BASE_FILES, output_dir)], file_suffix)

    LOGGER.info('Aligning and clipping dataset list')
    input_align_list = [
        args['lulc_raster_path'],
        args['dem_raster_path'],
        args['pawc_raster_path'],
        args['soil_depth_raster_path']]
    output_align_list = [
        file_registry['lulc_aligned_path'],
        file_registry['dem_aligned_path'],
        file_registry['pawc_aligned_path'],
        file_registry['soil_depth_aligned_path']]

    # make sure climate zone raster path isn't an empty string
    if ('climate_zone_raster_path' in args and
            args['climate_zone_raster_path']):
        input_align_list.append(args['climate_zone_raster_path'])
        output_align_list.append(file_registry['cz_aligned_raster_path'])

    precip_path_list = []
    et0_path_list = []

    et0_dir_list = [
        os.path.join(args['et0_dir'], f) for f in os.listdir(args['et0_dir'])]
    precip_dir_list = [
        os.path.join(args['precip_dir'], f) for f in os.listdir(
            args['precip_dir'])]

    for month_index in range(1, _N_MONTHS + 1):
        month_file_match = re.compile(r'.*[^\d]%d\.[^.]+$' % month_index)
        for data_type, dir_list, path_list in [
                ('et0', et0_dir_list, et0_path_list),
                ('Precip', precip_dir_list, precip_path_list)]:
            file_list = [
                month_file_path for month_file_path in dir_list
                if month_file_match.match(month_file_path)]
            if len(file_list) == 0:
                raise ValueError(
                    "No %s found for month %d" % (data_type, month_index))
            if len(file_list) > 1:
                raise ValueError(
                    "Ambiguous set of files found for month %d: %s" %
                    (month_index, file_list))
            path_list.append(file_list[0])

    input_align_list = precip_path_list + et0_path_list + input_align_list
    output_align_list = (
        file_registry['precip_path_aligned_list'] +
        file_registry['et0_path_aligned_list'] + output_align_list)

    align_index = len(input_align_list) - 1  # this aligns with the DEM
    interpolate_list = ['nearest'] * len(input_align_list)
    LOGGER.debug(input_align_list)
    LOGGER.debug(output_align_list)
    pygeoprocessing.align_dataset_list(
        input_align_list, output_align_list, interpolate_list, pixel_size,
        'intersection', align_index, aoi_uri=args['aoi_path'],
        assert_datasets_projected=True)

    if ('climate_zone_raster_path' not in args or
            not args['climate_zone_raster_path']):
        # Create a constant climate zone raster, there is only one entry so
        # use that CZ id for the fill value in the raster
        climate_zone_id = cz_rain_events_lookup.keys()[0]
        # the following makes sure we don't integer overflow
        if climate_zone_id > 0:
            climate_zone_nodata = climate_zone_id - 1
        else:
            climate_zone_nodata = climate_zone_id + 1
        pygeoprocessing.new_raster_from_base_uri(
            file_registry['dem_aligned_path'],
            file_registry['cz_aligned_raster_path'], 'GTiff',
            climate_zone_nodata, gdal.GDT_Int32, fill_value=climate_zone_id)

    # sometimes users input data where the DEM is defined in places where the
    # land cover isn't, mask those out
    LOGGER.info("Masking invalid lulc, dem, and possible soil group overlap")
    input_raster_path_list = [
        file_registry['dem_aligned_path'],
        file_registry['lulc_aligned_path']] + (
            file_registry['precip_path_aligned_list'])
    output_valid_raster_path_list = [
        file_registry['dem_valid_path'],
        file_registry['lulc_valid_path']] + (
            file_registry['valid_precip_path_list'])
    _mask_any_nodata(input_raster_path_list, output_valid_raster_path_list)

    LOGGER.info('Mapping LULC to root depth')
    root_depth_lookup = dict([
        (lucode, biophysical_table[lucode]['root_depth'])
        for lucode in biophysical_table])
    pygeoprocessing.reclassify_dataset_uri(
        file_registry['lulc_valid_path'], root_depth_lookup,
        file_registry['root_depth_path'], gdal.GDT_Float32, ROOT_DEPTH_NODATA,
        exception_flag='values_required', assert_dataset_projected=True)

    kc_lookup = {}
    for month_index in xrange(_N_MONTHS):
        LOGGER.info("For month %d: ", month_index)
        LOGGER.info('classify kc')
        kc_lookup = dict([
            (lucode, biophysical_table[lucode]['kc_%d' % (month_index+1)])
            for lucode in biophysical_table])
        pygeoprocessing.reclassify_dataset_uri(
            file_registry['lulc_valid_path'], kc_lookup,
            file_registry['kc_path_list'][month_index], gdal.GDT_Float32,
            KC_NODATA)
        month_label = MONTH_ID_TO_LABEL[month_index]
        climate_zone_rain_events_month = dict([
            (cz_id, cz_rain_events_lookup[cz_id][month_label]) for
            cz_id in cz_rain_events_lookup])
        pygeoprocessing.reclassify_dataset_uri(
            file_registry['cz_aligned_raster_path'],
            climate_zone_rain_events_month,
            file_registry['n_events_path_list'][month_index],
            gdal.GDT_Float32, N_EVENTS_NODATA)

        _calculate_aet_uphill(
            file_registry['valid_precip_path_list'][month_index],
            file_registry['kc_path_list'][month_index],
            file_registry['et0_path_aligned_list'][month_index],
            file_registry['pawc_aligned_path'],
            file_registry['n_events_path_list'][month_index],
            file_registry['root_depth_path'],
            file_registry['pet_path_aligned_list'][month_index],
            file_registry['wm_path_list'][month_index],
            file_registry['z_rm_path_list'][month_index],
            file_registry['aetm_path_list'][month_index])

        LOGGER.info("Calculate L1")
        _calculate_l1(
            file_registry['valid_precip_path_list'][month_index],
            file_registry['aetm_path_list'][month_index],
            file_registry['l1_path_list'][month_index]
            )

    LOGGER.info('flow direction')
    pygeoprocessing.routing.flow_direction_d_inf(
        file_registry['dem_valid_path'], file_registry['flow_dir_path'])

    LOGGER.info('flow weights')
    pygeoprocessing.routing.routing_core.calculate_flow_weights(
        file_registry['flow_dir_path'],
        file_registry['outflow_weights_path'],
        file_registry['outflow_direction_path'])

    LOGGER.info('flow accumulation')
    pygeoprocessing.routing.flow_accumulation(
        file_registry['flow_dir_path'], file_registry['dem_valid_path'],
        file_registry['flow_accum_path'])

    LOGGER.info('calculate slope')
    pygeoprocessing.calculate_slope(
        file_registry['dem_valid_path'], file_registry['slope_path'])

    LOGGER.info("Calculate TI")
    _calculate_ti(
        file_registry['flow_accum_path'], file_registry['slope_path'],
        file_registry['soil_depth_aligned_path'], file_registry['ti_path'])

    LOGGER.info("calculate subsidized area")

    def _l2_op(l1_array, pet_array):
        """Create L2 as -PET only where L1 is valid."""
        result = numpy.empty(l1_array.shape)
        result[:] = L2_NODATA
        valid_mask = (l1_array != L1_NODATA) & (pet_array != PET_NODATA)
        result[valid_mask] = -pet_array[valid_mask]
        return result

    for month_index in xrange(_N_MONTHS):
        LOGGER.info("For month %d: ", month_index)
        LOGGER.info("calculate L2_%d", month_index)
        _calculate_upstream_flow(
            file_registry['flow_dir_path'],
            file_registry['dem_valid_path'],
            file_registry['l1_path_list'][month_index],
            args['aoi_path'],
            file_registry['l1_upstream_path_list'][month_index])
        pygeoprocessing.vectorize_datasets(
            [file_registry['l1_path_list'][month_index],
             file_registry['pet_path_aligned_list'][month_index]], _l2_op,
            file_registry['l2_path_list'][month_index], gdal.GDT_Float32,
            L2_NODATA, pixel_size, "intersection",
            datasets_are_pre_aligned=True, vectorize_op=False)
        LOGGER.info("calculate subsided area for month %d", month_index)
        _calculate_subsidized_area(
            qb_value, file_registry['l1_path_list'][month_index],
            file_registry['l1_upstream_path_list'][month_index],
            file_registry['l2_path_list'][month_index],
            file_registry['pet_path_aligned_list'][month_index],
            file_registry['ti_path'], args['aoi_path'],
            file_registry['subsidized_path_list'][month_index])
        LOGGER.info("calculate L_%d", month_index)
        _calculate_l(
            file_registry['l1_path_list'][month_index],
            file_registry['l2_path_list'][month_index],
            file_registry['subsidized_path_list'][month_index],
            file_registry['l_path_list'][month_index])

    LOGGER.info("Annual calculation")

    def _sum_op(nodata_value):
        """Make a sum operation that operates on user defined nodata value."""
        def __sum_op(*arrays):
            """Sum values in arrays assuming they are nodata aligned."""
            result = numpy.empty(arrays[0].shape, dtype=numpy.float32)
            result[:] = nodata_value
            valid_mask = arrays[0] != nodata_value
            result[valid_mask] = numpy.sum(
                [_[valid_mask] for _ in arrays], axis=0)
            return result
        return __sum_op

    # Calc annual L1
    pygeoprocessing.vectorize_datasets(
        file_registry['l1_path_list'], _sum_op(L1_NODATA),
        file_registry['l1_annual_path'], gdal.GDT_Float32, L1_NODATA,
        pixel_size, "intersection", datasets_are_pre_aligned=True,
        vectorize_op=False)
    # calc annual L2
    pygeoprocessing.vectorize_datasets(
        file_registry['l2_path_list'], _sum_op(L2_NODATA),
        file_registry['l2_annual_path'], gdal.GDT_Float32, L2_NODATA,
        pixel_size, "intersection", datasets_are_pre_aligned=True,
        vectorize_op=False)
    # calc annual PET
    pygeoprocessing.vectorize_datasets(
        file_registry['pet_path_aligned_list'], _sum_op(PET_NODATA),
        file_registry['pet_annual_path'], gdal.GDT_Float32, PET_NODATA,
        pixel_size, "intersection", datasets_are_pre_aligned=True,
        vectorize_op=False)

    _calculate_upstream_flow(
        file_registry['flow_dir_path'], file_registry['dem_valid_path'],
        file_registry['l1_annual_path'], args['aoi_path'],
        file_registry['l1_annual_upstream_path'])
    _calculate_subsidized_area(
        qb_value, file_registry['l1_annual_path'],
        file_registry['l1_annual_upstream_path'],
        file_registry['l2_annual_path'], file_registry['pet_annual_path'],
        file_registry['ti_path'], args['aoi_path'],
        file_registry['subsidized_annual_path'])
    _calculate_l(
        file_registry['l1_annual_path'], file_registry['l2_annual_path'],
        file_registry['subsidized_annual_path'],
        file_registry['l_annual_path'])

    LOGGER.info('deleting temporary files')
    for file_id in _TMP_BASE_FILES:
        try:
            if isinstance(file_registry[file_id], basestring):
                os.remove(file_registry[file_id])
            elif isinstance(file_registry[file_id], list):
                for index in xrange(len(file_registry[file_id])):
                    os.remove(file_registry[file_id][index])
        except OSError:
            # Let it go.
            pass

    LOGGER.info('  (\\w/)  SWY Complete!')
    LOGGER.info('  (..  \\ ')
    LOGGER.info(' _/  )  \\______')
    LOGGER.info('(oo /\'\\        )`,')
    LOGGER.info(' `--\' (v  __( / ||')
    LOGGER.info('       |||  ||| ||')
    LOGGER.info('      //_| //_|')


def _mask_any_nodata(input_raster_path_list, output_raster_path_list):
    """Mask local pixel stacks that include nodata anywhere in the stack.

    Parameters:
        input_raster_path_list (list): list of input raster paths, all rasters
            are of the same projection, shape, and cell pixel_size
        output_raster_path_list (list): a parallel list to
            `input_raster_path_list` to hold the masked results of each input
            file

    Returns:
        None
    """
    base_nodata_list = [pygeoprocessing.get_nodata_from_uri(
        path) for path in input_raster_path_list]
    pixel_size = pygeoprocessing.get_cell_size_from_uri(
        input_raster_path_list[0])
    nodata_list = None
    for index in xrange(len(input_raster_path_list)):
        nodata_list = base_nodata_list[index:] + base_nodata_list[:index]

        def mask_if_not_both_valid(*value_list):
            """If values are nodata, nodata_list[0], else `value_list[0]`."""
            valid_mask = numpy.empty(value_list[0].shape, dtype=numpy.bool)
            valid_mask[:] = True
            for value_index in xrange(len(value_list)):
                valid_mask &= (
                    value_list[value_index] != nodata_list[value_index])
            return numpy.where(valid_mask, value_list[0], nodata_list[0])

        pygeoprocessing.vectorize_datasets(
            input_raster_path_list[index:]+input_raster_path_list[:index],
            mask_if_not_both_valid, output_raster_path_list[index],
            gdal.GDT_Float32, nodata_list[0], pixel_size, 'intersection',
            vectorize_op=False, datasets_are_pre_aligned=True)


def _calculate_aet_uphill(
        precip_path, kc_path, et0_path, pawc_path, n_events_path,
        root_depth_path, out_pet_path, out_wm_path, out_z_rm_path,
        out_aet_path):
    """Calculate uphill AET (case 1 in notes).

    Parameters:
        precip_path: (string) path to precipitation raster (mm).
        kc_path: (string) path to the per pixel Kc values.
        et0_path: (string) path to monthly potential evapotranspiration
            raster (mm).
        pawc_path: (string) path to the plant available water content raster
            (unitless).
        n_events_path: (string) path to raster that shows number of rain
            events per month on that pixel (unitless).
        root_depth_path: (string) path to root depth raster.
        out_pet_path: (string) path to output PET.
        out_wm_path: (string) path to output P/PET ratio.
        out_aet_path: (string) path to output actual evapotranspiration raster
            (mm).

    Returns:
        None.
    """
    LOGGER.info("calculate PET")
    et0_nodata = pygeoprocessing.get_nodata_from_uri(et0_path)
    kc_nodata = pygeoprocessing.get_nodata_from_uri(kc_path)
    pixel_size = pygeoprocessing.get_cell_size_from_uri(precip_path)

    def _pet_op(et0_array, kc_array):
        """Calculate PET."""
        result = numpy.empty(et0_array.shape, dtype=numpy.float32)
        result[:] = PET_NODATA
        valid_mask = ((et0_array != et0_nodata) & (kc_array != kc_nodata))
        result[valid_mask] = et0_array[valid_mask] * kc_array[valid_mask]
        return result

    pygeoprocessing.vectorize_datasets(
        [et0_path, kc_path], _pet_op, out_pet_path,
        gdal.GDT_Float32, PET_NODATA, pixel_size, 'intersection',
        vectorize_op=False, datasets_are_pre_aligned=True)

    LOGGER.info("calculate Wm")
    precip_nodata = pygeoprocessing.get_nodata_from_uri(precip_path)

    def _w_m_op(precip_array, pet_array):
        """Calculate W_m in equation 2a."""
        result = numpy.empty(precip_array.shape, dtype=numpy.float32)
        result[:] = AET_NODATA
        valid_mask = (
            (precip_array != precip_nodata) &
            (pet_array != PET_NODATA) &
            (pet_array != 0.0))

        result[valid_mask] = precip_array[valid_mask] / pet_array[valid_mask]
        result[valid_mask][pet_array[valid_mask] == 0] = 0.0
        return result

    pygeoprocessing.vectorize_datasets(
        [precip_path, out_pet_path], _w_m_op, out_wm_path,
        gdal.GDT_Float32, AET_NODATA, pixel_size, 'intersection',
        vectorize_op=False, datasets_are_pre_aligned=True)

    LOGGER.info("calculate Z*_r,m")
    root_depth_nodata = pygeoprocessing.get_nodata_from_uri(root_depth_path)
    n_events_nodata = pygeoprocessing.get_nodata_from_uri(n_events_path)
    pawc_nodata = pygeoprocessing.get_nodata_from_uri(pawc_path)

    def _z_rm_op(root_depth_array, pawc_array, n_events_array, precip_array):
        """Calculate Z^*_r,m from Eq. (2a)."""
        result = numpy.empty(precip_array.shape, dtype=numpy.float32)
        result[:] = AET_NODATA  # aet nodata is reasonable for z_rm
        valid_mask = (
            (root_depth_array != root_depth_nodata) &
            (pawc_array != pawc_nodata) &
            (n_events_array != n_events_nodata) &
            (precip_array != precip_nodata) &
            (precip_array != 0.0))
        result[valid_mask] = (
            root_depth_array[valid_mask] *
            pawc_array[valid_mask] *
            n_events_array[valid_mask] /
            precip_array[valid_mask])
        precip_zero_mask = precip_array == 0
        result[precip_zero_mask] = 0
        return result

    pygeoprocessing.vectorize_datasets(
        [root_depth_path, pawc_path, n_events_path, precip_path], _z_rm_op,
        out_z_rm_path, gdal.GDT_Float32, AET_NODATA, pixel_size,
        'intersection', vectorize_op=False, datasets_are_pre_aligned=True)

    LOGGER.info("calculating AET1")

    def _aet1_op(pet_array, root_depth_array, w_m_array, z_rm_array):
        """Calculate AET1."""
        result = numpy.empty(pet_array.shape, dtype=numpy.float32)
        result[:] = AET_NODATA
        valid_mask = (
            (pet_array != PET_NODATA) &
            (root_depth_array != root_depth_nodata) &
            (w_m_array != AET_NODATA) & (z_rm_array != AET_NODATA))

        try:
            # This section isolates positive versus negative exponents so that
            # numpy.exp doesn't overflow with a large positive value
            exponent = z_rm_array[valid_mask] * (1 - w_m_array[valid_mask])
            pos_exponent = exponent >= 0
            numerator = numpy.empty(exponent.shape, dtype=numpy.float32)
            denominator = numpy.empty(exponent.shape, dtype=numpy.float32)
            numerator[pos_exponent] = 1 - numpy.exp(-exponent[pos_exponent])
            numerator[~pos_exponent] = numpy.exp(exponent[~pos_exponent]) - 1

            denominator[pos_exponent] = (
                1 - w_m_array[valid_mask][pos_exponent] *
                numpy.exp(-exponent[pos_exponent]))
            denominator[~pos_exponent] = (
                numpy.exp(exponent[~pos_exponent]) -
                w_m_array[valid_mask][~pos_exponent])

            # this section calculates the fraction on the reduced set so it
            # can be copied back to the full result
            denom_nonzero_mask = denominator != 0
            intermediate_result = numpy.empty(
                exponent.shape, dtype=numpy.float32)
            intermediate_result[denom_nonzero_mask] = (
                pet_array[valid_mask][denom_nonzero_mask] *
                w_m_array[valid_mask][denom_nonzero_mask] * (
                    numerator[denom_nonzero_mask] /
                    denominator[denom_nonzero_mask]))
            intermediate_result[~denom_nonzero_mask] = 0.0
            result[valid_mask] = intermediate_result

            return result
        except RuntimeWarning:
            LOGGER.debug(numpy.sort(z_rm_array[valid_mask]))
            LOGGER.debug(numpy.sort(w_m_array[valid_mask]))
            raise

    pygeoprocessing.vectorize_datasets(
        [out_pet_path, root_depth_path, out_wm_path, out_z_rm_path], _aet1_op,
        out_aet_path, gdal.GDT_Float32, AET_NODATA, pixel_size,
        'intersection', vectorize_op=False, datasets_are_pre_aligned=True)


def _calculate_l1(precip_path, aet_path, l1_out_path):
    """Calculate L1 as P-AET1.

    Parameters:
        precip_path (string): path to precipitation raster.
        aet_path (string): path to AET raster.
        l1_out_path (string): path to output L1 raster.

    Returns:
        None.
    """
    precip_nodata = pygeoprocessing.get_nodata_from_uri(precip_path)
    aet_nodata = pygeoprocessing.get_nodata_from_uri(aet_path)

    def _l1_op(precip_array, aet_array):
        """Calculate L1."""
        valid_mask = (
            (precip_array != precip_nodata) &
            (aet_array != aet_nodata))
        result = numpy.empty(precip_array.shape, dtype=numpy.float32)
        result[:] = L1_NODATA
        result[valid_mask] = (
            precip_array[valid_mask] - aet_array[valid_mask])
        return result

    pixel_size = pygeoprocessing.get_cell_size_from_uri(precip_path)
    pygeoprocessing.vectorize_datasets(
        [precip_path, aet_path], _l1_op, l1_out_path, gdal.GDT_Float32,
        L1_NODATA, pixel_size, 'intersection', vectorize_op=False,
        datasets_are_pre_aligned=True)


def _calculate_ti(
        flow_accum_path, slope_path, soil_depth_path, ti_out_path):
    """Calculate topographic index from Eq. 5.

    Parameters:
        flow_accum_path (string): path to flow accumulation raster.
        slope_path (string): path to slope raster in units (1/1).
        soil_depth_path (string): path to raster of soil depth in mm.

    Returns:
        None.
    """
    flow_nodata = pygeoprocessing.get_nodata_from_uri(flow_accum_path)
    slope_nodata = pygeoprocessing.get_nodata_from_uri(slope_path)
    soil_depth_nodata = pygeoprocessing.get_nodata_from_uri(soil_depth_path)
    cell_size = pygeoprocessing.get_cell_size_from_uri(flow_accum_path)

    def _ti_op(flow_accum_array, slope_array, soil_depth_array):
        """Calculate ti.

        Returns:
            Eq 5: TI = ln( A / (tan(beta) D)
        """
        result = numpy.empty(flow_accum_array.shape, dtype=numpy.float32)
        zero_mask = (
            (slope_array == 0) |
            (soil_depth_array == 0) |
            (flow_accum_array == 0))
        valid_mask = (
            (flow_accum_array != flow_nodata) &
            (slope_array != slope_nodata) &
            (soil_depth_array != soil_depth_nodata) & ~zero_mask)
        result[:] = TI_NODATA
        result[zero_mask] = 0.0
        result[valid_mask] = numpy.log(
            (flow_accum_array[valid_mask] / cell_size) / (
                slope_array[valid_mask] * soil_depth_array[valid_mask]))
        return result

    pygeoprocessing.vectorize_datasets(
        [flow_accum_path, slope_path, soil_depth_path], _ti_op, ti_out_path,
        gdal.GDT_Float32, TI_NODATA, cell_size, 'intersection',
        vectorize_op=False, datasets_are_pre_aligned=True)


def _calculate_subsidized_area(
        qb_val, l1_path, l1_upstream_sum_path, l2_path, pet_path, ti_path,
        watershed_path, subsidized_out_path):
    """Calculated subsidized area such that Eq. 4 is balanced.

    Parameters:
        qb_val (float): Qb value to balance subsidized area to.
        l1_path (string): path to flow raster for upland regions.
        l1_upstream_sum_path (string): path to sum of L1 values upstream.
        l2_path (string): path to subsidized flow raster.
        pet_path (string): path to PET raster.
        ti_path (string): path to topographical index raster.
        watershed_path (string): path to a shapefile that defines the
            watersheds.
        subsidized_out_path (string): path to output raster that masks out
            the subsidized region.

    Returns:
        None.
    """
    (_, n_cols) = pygeoprocessing.get_row_col_from_uri(l1_path)
    array_sorter = _OutOfCoreNumpyArray(
        os.path.dirname(subsidized_out_path), 'ti_array')

    l1_raster = gdal.Open(l1_path)
    l1_band = l1_raster.GetRasterBand(1)
    l2_raster = gdal.Open(l2_path)
    l2_band = l2_raster.GetRasterBand(1)

    l1_upstream_raster = gdal.Open(l1_upstream_sum_path)
    l1_upstream_band = l1_upstream_raster.GetRasterBand(1)
    l1_upstream_nodata = pygeoprocessing.get_nodata_from_uri(
        l1_upstream_sum_path)

    pet_raster = gdal.Open(pet_path)
    pet_band = pet_raster.GetRasterBand(1)

    pygeoprocessing.new_raster_from_base_uri(
        ti_path, subsidized_out_path, 'GTiff', SUBSIDIZED_NODATA,
        gdal.GDT_Int32, fill_value=SUBSIDIZED_NODATA)

    l1_aggregate_values = pygeoprocessing.aggregate_raster_values_uri(
        l1_path, watershed_path)
    l2_aggregate_values = pygeoprocessing.aggregate_raster_values_uri(
        l2_path, watershed_path)

    assert(
        l1_aggregate_values.n_pixels[9999] ==
        l2_aggregate_values.n_pixels[9999])

    for block_info, ti_array in pygeoprocessing.iterblocks(ti_path):
        (x_indexes, y_indexes) = numpy.meshgrid(
            xrange(ti_array.shape[1]), xrange(ti_array.shape[0]))
        index_array = (
            x_indexes + block_info['xoff'] +
            (y_indexes + block_info['yoff']) * n_cols)
        pet_array = pet_band.ReadAsArray(**block_info)
        l2_array = l2_band.ReadAsArray(**block_info)
        l1_array = l1_band.ReadAsArray(**block_info)
        l1_upstream_array = l1_upstream_band.ReadAsArray(
            **block_info)

        valid_mask = (
            (l1_array != L1_NODATA) & (l2_array != L2_NODATA) &
            (ti_array != TI_NODATA) &
            (l1_upstream_array != l1_upstream_nodata) &
            (l1_upstream_array > pet_array))

        array_sorter.append(
            {
                'index_array':  index_array[valid_mask],
                'ti_array': ti_array[valid_mask],
                'l2_array': l2_array[valid_mask],
                'l1_array': l1_array[valid_mask],
            })

    # Initialize the sums
    l1_sum = l1_aggregate_values.total[9999]
    l2_sum = 0.0
    n_pixels = l1_aggregate_values.n_pixels[9999]

    for index_array, ti_array, l1_array, l2_array in zip(
            *[array_sorter.iterarray(_) for _ in [
                'index_array', 'ti_array', 'l1_array', 'l2_array']]):

        # find a cut in l1_array/l2_array that attempts to balance the
        # subsidized area
        low_enough = 0
        too_high = l1_array.size

        # This loop does a binary search for an optimal cut and can handle
        # the case where no elements would add to the solution.
        while True:
            cut_index = (too_high + low_enough) // 2
            running_sum = (
                l1_sum - numpy.sum(l1_array[0:cut_index]) +
                l2_sum + numpy.sum(l2_array[0:cut_index]) - qb_val * n_pixels)
            if running_sum >= 0:
                if too_high - low_enough > 1:
                    low_enough = cut_index
                else:
                    break
            else:
                too_high = cut_index
            if low_enough == too_high:
                cut_index = -1
                break

        if cut_index == -1:
            # there's nothing to add
            break

        l1_sum = l1_sum - numpy.sum(l1_array[0:cut_index])
        l2_sum = l2_sum + numpy.sum(l2_array[0:cut_index])

        subsidized_raster = gdal.Open(subsidized_out_path, gdal.GA_Update)
        subsidized_band = subsidized_raster.GetRasterBand(1)
        for block_info, subsidized_mask_array in (
                pygeoprocessing.iterblocks(subsidized_out_path)):
            (x_indexes, y_indexes) = numpy.meshgrid(
                xrange(subsidized_mask_array.shape[1]),
                xrange(subsidized_mask_array.shape[0]))
            block_indexes = (
                x_indexes + block_info['xoff'] +
                (y_indexes + block_info['yoff']) * n_cols)
            mask = numpy.in1d(
                block_indexes, index_array[0:cut_index]).reshape(
                    block_indexes.shape)
            subsidized_mask_array[mask] = 1
            subsidized_band.WriteArray(
                subsidized_mask_array, xoff=block_info['xoff'],
                yoff=block_info['yoff'])
        subsidized_band.FlushCache()
        subsidized_raster.FlushCache()


def _calculate_upstream_flow(
        flow_direction_path, dem_path, source_path, aoi_path,
        out_upstream_source_path):
    """Calculate upstream flow of source to pixel.

    Parameters:
        flow_direction_path (string): path to a raster that indicates the
            d-infinity flow direction per pixel.
        dem_path (string): path to DEM raster.
        aoi_path (string): path to AOI vector.
        out_upstream_source_path (string): path to output file that contains
            the sum of upstream flow from the `flow_direction_path` raster.

    Returns:
        None.
    """
    zero_absorption_source_path = pygeoprocessing.temporary_filename()
    loss_path = pygeoprocessing.temporary_filename()

    pygeoprocessing.make_constant_raster_from_base_uri(
        dem_path, 0.0, zero_absorption_source_path)

    pygeoprocessing.routing.route_flux(
        flow_direction_path, dem_path, source_path,
        zero_absorption_source_path, loss_path, out_upstream_source_path,
        'flux_only', aoi_uri=aoi_path, include_source=True)

    os.remove(zero_absorption_source_path)
    os.remove(loss_path)


def _calculate_l(l1_path, l2_path, subsidized_mask_path, l_out_path):
    """Calculate L by combining L1, L2, in correct subsidized areas.

    Parameters:
        l1_path (string): path to upland flow
        l2_path (string): path to subsidized flow
        subsidized_mask_path (string): path to raster mask whose pixels are
            1 where subsidized flow should relpace upland
        l_out_path (string): path to output flow.

    Returns:
        None.
    """
    def _l_op(l1_array, l2_array, subsidized_mask_array):
        """Combine L1 into L2 where subsidized area != 1."""
        valid_mask = (l1_array != L1_NODATA) & (l2_array != L2_NODATA)

        result = numpy.empty(l1_array.shape)
        result[:] = L_NODATA
        result[valid_mask] = l1_array[valid_mask]

        l2_mask = (subsidized_mask_array == 1) & valid_mask
        result[l2_mask] = l2_array[l2_mask]
        return result

    cell_size = pygeoprocessing.get_cell_size_from_uri(l1_path)
    pygeoprocessing.vectorize_datasets(
        [l1_path, l2_path, subsidized_mask_path], _l_op, l_out_path,
        gdal.GDT_Float32, L_NODATA, cell_size, 'intersection',
        vectorize_op=False, datasets_are_pre_aligned=True)


class _OutOfCoreNumpyArray(object):
    """Abstraction of a numpy array that can sort out of core."""

    _MAX_SIZE = 2**22

    def __init__(self, working_dir, primary_key):
        """Construct an empty out of core array."""
        self.array_dict = collections.defaultdict(lambda: numpy.array([]))
        self.working_dir = working_dir
        self.primary_key = primary_key
        self.array_files_dict = collections.defaultdict(list)

    def __del__(self):
        for filename_list in self.array_files_dict.itervalues():
            for filename in filename_list:
                os.remove(filename)

    def append(self, array_dict):
        """Append array to current array."""
        first_key = array_dict.iterkeys().next()
        if ((len(self.array_dict[first_key]) != 0) and
                (len(self.array_dict[first_key]) +
                 len(array_dict[first_key])) > self._MAX_SIZE):
            # save off all the current arrays to files
            self._sort()
        for key, array in array_dict.iteritems():
            self.array_dict[key] = numpy.append(array, self.array_dict[key])

    @staticmethod
    def _read_buffer(index_filename, data_filename, chunk_size):
        """Yield values out in sorted order of the primary key.

        Parameters:
            index_filename (string): a path to a file which contains sorted
                indexes.
            data_filename (string): a path to a file of binary floats the same
                length as the file pointed to by
                index_filename.
            chunk_size (int): how many floats to read from the file at once.

        Yields:
            an  (index, value) tuple, where the index can be used to
            globally sort the order in which value should appear.
        """
        current_seek_location = 0
        buffer_size = chunk_size * 4
        while True:
            data_file = open(data_filename, 'rb')
            data_file.seek(current_seek_location, 0)
            data_buffer = data_file.read(buffer_size)
            data_file.close()

            index_file = open(index_filename, 'rb')
            index_file.seek(current_seek_location, 0)
            index_buffer = index_file.read(buffer_size)
            index_file.close()
            if data_buffer == '':
                break
            current_seek_location += len(data_buffer)
            for (index, value) in zip(
                    struct.unpack('f'*(len(index_buffer) / 4), index_buffer),
                    struct.unpack('f'*(len(data_buffer) / 4), data_buffer)):
                # yield a tuple so that the first index is the order
                yield (index, value)

    def iterarray(self, key):
        """Iterate over a keyed array in memory chunks."""
        self._sort()
        chunk_size = 2**12

        array_iterator = heapq.merge(*[
            _OutOfCoreNumpyArray._read_buffer(
                index_filename, data_filename, chunk_size)
            for index_filename, data_filename in zip(
                self.array_files_dict[self.primary_key],
                self.array_files_dict[key])])
        while True:
            # this selects the second value from the iterator which is the
            # value, not the sorted index
            array = numpy.array(
                tuple(
                    [_[1] for _ in itertools.islice(
                        array_iterator, chunk_size)]))
            if len(array) > 0:
                yield array
            else:
                break

    def _sort(self):
        """Out of core argsort on array."""
        if len(self.array_dict[self.primary_key]) == 0:
            return
        argsort = numpy.argsort(-self.array_dict[self.primary_key])
        for key, array in self.array_dict.iteritems():
            # this section ensures that the primary key is storted as it is
            # sorted with a negative above
            if key == self.primary_key:
                scale = -1.0
            else:
                scale = 1.0
            file_path = os.path.join(self.working_dir, str(uuid.uuid4()))
            self.array_files_dict[key].append(file_path)
            with open(file_path, 'wb') as out_file:
                # This sorts and writes the file array in one step
                out_file.write(
                    struct.pack('f'*len(array), *(scale * array[argsort])))
                self.array_dict[key] = numpy.array([])
