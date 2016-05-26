"""InVEST Seasonal Water Yield Model."""

import itertools
import heapq
import struct
import collections
import os
import logging
import re
import fractions
import uuid
import warnings

import scipy.special
import numpy
from osgeo import gdal
from osgeo import ogr
import pygeoprocessing
import pygeoprocessing.routing
import pygeoprocessing.routing.routing_core
from .. import utils

import seasonal_water_yield_core  #pylint: disable=import-error

logging.basicConfig(format='%(asctime)s %(name)-20s %(levelname)-8s \
%(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

LOGGER = logging.getLogger(
    'natcap.invest.seasonal_water_yield.seasonal_water_yield')

_N_MONTHS = 12
MONTH_ID_TO_LABEL = [
    'jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct',
    'nov', 'dec']

_OUTPUT_BASE_FILES = {
    'aggregate_vector_path': 'aggregated_results.shp',
    'annual_precip_path': 'P.tif',
    'cn_path': 'CN.tif',
    'l_avail_path': 'L_avail.tif',
    'l_path': 'L.tif',
    'l_sum_path': 'L_sum.tif',
    'l_sum_avail_path': 'L_sum_avail.tif',
    'qf_path': 'QF.tif',
    'b_sum_path': 'B_sum.tif',
    'b_path': 'B.tif',
    'vri_path': 'Vri.tif',
    }

_INTERMEDIATE_BASE_FILES = {
    'aet_path': 'aet.tif',
    'aetm_path_list': ['aetm_%d.tif' % (_+1) for _ in xrange(_N_MONTHS)],
    'flow_dir_path': 'flow_dir.tif',
    'qfm_path_list': ['qf_%d.tif' % (_+1) for _ in xrange(_N_MONTHS)],
    'stream_path': 'stream.tif',
    'ti_path': 'ti.tif',
}

_TMP_BASE_FILES = {
    'outflow_direction_path': 'outflow_direction.tif',
    'outflow_weights_path': 'outflow_weights.tif',
    'si_path': 'Si.tif',
    'lulc_aligned_path': 'lulc_aligned.tif',
    'dem_aligned_path': 'dem_aligned.tif',
    'lulc_valid_path': 'lulc_valid.tif',
    'dem_valid_path': 'dem_valid.tif',
    'slope_path': 'slope.tif',
    'loss_path': 'loss.tif',
    'soil_depth_aligned_path': 'soil_depth.tif',
    'zero_absorption_source_path': 'zero_absorption.tif',
    'soil_group_aligned_path': 'soil_group_aligned.tif',
    'soil_group_valid_path': 'soil_group_valid.tif',
    'flow_accum_path': 'flow_accum.tif',
    'root_depth_path': 'root_depth.tif',
    'precip_path_aligned_list': [
        'prcp_a%d.tif' % x for x in xrange(_N_MONTHS)],
    'n_events_path_list': ['n_events_%d.tif' % x for x in xrange(_N_MONTHS)],
    'et0_path_aligned_list': ['et0_a_%d.tif' % x for x in xrange(_N_MONTHS)],
    'pet_path_aligned_list': ['pet_%d.tif' % x for x in xrange(_N_MONTHS)],
    'kc_path_list': ['kc_%d.tif' % x for x in xrange(_N_MONTHS)],
    'z_rm_path_list': ['z_rm_%d.tif' % x for x in xrange(_N_MONTHS)],
    'pawc_aligned_path': 'pawc_aligned.tif',
    'wm_path_list': ['wm_%d.tif' % x for x in xrange(_N_MONTHS)],
    'l1_path_list': ['l1_%d.tif' % x for x in xrange(_N_MONTHS)],
    'l_aligned_path': 'l_aligned.tif',
    'cz_aligned_raster_path': 'cz_aligned.tif',
    'subsidized_out_path_list': [
        'subsidized_%d.tif' % x for x in xrange(_N_MONTHS)],
    'temporary_subwatershed_path': 'temporary_subwatershed.shp',
    'l1_upstream_path_list': [
        'l1_upstream_sum_%d.tif' % x for x in xrange(_N_MONTHS)],
    }

ROOT_DEPTH_NODATA = -1.0
N_EVENTS_NODATA = -1
QF_NODATA = -1.0
KC_NODATA = -1.0
SI_NODATA = -1.0
CN_NODATA = -1.0
AET_NODATA = -1.0
L1_NODATA = -1.0
TI_NODATA = -1.0


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
        args['threshold_flow_accumulation'] (number): used when classifying
            stream pixels from the DEM by thresholding the number of upstream
            cells that must flow into a cell before it's considered
            part of a stream.
        args['et0_dir'] (string): required if
            args['user_defined_local_recharge'] is False.  Path to a directory
            that contains rasters of monthly reference evapotranspiration;
            units in mm.
        args['precip_dir'] (string): required if
            args['user_defined_local_recharge'] is False. A path to a
            directory that contains rasters of monthly precipitation; units in
            mm.
        args['dem_raster_path'] (string): a path to a digital elevation
            raster.
        args['pawc_raster_path'] (string): a path to the plant available water
            content raster (PAWC) where each pixel indicates the plant
            available water content as a fraction between 0 and 1.
        args['lulc_raster_path'] (string): a path to a land cover raster used
            to classify biophysical properties of pixels.
        args['soil_group_path'] (string): required if
            args['user_defined_local_recharge'] is  False. A path to a raster
            indicating SCS soil groups where integer values are mapped to soil
            types:
                1: A
                2: B
                3: C
                4: D
        args['soil_depth_raster_path'] (string): path to raster that indicates
            soil depth per pixel in mm.
        args['aoi_path'] (string): path to a vector that indicates the area
            over which the model should be run, as well as the area in which
            to aggregate over when calculating the output Qb.
        args['biophysical_table_path'] (string): path to a CSV table that maps
            landcover codes paired with soil group types to curve numbers as
            well as Kc values.  Headers must include 'lucode', 'CN_A', 'CN_B',
            'CN_C', 'CN_D', 'Kc_1', 'Kc_2', 'Kc_3', 'Kc_4', 'Kc_5', 'Kc_6',
            'Kc_7', 'Kc_8', 'Kc_9', 'Kc_10', 'Kc_11', 'Kc_12'.
        args['rain_events_table_path'] (string): Not required if
            args['user_defined_local_recharge'] is True or
            args['user_defined_climate_zones'] is True.  Path to a CSV table
            that has headers 'month' (1-12) and 'events' (int >= 0) that
            indicates the number of rain events per month
        args['alpha_m'] (float or string): required if args['monthly_alpha']
            is false.  Is the proportion of upslope annual available local
            recharge that is available in month m.
        args['beta_i'] (float or string): is the fraction of the upgradient
            subsidy that is available for downgradient evapotranspiration.
        args['gamma'] (float or string): is the fraction of pixel local
            recharge that is available to downgradient pixels.
        args['user_defined_local_recharge'] (boolean): if True, indicates user
            will provide pre-defined local recharge raster layer
        args['l_path'] (string): required if
            args['user_defined_local_recharge'] is True.  If provided pixels
            indicate the amount of local recharge; units in mm.
        args['user_defined_climate_zones'] (boolean): if True, user provides
            a climate zone rain events table and a climate zone raster map in
            lieu of a global rain events table.
        args['climate_zone_table_path'] (string): required if
            args['user_defined_climate_zones'] is True. Contains monthly
            precipitation events per climate zone.  Fields must be:
            "cz_id", "jan", "feb", "mar", "apr", "may", "jun", "jul",
            "aug", "sep", "oct", "nov", "dec".
        args['climate_zone_raster_path'] (string): required if
            args['user_defined_climate_zones'] is True, pixel values
            correspond to the "cz_id" values defined in
            args['climate_zone_table_path']
        args['monthly_alpha'] (boolean): if True, use the alpha
        args['monthly_alpha_path'] (string): required if args['monthly_alpha']
            is True.
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

    # fail early on a missing required rain events table
    if (not args['user_defined_local_recharge'] and
            not args['user_defined_climate_zones']):
        rain_events_lookup = (
            pygeoprocessing.get_lookup_from_table(
                args['rain_events_table_path'], 'month'))

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

    biophysical_table = pygeoprocessing.get_lookup_from_table(
        args['biophysical_table_path'], 'lucode')

    if args['monthly_alpha']:
        # parse out the alpha lookup table of the form (month_id: alpha_val)
        alpha_month = dict(
            (key, val['alpha']) for key, val in
            pygeoprocessing.get_lookup_from_table(
                args['monthly_alpha_path'], 'month').iteritems())
    else:
        # make all 12 entries equal to args['alpha_m']
        alpha_m = float(fractions.Fraction(args['alpha_m']))
        alpha_month = dict(
            (month_index+1, alpha_m) for month_index in xrange(12))

    beta_i = float(fractions.Fraction(args['beta_i']))
    gamma = float(fractions.Fraction(args['gamma']))
    threshold_flow_accumulation = float(args['threshold_flow_accumulation'])
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

    LOGGER.info('Checking that the AOI is not the output aggregate vector')
    if (os.path.normpath(args['aoi_path']) ==
            os.path.normpath(file_registry['aggregate_vector_path'])):
        raise ValueError(
            "The input AOI is the same as the output aggregate vector, "
            "please choose a different workspace or move the AOI file "
            "out of the current workspace %s" %
            file_registry['aggregate_vector_path'])

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

    if not args['user_defined_local_recharge']:
        precip_path_list = []
        et0_path_list = []

        et0_dir_list = [
            os.path.join(args['et0_dir'], f) for f in os.listdir(
                args['et0_dir'])]
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

        input_align_list = (
            precip_path_list + [args['soil_group_path']] + et0_path_list +
            input_align_list)
        output_align_list = (
            file_registry['precip_path_aligned_list'] +
            [file_registry['soil_group_aligned_path']] +
            file_registry['et0_path_aligned_list'] + output_align_list)

    align_index = len(input_align_list) - 1  # this aligns with the DEM
    if args['user_defined_local_recharge']:
        input_align_list.append(args['l_path'])
        output_align_list.append(file_registry['l_aligned_path'])
    elif args['user_defined_climate_zones']:
        input_align_list.append(args['climate_zone_raster_path'])
        output_align_list.append(
            file_registry['cz_aligned_raster_path'])
    interpolate_list = ['nearest'] * len(input_align_list)

    pygeoprocessing.align_dataset_list(
        input_align_list, output_align_list, interpolate_list, pixel_size,
        'intersection', align_index, aoi_uri=args['aoi_path'],
        assert_datasets_projected=True)

    # sometimes users input data where the DEM is defined in places where the
    # land cover isn't, mask those out
    LOGGER.info("Masking invalid lulc, dem, and possible soil group overlap")
    input_raster_path_list = [
        file_registry['dem_aligned_path'],
        file_registry['lulc_aligned_path']]
    output_valid_raster_path_list = [
        file_registry['dem_valid_path'],
        file_registry['lulc_valid_path']]
    if not args['user_defined_local_recharge']:
        input_raster_path_list.append(file_registry['soil_group_aligned_path'])
        output_valid_raster_path_list.append(
            file_registry['soil_group_valid_path'])
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
    LOGGER.info('classify kc')
    for month_index in xrange(12):
        kc_lookup = dict([
            (lucode, biophysical_table[lucode]['kc_%d' % (month_index+1)])
            for lucode in biophysical_table])
        pygeoprocessing.reclassify_dataset_uri(
            file_registry['lulc_valid_path'], kc_lookup,
            file_registry['kc_path_list'][month_index], gdal.GDT_Float32,
            KC_NODATA)

    # user didn't predefine local recharge so calculate it
    LOGGER.info('loading number of monthly events')
    for month_id in xrange(_N_MONTHS):
        if args['user_defined_climate_zones']:
            cz_rain_events_lookup = (
                pygeoprocessing.get_lookup_from_table(
                    args['climate_zone_table_path'], 'cz_id'))
            month_label = MONTH_ID_TO_LABEL[month_id]
            climate_zone_rain_events_month = dict([
                (cz_id, cz_rain_events_lookup[cz_id][month_label]) for
                cz_id in cz_rain_events_lookup])
            pygeoprocessing.reclassify_dataset_uri(
                file_registry['cz_aligned_raster_path'],
                climate_zone_rain_events_month,
                file_registry['n_events_path_list'][month_id],
                gdal.GDT_Float32, N_EVENTS_NODATA)
        else:
            # rain_events_lookup defined near entry point of execute
            n_events = rain_events_lookup[month_id+1]['events']
            pygeoprocessing.make_constant_raster_from_base_uri(
                file_registry['dem_valid_path'], n_events,
                file_registry['n_events_path_list'][month_id])

    for month_index in xrange(_N_MONTHS):
        LOGGER.info("For month %d: ", month_index)
        _calculate_aet_uphill(
            file_registry['precip_path_aligned_list'][month_index],
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
            file_registry['precip_path_aligned_list'][month_index],
            file_registry['aetm_path_list'][month_index],
            file_registry['l1_path_list'][month_index]
            )

    LOGGER.info('flow direction')
    pygeoprocessing.routing.flow_direction_d_inf(
        file_registry['dem_valid_path'],
        file_registry['flow_dir_path'])

    LOGGER.info('flow weights')
    pygeoprocessing.routing.routing_core.calculate_flow_weights(
        file_registry['flow_dir_path'],
        file_registry['outflow_weights_path'],
        file_registry['outflow_direction_path'])

    LOGGER.info('flow accumulation')
    pygeoprocessing.routing.flow_accumulation(
        file_registry['flow_dir_path'],
        file_registry['dem_valid_path'],
        file_registry['flow_accum_path'])

    LOGGER.info('calculate slope')
    pygeoprocessing.calculate_slope(
        file_registry['dem_valid_path'], file_registry['slope_path'])

    LOGGER.info("Calculate TI")
    _calculate_ti(
        file_registry['flow_accum_path'], file_registry['slope_path'],
        file_registry['soil_depth_aligned_path'], file_registry['ti_path'])

    LOGGER.info("calculate subsidized area")
    for month_index in xrange(_N_MONTHS):
        LOGGER.info("For month %d: ", month_index)

        _calculate_upstream_flow(
            file_registry['flow_dir_path'],
            file_registry['dem_valid_path'],
            file_registry['l1_path_list'][month_index],
            args['aoi_path'],
            file_registry['l1_upstream_path_list'][month_index])
        _calculate_subsidized_area(
            file_registry['l1_path_list'][month_index],
            file_registry['l1_upstream_path_list'][month_index],
            file_registry['pet_path_aligned_list'][month_index],
            file_registry['ti_path'], args['aoi_path'],
            file_registry['temporary_subwatershed_path'],
            file_registry['subsidized_out_path_list'][month_index])
    return


    LOGGER.info('stream thresholding')
    pygeoprocessing.routing.stream_threshold(
        file_registry['flow_accum_path'],
        threshold_flow_accumulation,
        file_registry['stream_path'])

    LOGGER.info('quick flow')
    if args['user_defined_local_recharge']:
        file_registry['l_path'] = file_registry['l_aligned_path']
        l_nodata = pygeoprocessing.get_nodata_from_uri(file_registry['l_path'])

        def l_avail_op(l_array):
            """Calculate equation [8] L_avail = max(gamma*L, 0)."""
            result = numpy.empty(l_array.shape)
            result[:] = l_nodata
            valid_mask = (l_array != l_nodata)
            valid_l_array = l_array[valid_mask]
            valid_l_array[valid_l_array < 0.0] = 0.0
            result[valid_mask] = valid_l_array * gamma
            return result
        pygeoprocessing.vectorize_datasets(
            [file_registry['l_path']], l_avail_op,
            file_registry['l_avail_path'], gdal.GDT_Float32, l_nodata,
            pixel_size, 'intersection', vectorize_op=False,
            datasets_are_pre_aligned=True)
    else:
        # user didn't predefine local recharge so calculate it
        LOGGER.info('loading number of monthly events')
        for month_id in xrange(_N_MONTHS):
            if args['user_defined_climate_zones']:
                cz_rain_events_lookup = (
                    pygeoprocessing.get_lookup_from_table(
                        args['climate_zone_table_path'], 'cz_id'))
                month_label = MONTH_ID_TO_LABEL[month_id]
                climate_zone_rain_events_month = dict([
                    (cz_id, cz_rain_events_lookup[cz_id][month_label]) for
                    cz_id in cz_rain_events_lookup])
                pygeoprocessing.reclassify_dataset_uri(
                    file_registry['cz_aligned_raster_path'],
                    climate_zone_rain_events_month,
                    file_registry['n_events_path_list'][month_id],
                    gdal.GDT_Float32, N_EVENTS_NODATA)
            else:
                # rain_events_lookup defined near entry point of execute
                n_events = rain_events_lookup[month_id+1]['events']
                pygeoprocessing.make_constant_raster_from_base_uri(
                    file_registry['dem_valid_path'], n_events,
                    file_registry['n_events_path_list'][month_id])

        LOGGER.info('calculate curve number')
        _calculate_curve_number_raster(
            file_registry['lulc_valid_path'],
            file_registry['soil_group_aligned_path'],
            biophysical_table, file_registry['cn_path'])

        LOGGER.info('calculate Si raster')
        _calculate_si_raster(
            file_registry['cn_path'], file_registry['stream_path'],
            file_registry['si_path'])

        for month_index in xrange(_N_MONTHS):
            LOGGER.info('calculate quick flow for month %d', month_index+1)
            _calculate_monthly_quick_flow(
                file_registry['precip_path_aligned_list'][month_index],
                file_registry['lulc_valid_path'], file_registry['cn_path'],
                file_registry['n_events_path_list'][month_index],
                file_registry['stream_path'],
                file_registry['qfm_path_list'][month_index],
                file_registry['si_path'])

        LOGGER.info('calculate QFi')

        def qfi_sum_op(*qf_values):
            """Sum the monthly qfis."""
            qf_sum = numpy.zeros(qf_values[0].shape)
            valid_mask = qf_values[0] != QF_NODATA
            valid_qf_sum = qf_sum[valid_mask]
            for index in xrange(len(qf_values)):
                valid_qf_sum += qf_values[index][valid_mask]
            qf_sum[:] = QF_NODATA
            qf_sum[valid_mask] = valid_qf_sum
            return qf_sum

        pygeoprocessing.vectorize_datasets(
            file_registry['qfm_path_list'], qfi_sum_op,
            file_registry['qf_path'], gdal.GDT_Float32, QF_NODATA,
            pixel_size, 'intersection', vectorize_op=False,
            datasets_are_pre_aligned=True)

        LOGGER.info('calculate local recharge')

        # call through to a cython function that does the necessary routing
        # between AET and L.sum.avail in equation [7], [4], and [3]
        seasonal_water_yield_core.calculate_local_recharge(
            file_registry['precip_path_aligned_list'],
            file_registry['et0_path_aligned_list'],
            file_registry['qfm_path_list'],
            file_registry['flow_dir_path'],
            file_registry['outflow_weights_path'],
            file_registry['outflow_direction_path'],
            file_registry['dem_valid_path'],
            file_registry['lulc_valid_path'], alpha_month,
            beta_i, gamma, file_registry['stream_path'],
            file_registry['l_path'],
            file_registry['l_avail_path'],
            file_registry['l_sum_avail_path'],
            file_registry['aet_path'], file_registry['kc_path_list'])

    # calculate Qb as the sum of local_recharge_avail over the AOI, Eq [9]
    qb_sum, qb_valid_count = _sum_valid(file_registry['l_path'])
    qb_result = qb_sum / qb_valid_count

    pixel_size = pygeoprocessing.get_cell_size_from_uri(
        file_registry['l_path'])
    ri_nodata = pygeoprocessing.get_nodata_from_uri(
        file_registry['l_path'])

    def vri_op(ri_array):
        """Calculate vri index [Eq 10]."""
        return numpy.where(
            ri_array != ri_nodata,
            ri_array / qb_result / qb_valid_count, ri_nodata)
    pygeoprocessing.vectorize_datasets(
        [file_registry['l_path']], vri_op,
        file_registry['vri_path'], gdal.GDT_Float32, ri_nodata,
        pixel_size, 'intersection', vectorize_op=False,
        datasets_are_pre_aligned=True)

    _aggregate_recharge(
        args['aoi_path'], file_registry['l_path'],
        file_registry['vri_path'],
        file_registry['aggregate_vector_path'])

    LOGGER.info('calculate L_sum')  # Eq. [12]
    pygeoprocessing.make_constant_raster_from_base_uri(
        file_registry['dem_valid_path'], 0.0,
        file_registry['zero_absorption_source_path'])
    pygeoprocessing.routing.route_flux(
        file_registry['flow_dir_path'],
        file_registry['dem_valid_path'],
        file_registry['l_path'],
        file_registry['zero_absorption_source_path'],
        file_registry['loss_path'],
        file_registry['l_sum_path'], 'flux_only',
        stream_uri=file_registry['stream_path'])

    LOGGER.info('calculate B_sum')
    seasonal_water_yield_core.route_baseflow_sum(
        file_registry['dem_valid_path'],
        file_registry['l_path'],
        file_registry['l_avail_path'],
        file_registry['l_sum_path'],
        file_registry['outflow_direction_path'],
        file_registry['outflow_weights_path'],
        file_registry['stream_path'],
        file_registry['b_sum_path'])

    LOGGER.info('calculate B')
    b_sum_nodata = pygeoprocessing.get_nodata_from_uri(
        file_registry['b_sum_path'])
    ri_nodata = b_sum_nodata

    def op_b(b_sum, l_avail, l_sum):
        """Calculate B=B_sum*Lavail/L_sum."""
        valid_mask = ((b_sum != b_sum_nodata) & (l_sum != 0))
        result = numpy.empty(b_sum.shape)
        result[:] = b_sum_nodata
        result[valid_mask] = (
            b_sum[valid_mask] * l_avail[valid_mask] / l_sum[valid_mask])
        return result

    pygeoprocessing.vectorize_datasets(
        [file_registry['b_sum_path'],
         file_registry['l_avail_path'],
         file_registry['l_sum_path']], op_b,
        file_registry['b_path'],
        gdal.GDT_Float32, b_sum_nodata, pixel_size, 'intersection',
        vectorize_op=False, datasets_are_pre_aligned=True)

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


def _calculate_monthly_quick_flow(
        precip_path, lulc_raster_path, cn_path, n_events_raster_path,
        stream_path, qf_monthly_path, si_path):
    """Calculate quick flow for a month.

    Parameters:
        precip_path (string): path to file that correspond to monthly
            precipitation
        lulc_raster_path (string): path to landcover raster
        cn_path (string): path to curve number raster
        n_events_raster_path (string): a path to a raster where each pixel
            indicates the number of rain events.
        stream_path (string): path to stream mask raster where 1 indicates a
            stream pixel, 0 is a non-stream but otherwise valid area from the
            original DEM, and nodata indicates areas outside the valid DEM.
        qf_monthly_path_list (list of string): list of paths to output monthly
            rasters.
        si_path (string): list to output raster for potential maximum retention

    Returns:
        None
    """
    cn_nodata = pygeoprocessing.get_nodata_from_uri(cn_path)

    def si_op(ci_array, stream_array):
        """Potential maximum retention."""
        si_array = 1000.0 / ci_array - 10
        si_array = numpy.where(ci_array != cn_nodata, si_array, SI_NODATA)
        si_array[stream_array == 1] = 0
        return si_array

    pixel_size = pygeoprocessing.get_cell_size_from_uri(
        lulc_raster_path)
    pygeoprocessing.vectorize_datasets(
        [cn_path, stream_path], si_op, si_path, gdal.GDT_Float32,
        SI_NODATA, pixel_size, 'intersection', vectorize_op=False,
        datasets_are_pre_aligned=True)

    p_nodata = pygeoprocessing.get_nodata_from_uri(precip_path)
    n_events_nodata = pygeoprocessing.get_nodata_from_uri(n_events_raster_path)

    def qf_op(p_im, s_i, n_events, stream_array):
        """Calculate quick flow as in Eq [1] in user's guide.

        Parameters:
            p_im (numpy.array): precipitation at pixel i on month m
            s_i (numpy.array): factor that is 1000/CN_i - 10
                (Equation 1b from user's guide)
            n_events (numpy.array): number of rain events on the pixel
            stream_mask (numpy.array): 1 if stream, otherwise not a stream
                pixel.

        Returns:
            quick flow (numpy.array)
        """
        valid_mask = (
            (p_im != p_nodata) & (s_i != SI_NODATA) & (p_im != 0.0) &
            (stream_array != 1) & (n_events != n_events_nodata) &
            (n_events > 0))
        valid_n_events = n_events[valid_mask]
        valid_si = s_i[valid_mask]

        # a_im is the mean rain depth on a rainy day at pixel i on month m
        # the 25.4 converts inches to mm since Si is in inches
        a_im = numpy.empty(valid_n_events.shape)
        a_im = p_im[valid_mask] / valid_n_events / 25.4
        qf_im = numpy.empty(p_im.shape)
        qf_im[:] = QF_NODATA

        # Precompute the last two terms in quickflow so we can handle a
        # numerical instability when s_i is large and/or a_im is small
        # on large valid_si/a_im this number will be zero and the latter
        # exponent will also be zero because of a divide by zero. rather than
        # raise that numerical warning, just handle it manually
        E1 = scipy.special.expn(1, valid_si / a_im)  #pylint: disable=invalid-name,no-member
        nonzero_e1_mask = E1 != 0
        exp_result = numpy.zeros(valid_si.shape)
        exp_result[nonzero_e1_mask] = numpy.exp(
            (0.8 * valid_si[nonzero_e1_mask]) / a_im[nonzero_e1_mask] +
            numpy.log(E1[nonzero_e1_mask]))

        # qf_im is the quickflow at pixel i on month m Eq. [1]
        qf_im[valid_mask] = (25.4 * valid_n_events * (
            (a_im - valid_si) * numpy.exp(-0.2 * valid_si / a_im) +
            valid_si ** 2 / a_im * exp_result))

        # if precip is 0, then QF should be zero
        qf_im[(p_im == 0) | (n_events == 0)] = 0.0
        # if we're on a stream, set quickflow to the precipitation
        qf_im[stream_array == 1] = p_im[stream_array == 1]
        return qf_im

    pygeoprocessing.vectorize_datasets(
        [precip_path, si_path, n_events_raster_path, stream_path], qf_op,
        qf_monthly_path, gdal.GDT_Float32, QF_NODATA, pixel_size,
        'intersection', vectorize_op=False, datasets_are_pre_aligned=True)


def _calculate_curve_number_raster(
        lulc_raster_path, soil_group_path, biophysical_table, cn_path):
    """Calculate the CN raster from the landcover and soil group rasters.

    Parameters:
        lulc_raster_path (string): path to landcover raster
        soil_group_path (string): path to raster indicating soil group where
            pixel values are in [1,2,3,4]
        biophysical_table (dict): maps landcover IDs to dictionaries that
            contain at least the keys 'cn_a', 'cn_b', 'cn_c', 'cn_d', that
            map to the curve numbers for that landcover and soil type.
        cn_path (string): path to output curve number raster to be output
            which will be the dimensions of the intersection of
            `lulc_raster_path` and `soil_group_path` the cell size of
            `lulc_raster_path`.

    Returns:
        None
    """
    soil_nodata = pygeoprocessing.get_nodata_from_uri(soil_group_path)
    map_soil_type_to_header = {
        1: 'cn_a',
        2: 'cn_b',
        3: 'cn_c',
        4: 'cn_d',
    }
    lulc_to_soil = {}
    lulc_nodata = pygeoprocessing.get_nodata_from_uri(lulc_raster_path)
    for soil_id, soil_column in map_soil_type_to_header.iteritems():
        lulc_to_soil[soil_id] = {
            'lulc_values': [],
            'cn_values': []
        }
        for lucode in sorted(biophysical_table.keys() + [lulc_nodata]):
            if lucode != lulc_nodata:
                lulc_to_soil[soil_id]['cn_values'].append(
                    biophysical_table[lucode][soil_column])
                lulc_to_soil[soil_id]['lulc_values'].append(lucode)
            else:
                # handle the lulc nodata with cn nodata
                lulc_to_soil[soil_id]['lulc_values'].append(lulc_nodata)
                lulc_to_soil[soil_id]['cn_values'].append(CN_NODATA)

        # Making the array an int64 to make sure it's big enough to handle
        # both signed and unsigned int32 values
        lulc_to_soil[soil_id]['lulc_values'] = (
            numpy.array(lulc_to_soil[soil_id]['lulc_values'],
                        dtype=numpy.int64))
        lulc_to_soil[soil_id]['cn_values'] = (
            numpy.array(lulc_to_soil[soil_id]['cn_values'],
                        dtype=numpy.float32))

    def cn_op(lulc_array, soil_group_array):
        """Map lulc code and soil to a curve number."""
        cn_result = numpy.empty(lulc_array.shape)
        cn_result[:] = CN_NODATA
        for soil_group_id in numpy.unique(soil_group_array):
            if soil_group_id == soil_nodata:
                continue
            current_soil_mask = (soil_group_array == soil_group_id)
            index = numpy.digitize(
                lulc_array.ravel(),
                lulc_to_soil[soil_group_id]['lulc_values'], right=True)
            cn_values = (
                lulc_to_soil[soil_group_id]['cn_values'][index]).reshape(
                    lulc_array.shape)
            cn_result[current_soil_mask] = cn_values[current_soil_mask]
        return cn_result

    pixel_size = pygeoprocessing.get_cell_size_from_uri(lulc_raster_path)
    pygeoprocessing.vectorize_datasets(
        [lulc_raster_path, soil_group_path], cn_op, cn_path, gdal.GDT_Float32,
        CN_NODATA, pixel_size, 'intersection', vectorize_op=False,
        datasets_are_pre_aligned=True)


def _calculate_si_raster(cn_path, stream_path, si_path):
    """Calculate the S factor of the quickflow equation [1].

    Parameters:
        cn_path (string): path to curve number raster
        stream_path (string): path to a stream raster (0, 1)
        si_path (string): path to output s_i raster

    Returns:
        None
    """
    cn_nodata = pygeoprocessing.get_nodata_from_uri(cn_path)

    def si_op(ci_factor, stream_mask):
        """Calculate si factor."""
        valid_mask = (ci_factor != cn_nodata) & (ci_factor > 0)
        si_array = numpy.empty(ci_factor.shape)
        si_array[:] = SI_NODATA
        # multiply by the stream mask != 1 so we get 0s on the stream and
        # unaffected results everywhere else
        si_array[valid_mask] = (
            (1000.0 / ci_factor[valid_mask] - 10) * (
                stream_mask[valid_mask] != 1))
        return si_array

    pixel_size = pygeoprocessing.get_cell_size_from_uri(cn_path)
    pygeoprocessing.vectorize_datasets(
        [cn_path, stream_path], si_op, si_path, gdal.GDT_Float32,
        SI_NODATA, pixel_size, 'intersection', vectorize_op=False,
        datasets_are_pre_aligned=True)


def _aggregate_recharge(
        aoi_path, l_path, vri_path, aggregate_vector_path):
    """Aggregate recharge values for the provided watersheds/AOIs.

    Generates a new shapefile that's a copy of 'aoi_path' in sum values from L
    and Vri.

    Parameters:
        aoi_path (string): path to shapefile that will be used to
            aggregate rasters
        l_path (string): path to (L) local recharge raster
        vri_path (string): path to Vri raster
        aggregate_vector_path (string): path to shapefile that will be created
            by this function as the aggregating output.  will contain fields
            'l_sum' and 'vri_sum' per original feature in `aoi_path`.  If this
            file exists on disk prior to the call it is overwritten with
            the result of this call.

    Returns:
        None
    """
    if os.path.exists(aggregate_vector_path):
        LOGGER.warn(
            '%s exists, deleting and writing new output',
            aggregate_vector_path)
        os.remove(aggregate_vector_path)

    esri_driver = ogr.GetDriverByName('ESRI Shapefile')
    original_aoi_vector = ogr.Open(aoi_path)

    esri_driver.CopyDataSource(
        original_aoi_vector, aggregate_vector_path)
    esri_driver = None
    ogr.DataSource.__swig_destroy__(original_aoi_vector)
    original_aoi_vector = None
    aggregate_vector = ogr.Open(aggregate_vector_path, 1)
    aggregate_layer = aggregate_vector.GetLayer()

    # make an identifying id per polygon that can be used for aggregation
    while True:
        serviceshed_defn = aggregate_layer.GetLayerDefn()
        poly_id_field = str(uuid.uuid4())[-8:]
        if serviceshed_defn.GetFieldIndex(poly_id_field) == -1:
            break
    layer_id_field = ogr.FieldDefn(poly_id_field, ogr.OFTInteger)
    aggregate_layer.CreateField(layer_id_field)
    for poly_index, poly_feat in enumerate(aggregate_layer):
        poly_feat.SetField(poly_id_field, poly_index)
        aggregate_layer.SetFeature(poly_feat)
    aggregate_layer.SyncToDisk()

    for raster_path, aggregate_field_id, op_type in [
            (l_path, 'qb', 'mean'), (vri_path, 'vri_sum', 'sum')]:

        # aggregate carbon stocks by the new ID field
        aggregate_stats = pygeoprocessing.aggregate_raster_values_uri(
            raster_path, aggregate_vector_path,
            shapefile_field=poly_id_field, ignore_nodata=True,
            threshold_amount_lookup=None, ignore_value_list=[],
            process_pool=None, all_touched=False)

        aggregate_field = ogr.FieldDefn(aggregate_field_id, ogr.OFTReal)
        aggregate_layer.CreateField(aggregate_field)

        aggregate_layer.ResetReading()
        for poly_index, poly_feat in enumerate(aggregate_layer):
            if op_type == 'mean':
                value = (aggregate_stats.total[poly_index] /
                         aggregate_stats.n_pixels[poly_index])
            elif op_type == 'sum':
                value = aggregate_stats.total[poly_index]
            poly_feat.SetField(aggregate_field_id, value)
            aggregate_layer.SetFeature(poly_feat)

    # don't need a random poly id anymore
    aggregate_layer.DeleteField(
        serviceshed_defn.GetFieldIndex(poly_id_field))
    aggregate_layer.SyncToDisk()
    aggregate_layer = None
    ogr.DataSource.__swig_destroy__(aggregate_vector)
    aggregate_vector = None


def _sum_valid(raster_path):
    """Calculate the sum of the non-nodata pixels in the raster.

    Parameters:
        raster_path (string): path to raster on disk

    Returns:
        (sum, n_pixels) tuple where sum is the sum of the non-nodata pixels
        and n_pixels is the count of them
    """
    raster_sum = 0
    raster_count = 0
    raster_nodata = pygeoprocessing.get_nodata_from_uri(raster_path)

    for _, block in pygeoprocessing.iterblocks(raster_path, band_list=[1]):
        valid_mask = block != raster_nodata
        raster_sum += numpy.sum(block[valid_mask])
        raster_count += numpy.count_nonzero(valid_mask)
    return raster_sum, raster_count


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
    pet_nodata = et0_nodata  # reasonable to set PET nodata to et0 nodata
    pixel_size = pygeoprocessing.get_cell_size_from_uri(precip_path)

    def _pet_op(et0_array, kc_array):
        """Calculate PET."""
        result = numpy.empty(et0_array.shape, dtype=numpy.float32)
        result[:] = pet_nodata
        valid_mask = ((et0_array != et0_nodata) & (kc_array != kc_nodata))
        result[valid_mask] = et0_array[valid_mask] * kc_array[valid_mask]
        return result

    pygeoprocessing.vectorize_datasets(
        [et0_path, kc_path], _pet_op, out_pet_path,
        gdal.GDT_Float32, pet_nodata, pixel_size, 'intersection',
        vectorize_op=False, datasets_are_pre_aligned=True)

    LOGGER.info("calculate Wm")
    precip_nodata = pygeoprocessing.get_nodata_from_uri(precip_path)

    def _w_m_op(precip_array, pet_array):
        """Calculate W_m in equation 2a."""
        result = numpy.empty(precip_array.shape, dtype=numpy.float32)
        result[:] = AET_NODATA
        valid_mask = (
            (precip_array != precip_nodata) &
            (pet_array != pet_nodata) &
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
            (pet_array != pet_nodata) &
            (root_depth_array != root_depth_nodata) &
            (w_m_array != AET_NODATA) &
            (z_rm_array != AET_NODATA))

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
        l1_path, l1_upstream_sum_path, pet_path, ti_path, subwatershed_path,
        temporary_subwatershed_path, subsidized_out_path):
    """Calculated subsidized area such that Eq. 4 is balanced.

    Parameters:
        l1_path (string): path to flow raster for upland regions.
        l1_upstream_sum_path (string): path to sum of L1 values upstream.
        pet_path (string): path to PET raster to use as a -1 * for flow
            raster for subsidized regions.
        ti_path (string): path to topographical index raster.
        subwatershed_path (string): path to a shapefile that defines the
            subwatersheds.
        temporary_subwatershed_path (string): path to a shapefile that can
            be used to mask out subsets of the subwatershed; this file
            will be overwritten if exists and deleted at the end of the
            call.
        subsidized_out_path (string): path to output raster that masks out
            the subsidized region.

    Returns:
        None.
    """
    (_, n_cols) = pygeoprocessing.get_row_col_from_uri(l1_path)
    array_sorter = _OutOfCoreNumpyArray(
        os.path.dirname(subsidized_out_path), 'ti_array')

    l1_raster = gdal.Open(l1_path)
    pet_raster = gdal.Open(pet_path)
    l1_band = l1_raster.GetRasterBand(1)
    pet_band = pet_raster.GetRasterBand(1)

    l1_upstream_raster = gdal.Open(l1_upstream_sum_path)
    l1_upstream_band = l1_upstream_raster.GetRasterBand(1)

    l1_nodata = pygeoprocessing.get_nodata_from_uri(l1_path)
    l1_upstream_nodata = pygeoprocessing.get_nodata_from_uri(
        l1_upstream_sum_path)
    pet_nodata = pygeoprocessing.get_nodata_from_uri(pet_path)
    ti_nodata = pygeoprocessing.get_nodata_from_uri(ti_path)

    pygeoprocessing.new_raster_from_base_uri(
        ti_path, subsidized_out_path, 'GTiff', TI_NODATA, gdal.GDT_Int32,
        fill_value=TI_NODATA)

    # create a temporary shapefile that can be used to rasterize each
    # subwatershed
    # make a shapefile that non-overlapping layers can be added to
    driver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(temporary_subwatershed_path):
        os.remove(temporary_subwatershed_path)
    temporary_subwatershed_vector = driver.CreateDataSource(
        temporary_subwatershed_path)
    spat_ref = pygeoprocessing.get_spatial_ref_uri(subwatershed_path)
    subset_layer = temporary_subwatershed_vector.CreateLayer(
        'subset_layer', spat_ref, ogr.wkbPolygon)
    subwatershed_vector = ogr.Open(subwatershed_path)
    subwatershed_layer = subwatershed_vector.GetLayer()
    defn = subwatershed_layer.GetLayerDefn()

    # For every field, create a duplicate field and add it to the new
    # subset_layer layer
    defn.GetFieldCount()
    for fld_index in range(defn.GetFieldCount()):
        original_field = defn.GetFieldDefn(fld_index)
        output_field = ogr.FieldDefn(
            original_field.GetName(), original_field.GetType())
        subset_layer.CreateField(output_field)

    for feature in subwatershed_layer:
        subset_layer.CreateFeature(feature)
        subset_layer.SyncToDisk()

        l1_sum = pygeoprocessing.aggregate_raster_values_uri(
            l1_path, subwatershed_path).total[9999]
        l2_sum = -pygeoprocessing.aggregate_raster_values_uri(
            pet_path, subwatershed_path).total[9999]

        LOGGER.debug(l1_sum)
        LOGGER.debug(l2_sum)

        mask_path = "%s.tif" % feature.GetFID()
        pygeoprocessing.new_raster_from_base_uri(
            ti_path, mask_path, 'GTiff', TI_NODATA, gdal.GDT_Int32,
            fill_value=TI_NODATA)
        pygeoprocessing.rasterize_layer_uri(
            mask_path, temporary_subwatershed_path, burn_values=[1])

        subset_layer.DeleteFeature(feature.GetFID())

        watershed_mask_raster = gdal.Open(mask_path)
        watershed_mask_band = watershed_mask_raster.GetRasterBand(1)

        for block_info, ti_array in pygeoprocessing.iterblocks(ti_path):
            (x_indexes, y_indexes) = numpy.meshgrid(
                xrange(ti_array.shape[1]), xrange(ti_array.shape[0]))
            index_array = (
                x_indexes + block_info['xoff'] +
                (y_indexes + block_info['yoff']) * n_cols)
            pet_array = pet_band.ReadAsArray(**block_info)
            l1_array = l1_band.ReadAsArray(**block_info)
            watershed_mask_array = watershed_mask_band.ReadAsArray(
                **block_info)
            l1_upstream_array = l1_upstream_band.ReadAsArray(
                **block_info)

            valid_mask = (
                (l1_array != l1_nodata) &
                (pet_array != pet_nodata) &
                (ti_array != ti_nodata) &
                (watershed_mask_array != TI_NODATA) &
                (l1_upstream_array != l1_upstream_nodata) &
                (l1_upstream_array > pet_array))

            array_sorter.append(
                {
                    'index_array':  index_array[valid_mask],
                    'ti_array': ti_array[valid_mask],
                    'pet_array': pet_array[valid_mask],
                    'l1_array': l1_array[valid_mask],
                })

        for sorted_indexes in array_sorter.iterarray('index_array'):
            subsidized_raster = gdal.Open(subsidized_out_path, gdal.GA_Update)
            subsidized_band = subsidized_raster.GetRasterBand(1)
            for block_info, subsidized_mask_array in (
                    pygeoprocessing.iterblocks(subsidized_out_path)):
                (x_indexes, y_indexes) = numpy.meshgrid(
                    xrange(subsidized_mask_array.shape[1]),
                    xrange(subsidized_mask_array.shape[0]))
                index_array = (
                    x_indexes + block_info['xoff'] +
                    (y_indexes + block_info['yoff']) * n_cols)
                mask = numpy.in1d(index_array, sorted_indexes).reshape(
                    index_array.shape)
                subsidized_mask_array[mask] = 1
                subsidized_band.WriteArray(
                    subsidized_mask_array, xoff=block_info['xoff'],
                    yoff=block_info['yoff'])
            subsidized_band.FlushCache()
            subsidized_raster.FlushCache()
            subsidized_band = None
            subsidized_raster = None

        LOGGER.error(
            "breaking because we don't handle overlapping watersheds")
        break

    return


class _OutOfCoreNumpyArray(object):
    """Abstraction of a numpy array that can sort out of core."""

    _MAX_SIZE = 40000

    def __init__(self, working_dir, primary_key):
        """Construct an empty out of core array."""
        self.array_dict = collections.defaultdict(lambda: numpy.array([]))
        self.working_dir = working_dir
        self.primary_key = primary_key

        # keep track of sorted files
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
        chunk_size = 1000

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
