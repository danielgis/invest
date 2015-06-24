"""Module for the execution of the biophysical component of the InVEST Nutrient
Deposition model."""

import logging
import os
import shutil
import math

from osgeo import gdal
from osgeo import ogr
import numpy

import pygeoprocessing.geoprocessing
import pygeoprocessing.routing
import pygeoprocessing.routing.routing_core


LOGGER = logging.getLogger('nutrient')
logging.basicConfig(format='%(asctime)s %(name)-15s %(levelname)-8s \
    %(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')


def execute(args):
    """
        Nutrient delivery ratio model:

        args - a python dictionary with the following entries:
            'workspace_dir' - a string uri pointing to the current workspace.
            'dem_uri' - a string uri pointing to the Digital Elevation Map
                (DEM), a GDAL raster on disk.
            'pixel_yield_uri' - a string uri pointing to the water yield raster
                output from the InVEST Water Yield model.
            'lulc_uri' - a string uri pointing to the landcover GDAL raster.
            'watersheds_uri' - a string uri pointing to an OGR shapefile on
                disk representing the user's watersheds.
            'biophysical_table_uri' - a string uri to a supported table on disk
                containing nutrient retention values. (SAY WHAT VALUES ARE)
            'calc_p' - True if phosphorous is meant to be modeled, if True then
                biophyscial table must have p fields in them.
            'calc_n' - True if nitrogen is meant to be modeled, if True then
                biophyscial table must have n fields in them.
            'subsurface_critical_length_n' - the subsurface flow critical length
                for nitrogen
            'subsurface_critical_length_p' - the subsurface flow critical length
                for phosphorous
            'subsurface_eff_n' - the maximum retention efficiency that soil can
                reach for nitrogen
            'subsurface_eff_p' - the maximum retention efficiency that soil can
                reach for phosphorous
            'results_suffix' - (optional) a text field to append to all output files.
            'threshold_flow_accumulation' - a number representing the flow accumulation.
            '_prepare' - (optional) The preprocessed set of data created by the
                ndr._prepare call.  This argument could be used in cases where the
                call to this function is scripted and can save a significant amount
                of runtime.

        returns nothing.
    """

    def _validate_inputs(nutrients_to_process, lucode_to_parameters):
        """Validation helper method to check that table headers are included
            that are necessary depending on the nutrient type requested by
            the user"""

        #Make sure all the nutrient inputs are good
        if len(nutrients_to_process) == 0:
            raise ValueError("Neither phosphorous nor nitrogen was selected"
                             " to be processed.  Choose at least one.")

        #Build up a list that'll let us iterate through all the input tables
        #and check for the required rows, and report errors if something
        #is missing.
        row_header_table_list = []

        lu_parameter_row = lucode_to_parameters.values()[0]
        row_header_table_list.append(
            (lu_parameter_row, ['load_', 'eff_', 'crit_len_'],
             args['biophysical_table_uri']))

        missing_headers = []
        for row, header_prefixes, table_type in row_header_table_list:
            for nutrient_id in nutrients_to_process:
                for header_prefix in header_prefixes:
                    header = header_prefix + nutrient_id
                    if header not in row:
                        missing_headers.append(
                            "Missing header %s from %s" % (header, table_type))

        if len(missing_headers) > 0:
            raise ValueError('\n'.join(missing_headers))


    if not args['calc_p'] and not args['calc_n']:
        raise Exception('Neither "Calculate Nitrogen" nor "Calculate Phosporus" is selected.  At least one must be selected.')

    #Load all the tables for preprocessing
    workspace = args['workspace_dir']
    output_dir = os.path.join(workspace, 'output')
    intermediate_dir = os.path.join(workspace, 'intermediate')

    try:
        file_suffix = args['results_suffix']
        if not file_suffix.startswith('_'):
            file_suffix = '_' + file_suffix
    except KeyError:
        file_suffix = ''

    for folder in [workspace, output_dir, intermediate_dir]:
        if not os.path.exists(folder):
            os.makedirs(folder)

    #Build up a list of nutrients to process based on what's checked on
    nutrients_to_process = []
    for nutrient_id in ['n', 'p']:
        if args['calc_' + nutrient_id]:
            nutrients_to_process.append(nutrient_id)
    lucode_to_parameters = pygeoprocessing.geoprocessing.get_lookup_from_csv(
        args['biophysical_table_uri'], 'lucode')

    _validate_inputs(nutrients_to_process, lucode_to_parameters)

    if '_prepare' in args:
        preprocessed_data = args['_prepare']
    else:
        preprocessed_data = _prepare(**args)

    aligned_dem_uri = preprocessed_data['aligned_dem_uri']
    thresholded_slope_uri = preprocessed_data['thresholded_slope_uri']
    flow_accumulation_uri = preprocessed_data['flow_accumulation_uri']
    flow_direction_uri = preprocessed_data['flow_direction_uri']

    dem_pixel_size = pygeoprocessing.geoprocessing.get_cell_size_from_uri(
        args['dem_uri'])
    #Pixel size is in m^2, so square and divide by 10000 to get cell size in Ha
    cell_area_ha = dem_pixel_size ** 2 / 10000.0
    out_pixel_size = dem_pixel_size

    #Align all the input rasters
    dem_uri = pygeoprocessing.geoprocessing.temporary_filename()
    lulc_uri = pygeoprocessing.geoprocessing.temporary_filename()
    pygeoprocessing.geoprocessing.align_dataset_list(
        [args['dem_uri'], args['lulc_uri']],
        [dem_uri, lulc_uri], ['nearest'] * 2,
        out_pixel_size, 'intersection', dataset_to_align_index=0,
        aoi_uri=args['watersheds_uri'])

    nodata_landuse = pygeoprocessing.geoprocessing.get_nodata_from_uri(lulc_uri)
    nodata_load = -1.0

    #classify streams from the flow accumulation raster
    LOGGER.info("Classifying streams from flow accumulation raster")
    stream_uri = os.path.join(intermediate_dir, 'stream%s.tif' % file_suffix)
    pygeoprocessing.routing.stream_threshold(
        flow_accumulation_uri,
        float(args['threshold_flow_accumulation']), stream_uri)
    nodata_stream = pygeoprocessing.geoprocessing.get_nodata_from_uri(stream_uri)

    def map_load_function(load_type):
        """Function generator to map arbitrary nutrient type"""
        def map_load(lucode_array):
            """converts unit load to total load & handles nodata"""
            result = numpy.empty(lucode_array.shape)
            result[:] = nodata_load
            for lucode in numpy.unique(lucode_array):
                if lucode != nodata_landuse:
                    result[lucode_array == lucode] = (
                        lucode_to_parameters[lucode][load_type] * cell_area_ha)
            return result
        return map_load
    def map_eff_function(load_type):
        """Function generator to map arbitrary efficiency type"""
        def map_eff(lucode_array, stream_array):
            """maps efficiencies from lulcs, handles nodata, and is aware that
                streams have no retention"""
            result = numpy.empty(lucode_array.shape, dtype=numpy.float32)
            result[:] = nodata_load
            for lucode in numpy.unique(lucode_array):
                if lucode == nodata_landuse:
                    continue
                mask = (lucode_array == lucode) & (stream_array != nodata_stream)
                result[mask] = lucode_to_parameters[lucode][load_type] * (1 - stream_array[mask])
            return result
        return map_eff

    #Build up the load and efficiency rasters from the landcover map
    load_uri = {}
    load_subsurface_uri = {}
    eff_uri = {}
    for nutrient in nutrients_to_process:
        load_uri[nutrient] = os.path.join(
            intermediate_dir, 'load_%s%s.tif' % (nutrient, file_suffix))
        pygeoprocessing.geoprocessing.vectorize_datasets(
            [lulc_uri], map_load_function('load_%s' % nutrient),
            load_uri[nutrient], gdal.GDT_Float32, nodata_load, out_pixel_size,
            "intersection", vectorize_op=False)

        load_subsurface_uri[nutrient] = os.path.join(
            intermediate_dir, 'load_subsurface_%s%s.tif' % (nutrient, file_suffix))
        #pygeoprocessing.geoprocessing.vectorize_datasets(
        #    [lulc_uri], map_load_function('load_subsurface_%s' % nutrient),
        #    load_subsurface_uri[nutrient], gdal.GDT_Float32, nodata_load, out_pixel_size,
        #    "intersection", vectorize_op=False)

        eff_uri[nutrient] = os.path.join(
            intermediate_dir, 'eff_%s%s.tif' % (nutrient, file_suffix))
        pygeoprocessing.geoprocessing.vectorize_datasets(
            [lulc_uri, stream_uri], map_eff_function('eff_%s' % nutrient),
            eff_uri[nutrient], gdal.GDT_Float32, nodata_load, out_pixel_size,
            "intersection", vectorize_op=False)

    field_summaries = {}
    field_header_order = []

    watershed_output_datasource_uri = os.path.join(
        output_dir, 'watershed_results_ndr%s.shp' % file_suffix)
    #If there is already an existing shapefile with the same name and path,
    #delete it then copy the input shapefile into the designated output folder
    if os.path.isfile(watershed_output_datasource_uri):
        os.remove(watershed_output_datasource_uri)
    esri_driver = ogr.GetDriverByName('ESRI Shapefile')
    original_datasource = ogr.Open(args['watersheds_uri'])
    output_datasource = esri_driver.CopyDataSource(
        original_datasource, watershed_output_datasource_uri)
    output_layer = output_datasource.GetLayer()

    add_fields_to_shapefile('ws_id', field_summaries, output_layer, field_header_order)
    field_header_order = []

    export_uri = {}
    field_summaries = {}

    #Calculate the W factor
    LOGGER.info('calculate per pixel W')
    original_w_factor_uri = os.path.join(
        intermediate_dir, 'w_factor%s.tif' % file_suffix)
    thresholded_w_factor_uri = os.path.join(
        intermediate_dir, 'thresholded_w_factor%s.tif' % file_suffix)

    #map lulc to biophysical table
    lulc_to_c = dict([
        (lulc_code, float(table['usle_c'])) for
        (lulc_code, table) in lucode_to_parameters.items()])
    lulc_nodata = pygeoprocessing.geoprocessing.get_nodata_from_uri(lulc_uri)
    w_nodata = -1.0

    pygeoprocessing.geoprocessing.reclassify_dataset_uri(
        lulc_uri, lulc_to_c, original_w_factor_uri, gdal.GDT_Float32,
        w_nodata, exception_flag='values_required')
    def threshold_w(w_val):
        '''Threshold w to 0.001'''
        w_val_copy = w_val.copy()
        nodata_mask = w_val == w_nodata
        w_val_copy[w_val < 0.001] = 0.001
        w_val_copy[nodata_mask] = w_nodata
        return w_val_copy
    pygeoprocessing.geoprocessing.vectorize_datasets(
        [original_w_factor_uri], threshold_w, thresholded_w_factor_uri,
        gdal.GDT_Float32, w_nodata, out_pixel_size, "intersection",
        dataset_to_align_index=0, vectorize_op=False)

    #calculate W_bar
    zero_absorption_source_uri = pygeoprocessing.geoprocessing.temporary_filename()
    loss_uri = pygeoprocessing.geoprocessing.temporary_filename()
    #need this for low level route_flux function
    pygeoprocessing.geoprocessing.make_constant_raster_from_base_uri(
        aligned_dem_uri, 0.0, zero_absorption_source_uri)

    flow_accumulation_nodata = pygeoprocessing.geoprocessing.get_nodata_from_uri(
        flow_accumulation_uri)

    w_accumulation_uri = os.path.join(intermediate_dir, 'w_accumulation%s.tif' % file_suffix)
    s_accumulation_uri = os.path.join(intermediate_dir, 's_accumulation%s.tif' % file_suffix)
    for factor_uri, accumulation_uri in [
            (thresholded_w_factor_uri, w_accumulation_uri),
            (thresholded_slope_uri, s_accumulation_uri)]:
        LOGGER.info("calculating %s", accumulation_uri)
        pygeoprocessing.routing.route_flux(
            flow_direction_uri, aligned_dem_uri, factor_uri,
            zero_absorption_source_uri, loss_uri, accumulation_uri, 'flux_only',
            aoi_uri=args['watersheds_uri'])

    LOGGER.info("calculating w_bar")

    w_bar_uri = os.path.join(intermediate_dir, 'w_bar%s.tif' % file_suffix)
    w_bar_nodata = pygeoprocessing.geoprocessing.get_nodata_from_uri(w_accumulation_uri)
    s_bar_uri = os.path.join(intermediate_dir, 's_bar%s.tif' % file_suffix)
    s_bar_nodata = pygeoprocessing.geoprocessing.get_nodata_from_uri(s_accumulation_uri)
    for bar_nodata, accumulation_uri, bar_uri in [
            (w_bar_nodata, w_accumulation_uri, w_bar_uri),
            (s_bar_nodata, s_accumulation_uri, s_bar_uri)]:
        LOGGER.info("calculating %s", accumulation_uri)
        def bar_op(base_accumulation, flow_accumulation):
            return numpy.where(
                (base_accumulation != bar_nodata) & (flow_accumulation != flow_accumulation_nodata),
                base_accumulation / flow_accumulation, bar_nodata)
        pygeoprocessing.geoprocessing.vectorize_datasets(
            [accumulation_uri, flow_accumulation_uri], bar_op, bar_uri,
            gdal.GDT_Float32, bar_nodata, out_pixel_size, "intersection",
            dataset_to_align_index=0, vectorize_op=False)

    LOGGER.info('calculating d_up')
    d_up_uri = os.path.join(intermediate_dir, 'd_up%s.tif' % file_suffix)
    cell_area = out_pixel_size ** 2
    d_up_nodata = -1.0
    def d_up(w_bar, s_bar, flow_accumulation):
        """Calculate the d_up index
            w_bar * s_bar * sqrt(upstream area) """
        d_up_array = w_bar * s_bar * numpy.sqrt(flow_accumulation * cell_area)
        return numpy.where(
            (w_bar != w_bar_nodata) & (s_bar != s_bar_nodata) &
            (flow_accumulation != flow_accumulation_nodata), d_up_array,
            d_up_nodata)
    pygeoprocessing.geoprocessing.vectorize_datasets(
        [w_bar_uri, s_bar_uri, flow_accumulation_uri], d_up, d_up_uri,
        gdal.GDT_Float32, d_up_nodata, out_pixel_size, "intersection",
        dataset_to_align_index=0, vectorize_op=False)

    LOGGER.info('calculate WS factor')
    ws_factor_inverse_uri = os.path.join(
        intermediate_dir, 'ws_factor_inverse%s.tif' % file_suffix)
    ws_nodata = -1.0
    slope_nodata = pygeoprocessing.geoprocessing.get_nodata_from_uri(
        thresholded_slope_uri)

    def ws_op(w_factor, s_factor):
        #calculating the inverse so we can use the distance to stream factor function
        return numpy.where(
            (w_factor != w_nodata) & (s_factor != slope_nodata),
            1.0 / (w_factor * s_factor), ws_nodata)

    pygeoprocessing.geoprocessing.vectorize_datasets(
        [thresholded_w_factor_uri, thresholded_slope_uri], ws_op, ws_factor_inverse_uri,
        gdal.GDT_Float32, ws_nodata, out_pixel_size, "intersection",
        dataset_to_align_index=0, vectorize_op=False)

    LOGGER.info('calculating d_dn')
    d_dn_uri = os.path.join(intermediate_dir, 'd_dn%s.tif' % file_suffix)
    pygeoprocessing.routing.routing_core.distance_to_stream(
        flow_direction_uri, stream_uri, d_dn_uri, factor_uri=ws_factor_inverse_uri)

    LOGGER.info('calculating downstream distance')
    downstream_distance_uri = os.path.join(
        intermediate_dir, 'downstream_distance%s.tif' % file_suffix)
    pygeoprocessing.routing.routing_core.distance_to_stream(
        flow_direction_uri, stream_uri, downstream_distance_uri)
    downstream_distance_nodata = pygeoprocessing.geoprocessing.get_nodata_from_uri(
        downstream_distance_uri)

    LOGGER.info('calculate ic')
    ic_factor_uri = os.path.join(intermediate_dir, 'ic_factor%s.tif' % file_suffix)
    ic_nodata = -9999.0
    d_up_nodata = pygeoprocessing.geoprocessing.get_nodata_from_uri(d_up_uri)
    d_dn_nodata = pygeoprocessing.geoprocessing.get_nodata_from_uri(d_dn_uri)
    def ic_op(d_up, d_dn):
        nodata_mask = (d_up == d_up_nodata) | (d_dn == d_dn_nodata) | (d_up == 0) | (d_dn == 0)
        return numpy.where(
            nodata_mask, ic_nodata, numpy.log10(d_up/d_dn))
    pygeoprocessing.geoprocessing.vectorize_datasets(
        [d_up_uri, d_dn_uri], ic_op, ic_factor_uri,
        gdal.GDT_Float32, ic_nodata, out_pixel_size, "intersection",
        dataset_to_align_index=0, vectorize_op=False)

    ic_min, ic_max, _, _ = pygeoprocessing.geoprocessing.get_statistics_from_uri(ic_factor_uri)
    ic_0_param = (ic_min + ic_max) / 2.0
    k_param = float(args['k_param'])

    lulc_mask_uri = pygeoprocessing.geoprocessing.temporary_filename()
    current_l_lulc_uri = pygeoprocessing.geoprocessing.temporary_filename()
    l_lulc_temp_uri = pygeoprocessing.geoprocessing.temporary_filename()

    for nutrient in nutrients_to_process:
        #calculate l for each lulc type
        LOGGER.info('calculating l lulc raster set')
        l_lulc_uri = os.path.join(
            intermediate_dir, 'l_lulc_%s%s.tif' % (nutrient, file_suffix))
        l_lulc_nodata = -1.0
        pygeoprocessing.geoprocessing.new_raster_from_base_uri(
            lulc_uri, l_lulc_uri, 'GTiff', l_lulc_nodata, gdal.GDT_Float32,
            fill_value=l_lulc_nodata)

        for lulc_code in get_unique_lulc_codes(lulc_uri):
            #no reason to process the nodata area, we already mask that out
            if lulc_code == lulc_nodata:
                continue
            LOGGER.info('processing lulc code %d for l_lulc' % lulc_code)
            mask_nodata = 2
            def mask_lulc_type(lulc_array):
                result = numpy.zeros(lulc_array.shape, dtype=numpy.int8)
                result[lulc_array == lulc_nodata] = mask_nodata
                result[lulc_array == lulc_code] = 1
                return result

            pygeoprocessing.geoprocessing.vectorize_datasets(
                [lulc_uri], mask_lulc_type, lulc_mask_uri, gdal.GDT_Byte,
                mask_nodata, out_pixel_size, 'intersection', vectorize_op=False)

            pygeoprocessing.routing.routing_core.distance_to_stream(
                flow_direction_uri, stream_uri, current_l_lulc_uri,
                factor_uri=lulc_mask_uri)

            current_l_lulc_nodata = pygeoprocessing.geoprocessing.get_nodata_from_uri(current_l_lulc_uri)
            def add_to_l_lulc(current_l_lulc_array, lulc_mask_array, l_lulc_array):
                result = l_lulc_array.copy()
                mask = (lulc_mask_array == 1) & (current_l_lulc_array != current_l_lulc_nodata)
                result[mask] = current_l_lulc_array[mask]
                return result

            #copy the file to avoid aliasing
            shutil.copyfile(l_lulc_uri, l_lulc_temp_uri)
            pygeoprocessing.geoprocessing.vectorize_datasets(
                [current_l_lulc_uri, lulc_mask_uri, l_lulc_temp_uri],
                add_to_l_lulc, l_lulc_uri, gdal.GDT_Float32,
                l_lulc_nodata, out_pixel_size, 'intersection',
                vectorize_op=False)

        ndr_max_uri = os.path.join(
            intermediate_dir, 'ndr_max_%s%s.tif' % (nutrient, file_suffix))
        ndr_max_nodata = -1.0

        def calculate_ndr_max(l_lulc_array, lulc_array):
            critical_length_array = numpy.empty(lulc_array.shape, dtype=numpy.float32)
            critical_length_array[:] = 1
            eff_lulc_array = numpy.empty(lulc_array.shape, dtype=numpy.float32)
            for lucode in numpy.unique(lulc_array):
                if lucode == lulc_nodata:
                    continue
                mask = lulc_array == lucode
                eff_lulc_array[mask] = lucode_to_parameters[lucode]['eff_%s' % nutrient]
                critical_length_array[mask] = lucode_to_parameters[lucode]['crit_len_%s' % nutrient]
            return numpy.where(
                (lulc_array == lulc_nodata) | (l_lulc_array == l_lulc_nodata), ndr_max_nodata,
                1 - eff_lulc_array * (1 - numpy.exp(-5.0 * l_lulc_array / critical_length_array)))

        pygeoprocessing.geoprocessing.vectorize_datasets(
            [l_lulc_uri, lulc_uri],
            calculate_ndr_max, ndr_max_uri, gdal.GDT_Float32,
            ndr_max_nodata, out_pixel_size, 'intersection',
            vectorize_op=False)

        LOGGER.info('calculate NDR')
        ndr_uri = os.path.join(
            intermediate_dir, 'ndr_%s%s.tif' % (nutrient, file_suffix))
        ndr_nodata = -1.0
        def calculate_ndr(ndr_max_array, ic_array):
            return numpy.where(
                (ndr_max_array == ndr_nodata) | (ic_array == ic_nodata), ndr_nodata,
                ndr_max_array * 1.0 / (1.0 + numpy.exp((ic_0_param - ic_array) / k_param)))

        pygeoprocessing.geoprocessing.vectorize_datasets(
            [ndr_max_uri, ic_factor_uri], calculate_ndr, ndr_uri,
            gdal.GDT_Float32, ndr_nodata, out_pixel_size, 'intersection',
            vectorize_op=False)

        LOGGER.info('calculate subsurface NDR')
        ndr_subsurface_uri = os.path.join(
            intermediate_dir, 'ndr_subsurface_%s%s.tif' % (nutrient, file_suffix))

        subsurface_eff = float(args['subsurface_eff_' + nutrient])
        crit_subsurface_len = float(args['subsurface_critical_length_' + nutrient])

        def calculate_subsurface_ndr(downstream_distance):
            return numpy.where(
                (downstream_distance == downstream_distance_nodata), ndr_nodata,
                1 - subsurface_eff * (1.0 - numpy.exp(-5 * downstream_distance / crit_subsurface_len)))
        #pygeoprocessing.geoprocessing.vectorize_datasets(
        #    [downstream_distance_uri], calculate_subsurface_ndr, ndr_subsurface_uri,
        #    gdal.GDT_Float32, ndr_nodata, out_pixel_size, 'intersection',
        #    vectorize_op=False)

        export_uri[nutrient] = os.path.join(
            output_dir, '%s_export%s.tif' % (nutrient, file_suffix))

        load_nodata = pygeoprocessing.geoprocessing.get_nodata_from_uri(load_uri[nutrient])
        export_nodata = -1.0
        def calculate_export(load_array, ndr_array):
            return numpy.where(
                (load_array == load_nodata) | (ndr_array == ndr_nodata),
                export_nodata, load_array * ndr_array)

        pygeoprocessing.geoprocessing.vectorize_datasets(
            [load_uri[nutrient], ndr_uri],
            calculate_export,
            export_uri[nutrient], gdal.GDT_Float32,
            export_nodata, out_pixel_size, "intersection", vectorize_op=False)

        #Summarize the results in terms of watershed:
        LOGGER.info("Summarizing the results of nutrient %s", nutrient)
        load_tot = pygeoprocessing.geoprocessing.aggregate_raster_values_uri(
            load_uri[nutrient], args['watersheds_uri'], 'ws_id').total
        export_tot = pygeoprocessing.geoprocessing.aggregate_raster_values_uri(
            export_uri[nutrient], args['watersheds_uri'], 'ws_id').total

        field_summaries['%s_load_tot' % nutrient] = load_tot
        field_summaries['%s_exp_tot' % nutrient] = export_tot
        field_header_order = (
            [x % nutrient for x in ['%s_load_tot', '%s_exp_tot']] + field_header_order)

    LOGGER.info('Writing summaries to output shapefile')
    add_fields_to_shapefile('ws_id', field_summaries, output_layer, field_header_order)

    LOGGER.info('cleaning up temp files')
    for uri in [
            zero_absorption_source_uri, loss_uri, lulc_mask_uri,
            current_l_lulc_uri, l_lulc_temp_uri, dem_uri, lulc_uri]:
        os.remove(uri)


def add_fields_to_shapefile(
        key_field, field_summaries, output_layer, field_header_order=None):
    """Adds fields and their values indexed by key fields to an OGR
        layer open for writing.

        key_field - name of the key field in the output_layer that
            uniquely identifies each polygon.
        field_summaries - a dictionary indexed by the desired field
            name to place in the polygon that indexes to another
            dictionary indexed by key_field value to map to that
            particular polygon.  ex {'field_name_1': {key_val1: value,
            key_val2: value}, 'field_name_2': {key_val1: value, etc.
        output_layer - an open writable OGR layer
        field_header_order - a list of field headers in the order we
            wish them to appear in the output table, if None then
            random key order in field summaries is used.

        returns nothing"""
    if field_header_order == None:
        field_header_order = field_summaries.keys()

    for field_name in field_header_order:
        field_def = ogr.FieldDefn(field_name, ogr.OFTReal)
        output_layer.CreateField(field_def)

    #Initialize each feature field to 0.0
    for feature_id in xrange(output_layer.GetFeatureCount()):
        feature = output_layer.GetFeature(feature_id)
        for field_name in field_header_order:
            try:
                ws_id = feature.GetFieldAsInteger(key_field)
                feature.SetField(
                    field_name, float(field_summaries[field_name][ws_id]))
            except KeyError:
                LOGGER.warning('unknown field %s', field_name)
                feature.SetField(field_name, 0.0)
        #Save back to datasource
        output_layer.SetFeature(feature)

def get_unique_lulc_codes(dataset_uri):
    """Find all the values in the input raster and return a list of unique
        values in that raster

        dataset_uri - uri to a land cover map that has integer values

        returns a unique list of codes in dataset_uri"""

    dataset = gdal.Open(dataset_uri)
    dataset_band = dataset.GetRasterBand(1)
    block_size = dataset_band.GetBlockSize()

    n_rows, n_cols = dataset.RasterYSize, dataset.RasterXSize
    cols_per_block, rows_per_block = block_size[0], block_size[1]
    n_col_blocks = int(math.ceil(n_cols / float(cols_per_block)))
    n_row_blocks = int(math.ceil(n_rows / float(rows_per_block)))

    unique_codes = set()
    for row_block_index in xrange(n_row_blocks):
        row_offset = row_block_index * rows_per_block
        row_block_width = n_rows - row_offset
        if row_block_width > rows_per_block:
            row_block_width = rows_per_block

        for col_block_index in xrange(n_col_blocks):
            col_offset = col_block_index * cols_per_block
            col_block_width = n_cols - col_offset
            if col_block_width > cols_per_block:
                col_block_width = cols_per_block
            result = dataset_band.ReadAsArray(
                xoff=col_offset, yoff=row_offset,
                win_xsize=col_block_width,
                win_ysize=row_block_width)
            unique_codes.update(numpy.unique(result))

    return unique_codes


def _prepare(**args):
    """A function to preprocess the static data that goes into the NDR model
        that is unlikely to change when running a batch process.

        args['dem_uri'] - dem layer
        args['watersheds_uri'] - layer to AOI/watersheds

        return a dictionary with the keys:
            'aligned_dem_uri': aligned_dem_uri,
            'thresholded_slope_uri': thresholded_slope_uri,
            'flow_accumulation_uri': flow_accumulation_uri,
            'flow_direction_uri': flow_direction_uri
    """

    intermediate_dir = os.path.join(args['workspace_dir'], 'prepared_data')

    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)

    dem_pixel_size = pygeoprocessing.geoprocessing.get_cell_size_from_uri(args['dem_uri'])

    #Align all the input rasters
    aligned_dem_uri = pygeoprocessing.geoprocessing.temporary_filename()
    pygeoprocessing.geoprocessing.align_dataset_list(
        [args['dem_uri']], [aligned_dem_uri], ['nearest'], dem_pixel_size,
        'intersection', dataset_to_align_index=0,
        aoi_uri=args['watersheds_uri'])

    #Calculate flow accumulation
    LOGGER.info("calculating flow accumulation")
    flow_accumulation_uri = os.path.join(
        intermediate_dir, 'flow_accumulation.tif')
    flow_direction_uri = os.path.join(
        intermediate_dir, 'flow_direction.tif')

    pygeoprocessing.routing.flow_direction_d_inf(
        aligned_dem_uri, flow_direction_uri)
    pygeoprocessing.routing.flow_accumulation(
        flow_direction_uri, aligned_dem_uri, flow_accumulation_uri)

    #Calculate slope
    LOGGER.info("Calculating slope")
    original_slope_uri = os.path.join(intermediate_dir, 'slope.tif')
    thresholded_slope_uri = os.path.join(
        intermediate_dir, 'thresholded_slope.tif')
    pygeoprocessing.geoprocessing.calculate_slope(aligned_dem_uri, original_slope_uri)
    slope_nodata = pygeoprocessing.geoprocessing.get_nodata_from_uri(original_slope_uri)
    def threshold_slope(slope):
        '''Threshold slope between 0.001 and 1.0'''
        slope_copy = slope.copy()
        nodata_mask = slope == slope_nodata
        slope_copy[slope < 0.001] = 0.001
        slope_copy[slope > 1.0] = 1.0
        slope_copy[nodata_mask] = slope_nodata
        return slope_copy
    pygeoprocessing.geoprocessing.vectorize_datasets(
        [original_slope_uri], threshold_slope, thresholded_slope_uri,
        gdal.GDT_Float32, slope_nodata, dem_pixel_size, "intersection",
        dataset_to_align_index=0, vectorize_op=False)

    return {
        'aligned_dem_uri': aligned_dem_uri,
        'thresholded_slope_uri': thresholded_slope_uri,
        'flow_accumulation_uri': flow_accumulation_uri,
        'flow_direction_uri': flow_direction_uri
    }
