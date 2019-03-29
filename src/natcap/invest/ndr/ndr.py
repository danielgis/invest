"""InVEST Nutrient Delivery Ratio (NDR) module."""
from __future__ import absolute_import
import itertools
import logging
import os

from osgeo import gdal
from osgeo import ogr
import numpy
import taskgraph
import pygeoprocessing
import pygeoprocessing.routing

from .. import validation
from .. import utils
from natcap.invest.ndr import ndr_core

LOGGER = logging.getLogger('natcap.invest.ndr.ndr')
logging.getLogger('taskgraph').setLevel(logging.INFO)

_OUTPUT_BASE_FILES = {
    'n_export_path': 'n_export.tif',
    'p_export_path': 'p_export.tif',
    'watershed_results_ndr_path': 'watershed_results_ndr.shp',
    }

_INTERMEDIATE_BASE_FILES = {
    'ic_factor_path': 'ic_factor.tif',
    'load_n_path': 'load_n.tif',
    'load_p_path': 'load_p.tif',
    'modified_load_n_path': 'modified_load_n.tif',
    'modified_load_p_path': 'modified_load_p.tif',
    'modified_sub_load_n_path': 'modified_sub_load_n.tif',
    'modified_sub_load_p_path': 'modified_sub_load_p.tif',
    'ndr_n_path': 'ndr_n.tif',
    'ndr_p_path': 'ndr_p.tif',
    'runoff_proxy_index_path': 'runoff_proxy_index.tif',
    's_accumulation_path': 's_accumulation.tif',
    's_bar_path': 's_bar.tif',
    's_factor_inverse_path': 's_factor_inverse.tif',
    'stream_path': 'stream.tif',
    'sub_crit_len_n_path': 'sub_crit_len_n.tif',
    'sub_crit_len_p_path': 'sub_crit_len_p.tif',
    'sub_eff_n_path': 'sub_eff_n.tif',
    'sub_eff_p_path': 'sub_eff_p.tif',
    'sub_effective_retention_n_path': 'sub_effective_retention_n.tif',
    'sub_effective_retention_p_path': 'sub_effective_retention_p.tif',
    'sub_load_n_path': 'sub_load_n.tif',
    'sub_load_p_path': 'sub_load_p.tif',
    'sub_ndr_n_path': 'sub_ndr_n.tif',
    'sub_ndr_p_path': 'sub_ndr_p.tif',
    'crit_len_n_path': 'crit_len_n.tif',
    'crit_len_p_path': 'crit_len_p.tif',
    'd_dn_path': 'd_dn.tif',
    'd_up_path': 'd_up.tif',
    'eff_n_path': 'eff_n.tif',
    'eff_p_path': 'eff_p.tif',
    'effective_retention_n_path': 'effective_retention_n.tif',
    'effective_retention_p_path': 'effective_retention_p.tif',
    'flow_accumulation_path': 'flow_accumulation.tif',
    'flow_direction_path': 'flow_direction.tif',
    'thresholded_slope_path': 'thresholded_slope.tif',
    'processed_cell_path': 'processed_cell.tif',
    }

_CACHE_BASE_FILES = {
    'filled_dem_path': 'filled_dem.tif',
    'aligned_dem_path': 'aligned_dem.tif',
    'loss_path': 'loss.tif',
    'slope_path': 'slope.tif',
    'aligned_lulc_path': 'aligned_lulc.tif',
    'aligned_runoff_proxy_path': 'aligned_runoff_proxy.tif',
    'zero_absorption_source_path': 'zero_absorption_source.tif',
    'runoff_mean_pickle_path': 'runoff_mean.pickle',
    }

_TARGET_NODATA = -1


def execute(args):
    """Nutrient Delivery Ratio.

    Parameters:
        args['workspace_dir'] (string):  path to current workspace
        args['dem_path'] (string): path to digital elevation map raster
        args['lulc_path'] (string): a path to landcover map raster
        args['runoff_proxy_path'] (string): a path to a runoff proxy raster
        args['watersheds_path'] (string): path to the watershed shapefile
        args['biophysical_table_path'] (string): path to csv table on disk
            containing nutrient retention values.

            For each nutrient type [t] in args['calc_[t]'] that is true, must
            contain the following headers:

            'load_[t]', 'eff_[t]', 'crit_len_[t]'

            If args['calc_n'] is True, must also contain the header
            'proportion_subsurface_n' field.

        args['calc_p'] (boolean): if True, phosphorous is modeled,
            additionally if True then biophysical table must have p fields in
            them
        args['calc_n'] (boolean): if True nitrogen will be modeled,
            additionally biophysical table must have n fields in them.
        args['results_suffix'] (string): (optional) a text field to append to
            all output files
        args['threshold_flow_accumulation']: a number representing the flow
            accumulation in terms of upstream pixels.
        args['n_workers'] (int): if present, indicates how many worker
            processes should be used in parallel processing. -1 indicates
            single process mode, 0 is single process but non-blocking mode,
            and >= 1 is number of processes.

    Returns:
        None

    """
    def _validate_inputs(nutrients_to_process, lucode_to_parameters):
        """Validate common errors in inputs.

        Parameters:
            nutrients_to_process (list): list of 'n' and/or 'p'
            lucode_to_parameters (dictionary): biophysical input table mapping
                lucode to dictionary of table parameters.  Used to validate
                the correct columns are input

        Returns:
            None

        Raises:
            ValueError whenever a missing field in the parameter table is
            detected along with a message describing every missing field.

        """
        # Make sure all the nutrient inputs are good
        if len(nutrients_to_process) == 0:
            raise ValueError("Neither phosphorous nor nitrogen was selected"
                             " to be processed.  Choose at least one.")

        # Build up a list that'll let us iterate through all the input tables
        # and check for the required rows, and report errors if something
        # is missing.
        row_header_table_list = []

        lu_parameter_row = lucode_to_parameters.values()[0]
        row_header_table_list.append(
            (lu_parameter_row, ['load_', 'eff_', 'crit_len_'],
             args['biophysical_table_path']))

        missing_headers = []
        for row, header_prefixes, table_type in row_header_table_list:
            for nutrient_id in nutrients_to_process:
                for header_prefix in header_prefixes:
                    header = header_prefix + nutrient_id
                    if header not in row:
                        missing_headers.append(
                            "Missing header %s from %s" % (
                                header, table_type))

        # proportion_subsurface_n is a special case in which phosphorous does
        # not have an equivalent.
        if ('n' in nutrients_to_process and
                'proportion_subsurface_n' not in lu_parameter_row):
            missing_headers.append(
                "Missing header proportion_subsurface_n from " +
                args['biophysical_table_path'])

        if len(missing_headers) > 0:
            raise ValueError('\n'.join(missing_headers))

    # Load all the tables for preprocessing
    output_dir = os.path.join(args['workspace_dir'])
    intermediate_output_dir = os.path.join(
        args['workspace_dir'], 'intermediate_outputs')
    output_dir = os.path.join(args['workspace_dir'])
    cache_dir = os.path.join(intermediate_output_dir, 'cache_dir')
    for dir_path in [output_dir, intermediate_output_dir, cache_dir]:
        try:
            os.makedirs(dir_path)
        except OSError:
            pass

    n_workers = -1  # single process mode, but adjust if in args
    if 'n_workers' in args:
        n_workers = int(args['n_workers'])
    task_graph = taskgraph.TaskGraph(
        cache_dir, n_workers, reporting_interval=5.0)

    file_suffix = utils.make_suffix_string(args, 'results_suffix')
    f_reg = utils.build_file_registry(
        [(_OUTPUT_BASE_FILES, output_dir),
         (_INTERMEDIATE_BASE_FILES, intermediate_output_dir),
         (_CACHE_BASE_FILES, cache_dir)], file_suffix)

    # Build up a list of nutrients to process based on what's checked on
    nutrients_to_process = []
    for nutrient_id in ['n', 'p']:
        if args['calc_' + nutrient_id]:
            nutrients_to_process.append(nutrient_id)

    lucode_to_parameters = utils.build_lookup_from_csv(
        args['biophysical_table_path'], 'lucode')

    _validate_inputs(nutrients_to_process, lucode_to_parameters)
    dem_info = pygeoprocessing.get_raster_info(args['dem_path'])

    base_raster_list = [
        args['dem_path'], args['lulc_path'], args['runoff_proxy_path']]
    aligned_raster_list = [
        f_reg['aligned_dem_path'], f_reg['aligned_lulc_path'],
        f_reg['aligned_runoff_proxy_path']]
    align_raster_task = task_graph.add_task(
        func=pygeoprocessing.align_and_resize_raster_stack,
        args=(
            base_raster_list, aligned_raster_list,
            ['near']*len(base_raster_list), dem_info['pixel_size'],
            'intersection'),
        kwargs={
            'base_vector_path_list': [args['watersheds_path']],
            'vector_mask_options': {
                'mask_vector_path': args['watersheds_path']}},
        target_path_list=aligned_raster_list,
        task_name='align rasters')

    fill_pits_task = task_graph.add_task(
        func=pygeoprocessing.routing.fill_pits,
        args=(
            (f_reg['aligned_dem_path'], 1), f_reg['filled_dem_path']),
        kwargs={'working_dir': cache_dir},
        dependent_task_list=[align_raster_task],
        target_path_list=[f_reg['filled_dem_path']],
        task_name='fill pits')

    flow_dir_task = task_graph.add_task(
        func=pygeoprocessing.routing.flow_dir_mfd,
        args=(
            (f_reg['filled_dem_path'], 1), f_reg['flow_direction_path']),
        kwargs={'working_dir': cache_dir},
        dependent_task_list=[fill_pits_task],
        target_path_list=[f_reg['flow_direction_path']],
        task_name='flow dir')

    flow_accum_task = task_graph.add_task(
        func=pygeoprocessing.routing.flow_accumulation_mfd,
        args=(
            (f_reg['flow_direction_path'], 1),
            f_reg['flow_accumulation_path']),
        target_path_list=[f_reg['flow_accumulation_path']],
        dependent_task_list=[flow_dir_task],
        task_name='flow accum')

    stream_extraction_task = task_graph.add_task(
        func=pygeoprocessing.routing.extract_streams_mfd,
        args=(
            (f_reg['flow_accumulation_path'], 1),
            (f_reg['flow_direction_path'], 1),
            float(args['threshold_flow_accumulation']), f_reg['stream_path']),
        target_path_list=[f_reg['stream_path']],
        dependent_task_list=[flow_accum_task],
        task_name='stream extraction')

    calculate_slope_task = task_graph.add_task(
        func=pygeoprocessing.calculate_slope,
        args=((f_reg['filled_dem_path'], 1), f_reg['slope_path']),
        target_path_list=[f_reg['slope_path']],
        dependent_task_list=[stream_extraction_task],
        task_name='calculate slope')

    threshold_slope_task = task_graph.add_task(
        func=_slope_proportion_and_threshold,
        args=(f_reg['slope_path'], f_reg['thresholded_slope_path']),
        target_path_list=[f_reg['thresholded_slope_path']],
        dependent_task_list=[calculate_slope_task],
        task_name='threshold slope')

    runoff_proxy_index_task = task_graph.add_task(
        func=normalize_raster,
        args=((f_reg['aligned_runoff_proxy_path'], 1),
              f_reg['runoff_proxy_index_path']),
        target_path_list=[f_reg['runoff_proxy_index_path']],
        dependent_task_list=[align_raster_task],
        task_name='runoff proxy mean')

    for nutrient in nutrients_to_process:
        load_path = f_reg['load_%s_path' % nutrient]
        modified_load_path = f_reg['modified_load_%s_path' % nutrient]
        # Perrine says that 'n' is the only case where we could consider a
        # prop subsurface component.  So there's a special case for that.
        if nutrient == 'n':
            subsurface_proportion_type = 'proportion_subsurface_n'
        else:
            subsurface_proportion_type = None
        LOGGER.info("Mapping %s load to LULC", nutrient)
        load_task = task_graph.add_task(
            func=calculate_load,
            args=(
                f_reg['aligned_lulc_path'], lucode_to_parameters, nutrient,
                subsurface_proportion_type, load_path),
            dependent_task_list=[align_raster_task],
            target_path_list=[load_path],
            task_name='%s load calc' % nutrient)

        modified_load_task = task_graph.add_task(
            func=multiply_rasters,
            args=([load_path, f_reg['runoff_proxy_index_path']],
                  _TARGET_NODATA, modified_load_path),
            target_path_list=[modified_load_path],
            dependent_task_list=[load_task, runoff_proxy_index_task],
            task_name='modified load %s' % nutrient)

        sub_load_path = f_reg['sub_load_%s_path' % nutrient]
        subsurface_load_task = task_graph.add_task(
            func=map_subsurface_load,
            args=(f_reg['aligned_lulc_path'], lucode_to_parameters,
                  'load_%s' % nutrient, subsurface_proportion_type,
                  sub_load_path),
            target_path_list=[sub_load_path],
            dependent_task_list=[align_raster_task],
            task_name='map subsurface load %s' % nutrient)

        modified_sub_load_path = f_reg['modified_sub_load_%s_path' % nutrient]
        modified_subsurface_load_task = task_graph.add_task(
            func=multiply_rasters,
            args=(
                [sub_load_path, f_reg['runoff_proxy_index_path']],
                _TARGET_NODATA, modified_sub_load_path),
            target_path_list=[modified_sub_load_path],
            dependent_task_list=[subsurface_load_task],
            task_name='modified subsurface load %s' % nutrient)

        sub_eff_val = args['subsurface_eff_%s' % nutrient]
        sub_crit_len = args['subsurface_critical_length_%s' % nutrient]

        eff_path = f_reg['eff_%s_path' % nutrient]
        eff_task = task_graph.add_task(
            func=calculate_ret_eff,
            args=(
                f_reg['aligned_lulc_path'], f_reg['stream_path'],
                lucode_to_parameters, 'eff_%s' % nutrient, eff_path),
            target_path_list=[eff_path],
            dependent_task_list=[align_raster_task, stream_extraction_task],
            task_name='ret eff %s' % nutrient)

        crit_len_path = f_reg['crit_len_%s_path' % nutrient]
        crit_len_task = task_graph.add_task(
            func=calculate_ret_eff,
            args=(
                f_reg['aligned_lulc_path'], f_reg['stream_path'],
                lucode_to_parameters, 'crit_len_%s' % nutrient, crit_len_path),
            target_path_list=[eff_path],
            dependent_task_list=[align_raster_task, stream_extraction_task],
            task_name='ret eff %s' % nutrient)

    s_task = task_graph.add_task(
        func=pygeoprocessing.routing.flow_accumulation_mfd,
        args=((f_reg['flow_direction_path'], 1), f_reg['s_accumulation_path']),
        kwargs={
            'weight_raster_path_band': (f_reg['thresholded_slope_path'], 1)},
        target_path_list=[f_reg['s_accumulation_path']],
        dependent_task_list=[flow_dir_task],
        task_name='route s')

    task_graph.close()
    task_graph.join()
    return

    watershed_output_datasource_uri = os.path.join(
        output_dir, 'watershed_results_ndr%s.shp' % file_suffix)
    # If there is already an existing shapefile with the same name and path,
    # delete it then copy the input shapefile into the designated output folder
    if os.path.isfile(watershed_output_datasource_uri):
        os.remove(watershed_output_datasource_uri)
    original_datasource = gdal.OpenEx(args['watersheds_path'], gdal.OF_VECTOR)
    driver = gdal.GetDriverByName('ESRI Shapefile')
    output_datasource = driver.CreateCopy(
        watershed_output_datasource_uri, original_datasource)
    output_layer = output_datasource.GetLayer()

    # need this for low level route_flux function
    natcap.invest.pygeoprocessing_0_3_3.make_constant_raster_from_base_uri(
        f_reg['aligned_dem_path'], 0.0, f_reg['zero_absorption_source_path'])

    flow_accumulation_nodata = (
        natcap.invest.pygeoprocessing_0_3_3.get_nodata_from_uri(f_reg['flow_accumulation_path']))

    LOGGER.info("calculating %s", f_reg['s_accumulation_path'])
    natcap.invest.pygeoprocessing_0_3_3.routing.route_flux(
        f_reg['flow_direction_path'], f_reg['aligned_dem_path'],
        f_reg['thresholded_slope_path'], f_reg['zero_absorption_source_path'],
        f_reg['loss_path'], f_reg['s_accumulation_path'], 'flux_only',
        aoi_uri=args['watersheds_path'])

    s_bar_nodata = natcap.invest.pygeoprocessing_0_3_3.get_nodata_from_uri(
        f_reg['s_accumulation_path'])
    LOGGER.info("calculating %s", f_reg['s_bar_path'])

    def bar_op(base_accumulation, flow_accumulation):
        """Calculate bar operation."""
        valid_mask = (
            (base_accumulation != s_bar_nodata) &
            (flow_accumulation != flow_accumulation_nodata))
        result = numpy.empty(valid_mask.shape, dtype=numpy.float32)
        result[:] = s_bar_nodata
        result[valid_mask] = (
            base_accumulation[valid_mask] / flow_accumulation[valid_mask])
        return result

    natcap.invest.pygeoprocessing_0_3_3.vectorize_datasets(
        [f_reg['s_accumulation_path'], f_reg['flow_accumulation_path']],
        bar_op, f_reg['s_bar_path'], gdal.GDT_Float32, s_bar_nodata,
        out_pixel_size, "intersection", dataset_to_align_index=0,
        vectorize_op=False)

    LOGGER.info('calculating d_up')
    cell_area = out_pixel_size ** 2
    d_up_nodata = -1.0

    def d_up(s_bar, flow_accumulation):
        """Calculate d_up index."""
        valid_mask = (
            (s_bar != s_bar_nodata) &
            (flow_accumulation != flow_accumulation_nodata))
        result = numpy.empty(valid_mask.shape, dtype=numpy.float32)
        result[:] = d_up_nodata
        result[valid_mask] = (
            s_bar[valid_mask] * numpy.sqrt(
                flow_accumulation[valid_mask] * cell_area))
        return result

    natcap.invest.pygeoprocessing_0_3_3.vectorize_datasets(
        [f_reg['s_bar_path'], f_reg['flow_accumulation_path']], d_up,
        f_reg['d_up_path'], gdal.GDT_Float32, d_up_nodata, out_pixel_size,
        "intersection", dataset_to_align_index=0, vectorize_op=False)

    LOGGER.info('calculate inverse S factor')
    s_nodata = -1.0
    slope_nodata = natcap.invest.pygeoprocessing_0_3_3.get_nodata_from_uri(
        f_reg['thresholded_slope_path'])

    def s_inverse_op(s_factor):
        """Calculate inverse of S factor."""
        valid_mask = s_factor != slope_nodata
        result = numpy.empty(valid_mask.shape, dtype=numpy.float32)
        result[valid_mask] = 1.0 / s_factor[valid_mask]
        return result

    natcap.invest.pygeoprocessing_0_3_3.vectorize_datasets(
        [f_reg['thresholded_slope_path']], s_inverse_op,
        f_reg['s_factor_inverse_path'], gdal.GDT_Float32, s_nodata,
        out_pixel_size, "intersection", dataset_to_align_index=0,
        vectorize_op=False)

    LOGGER.info('calculating d_dn')
    natcap.invest.pygeoprocessing_0_3_3.routing.distance_to_stream(
        f_reg['flow_direction_path'], f_reg['stream_path'], f_reg['d_dn_path'],
        factor_uri=f_reg['s_factor_inverse_path'])

    LOGGER.info('calculate ic')
    ic_nodata = -9999.0
    d_up_nodata = natcap.invest.pygeoprocessing_0_3_3.get_nodata_from_uri(f_reg['d_up_path'])
    d_dn_nodata = natcap.invest.pygeoprocessing_0_3_3.get_nodata_from_uri(f_reg['d_dn_path'])

    def ic_op(d_up, d_dn):
        """Calculate IC0."""
        valid_mask = (
            (d_up != d_up_nodata) & (d_dn != d_dn_nodata) & (d_up != 0) &
            (d_dn != 0))
        result = numpy.empty(valid_mask.shape, dtype=numpy.float32)
        result[:] = ic_nodata
        result[valid_mask] = numpy.log10(d_up[valid_mask] / d_dn[valid_mask])
        return result

    natcap.invest.pygeoprocessing_0_3_3.vectorize_datasets(
        [f_reg['d_up_path'], f_reg['d_dn_path']], ic_op,
        f_reg['ic_factor_path'], gdal.GDT_Float32, ic_nodata, out_pixel_size,
        "intersection", dataset_to_align_index=0, vectorize_op=False)

    ic_min, ic_max, _, _ = (
        natcap.invest.pygeoprocessing_0_3_3.get_statistics_from_uri(f_reg['ic_factor_path']))
    ic_0_param = (ic_min + ic_max) / 2.0
    k_param = float(args['k_param'])

    # define some variables outside the loop for closure
    effective_retention_nodata = None
    ndr_nodata = None
    sub_effective_retention_nodata = None
    load_nodata = None
    export_nodata = None
    field_header_order = []
    field_summaries = {}
    for nutrient in nutrients_to_process:
        effective_retention_path = (
            f_reg['effective_retention_%s_path' % nutrient])
        LOGGER.info('calculate effective retention')
        eff_path = f_reg['eff_%s_path' % nutrient]
        crit_len_path = f_reg['crit_len_%s_path' % nutrient]
        ndr_core.ndr_eff_calculation(
            f_reg['flow_direction_path'], f_reg['stream_path'], eff_path,
            crit_len_path, effective_retention_path)
        effective_retention_nodata = (
            natcap.invest.pygeoprocessing_0_3_3.get_nodata_from_uri(
                effective_retention_path))
        LOGGER.info('calculate NDR')
        ndr_path = f_reg['ndr_%s_path' % nutrient]
        ndr_nodata = -1.0

        def calculate_ndr(effective_retention_array, ic_array):
            """Calculate NDR."""
            valid_mask = (
                (effective_retention_array != effective_retention_nodata) &
                (ic_array != ic_nodata))
            result = numpy.empty(valid_mask.shape, dtype=numpy.float32)
            result[:] = ndr_nodata
            result[valid_mask] = (
                (1.0 - effective_retention_array[valid_mask]) /
                (1.0 + numpy.exp(
                    (ic_0_param - ic_array[valid_mask]) / k_param)))
            return result

        natcap.invest.pygeoprocessing_0_3_3.vectorize_datasets(
            [effective_retention_path, f_reg['ic_factor_path']],
            calculate_ndr, ndr_path, gdal.GDT_Float32, ndr_nodata,
            out_pixel_size, 'intersection', vectorize_op=False)

        sub_effective_retention_path = (
            f_reg['sub_effective_retention_%s_path' % nutrient])
        LOGGER.info('calculate subsurface effective retention')
        sub_eff_path = f_reg['sub_eff_%s_path' % nutrient]
        ndr_core.ndr_eff_calculation(
            f_reg['flow_direction_path'], f_reg['stream_path'], sub_eff_path,
            sub_crit_len_path, sub_effective_retention_path)
        sub_effective_retention_nodata = (
            natcap.invest.pygeoprocessing_0_3_3.get_nodata_from_uri(
                sub_effective_retention_path))
        LOGGER.info('calculate sub NDR')
        sub_ndr_path = f_reg['sub_ndr_%s_path' % nutrient]
        ndr_nodata = -1.0

        def calculate_sub_ndr(sub_eff_ret_array):
            """Calculate subsurface NDR."""
            valid_mask = (
                sub_eff_ret_array !=
                sub_effective_retention_nodata)
            result = numpy.empty(valid_mask.shape, dtype=numpy.float32)
            result[:] = ndr_nodata
            result[valid_mask] = (1.0 - sub_eff_ret_array[valid_mask])
            return result

        natcap.invest.pygeoprocessing_0_3_3.vectorize_datasets(
            [sub_effective_retention_path], calculate_sub_ndr, sub_ndr_path,
            gdal.GDT_Float32, ndr_nodata, out_pixel_size, 'intersection',
            vectorize_op=False)

        export_path = f_reg['%s_export_path' % nutrient]

        load_nodata = natcap.invest.pygeoprocessing_0_3_3.get_nodata_from_uri(
            load_path)
        export_nodata = -1.0

        def calculate_export(
                modified_load_array, ndr_array, modified_sub_load_array,
                sub_ndr_array):
            """Combine NDR and subsurface NDR."""
            valid_mask = (
                (modified_load_array != load_nodata) &
                (ndr_array != ndr_nodata) &
                (modified_sub_load_array != load_nodata) &
                (sub_ndr_array != ndr_nodata))
            result = numpy.empty(valid_mask.shape, dtype=numpy.float32)
            result[:] = export_nodata
            result[valid_mask] = (
                modified_load_array[valid_mask] * ndr_array[valid_mask] +
                modified_sub_load_array[valid_mask] *
                sub_ndr_array[valid_mask])
            return result

        modified_load_path = f_reg['modified_load_%s_path' % nutrient]
        modified_sub_load_path = f_reg['modified_sub_load_%s_path' % nutrient]
        natcap.invest.pygeoprocessing_0_3_3.vectorize_datasets(
            [modified_load_path, ndr_path,
             modified_sub_load_path, sub_ndr_path], calculate_export,
            export_path, gdal.GDT_Float32, export_nodata,
            out_pixel_size, "intersection", vectorize_op=False)

        # summarize the results in terms of watershed:
        LOGGER.info("Summarizing the results of nutrient %s", nutrient)
        load_path = f_reg['load_%s_path' % nutrient]
        load_tot = natcap.invest.pygeoprocessing_0_3_3.aggregate_raster_values_uri(
            load_path, args['watersheds_path'], 'ws_id').total
        export_tot = natcap.invest.pygeoprocessing_0_3_3.aggregate_raster_values_uri(
            export_path, args['watersheds_path'], 'ws_id').total

        field_summaries['%s_load_tot' % nutrient] = load_tot
        field_summaries['%s_exp_tot' % nutrient] = export_tot
        field_header_order = (
            [x % nutrient for x in ['%s_load_tot', '%s_exp_tot']] +
            field_header_order)

    LOGGER.info('Writing summaries to output shapefile')
    _add_fields_to_shapefile(
        'ws_id', field_summaries, output_layer, field_header_order)

    LOGGER.info(r'NDR complete!')
    LOGGER.info(r'  _   _    ____    ____     ')
    LOGGER.info(r' | \ |"|  |  _"\U |  _"\ u  ')
    LOGGER.info(r'<|  \| |>/| | | |\| |_) |/  ')
    LOGGER.info(r'U| |\  |uU| |_| |\|  _ <    ')
    LOGGER.info(r' |_| \_|  |____/ u|_| \_\   ')
    LOGGER.info(r' ||   \\,-.|||_   //   \\_  ')
    LOGGER.info(r' (_")  (_/(__)_) (__)  (__) ')


def _slope_proportion_and_threshold(slope_path, target_threshold_slope_path):
    """Rescale slope to proportion and threshold to between 0.005 and 1.0.

    Parameters:
        slope_path (string): a raster with slope values in percent.
        target_threshold_slope_path (string): generated raster with slope
            values as a proportion (100% is 1.0) and thresholded to values
            between 0.005 and 1.0.

    Returns:
        None.

    """
    slope_nodata = pygeoprocessing.get_raster_info(slope_path)['nodata'][0]

    def _slope_proportion_and_threshold_op(slope):
        """Rescale and threshold slope between 0.005 and 1.0."""
        valid_mask = slope != slope_nodata
        result = numpy.empty(valid_mask.shape, dtype=numpy.float32)
        result[:] = slope_nodata
        slope_fraction = slope[valid_mask] / 100
        slope_fraction[slope_fraction < 0.005] = 0.005
        slope_fraction[slope_fraction > 1.0] = 1.0
        result[valid_mask] = slope_fraction
        return result

    pygeoprocessing.raster_calculator(
        [(slope_path, 1)], _slope_proportion_and_threshold_op,
        target_threshold_slope_path, gdal.GDT_Float32, slope_nodata)


def _add_fields_to_shapefile(
        key_field, field_summaries, output_layer, field_header_order):
    """Add fields and values to an OGR layer open for writing.

    Parameters:
        key_field (string): name of the key field in the output_layer that
            uniquely identifies each polygon.
        field_summaries (dict): index for the desired field name to place in
            the polygon that indexes to another dictionary indexed by
            key_field value to map to that particular polygon.
            ex {'field_name_1': {key_val1: value,
            key_val2: value}, 'field_name_2': {key_val1: value, etc.
        output_layer (ogr.Layer): an open writable OGR layer
        field_header_order (list of string): a list of field headers in the
            order to appear in the output table.

    Returns:
        None.
    """
    for field_name in field_header_order:
        field_def = ogr.FieldDefn(field_name, ogr.OFTReal)
        field_def.SetWidth(24)
        field_def.SetPrecision(11)
        output_layer.CreateField(field_def)

    # Initialize each feature field to 0.0
    for feature_id in xrange(output_layer.GetFeatureCount()):
        feature = output_layer.GetFeature(feature_id)
        for field_name in field_header_order:
            ws_id = feature.GetFieldAsInteger(key_field)
            feature.SetField(
                field_name, float(field_summaries[field_name][ws_id]))
        # Save back to datasource
        output_layer.SetFeature(feature)


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

    required_keys = [
        'workspace_dir',
        'dem_path',
        'lulc_path',
        'runoff_proxy_path',
        'watersheds_path',
        'biophysical_table_path',
        'calc_p',
        'calc_n',
        'threshold_flow_accumulation',
        'k_param']

    if 'calc_n' in args and args['calc_n']:
        required_keys.extend([
            'subsurface_critical_length_n', 'subsurface_eff_n'])

    if 'calc_p' in args and args['calc_p']:
        required_keys.extend([
            'subsurface_critical_length_p', 'subsurface_eff_p'])

    for key in required_keys:
        if limit_to is None or limit_to == key:
            if key not in args:
                missing_key_list.append(key)
            elif args[key] in ['', None]:
                no_value_list.append(key)

    if len(missing_key_list) > 0:
        # if there are missing keys, we have raise KeyError to stop hard
        raise KeyError(
            "The following keys were expected in `args` but were missing " +
            ', '.join(missing_key_list))

    if limit_to is None and (not args['calc_p'] and not args['calc_n']):
        validation_error_list.append(
            (['calc_p', 'calc_n', 'dem_path', 'lulc_path'],
             "At least nitrogen or phosphorous must be selected"))

    if len(no_value_list) > 0:
        validation_error_list.append(
            (no_value_list, 'parameter has no value'))

    file_type_list = [
        ('lulc_path', 'raster'),
        ('dem_path', 'raster'),
        ('runoff_proxy_path', 'raster'),
        ('biophysical_table_path', 'table'),
        ('watersheds_path', 'vector')]

    # check that existing/optional files are the correct types
    with utils.capture_gdal_logging():
        for key, key_type in file_type_list:
            if (limit_to is None or limit_to == key) and key in args:
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


def normalize_raster(base_raster_path_band, target_normalized_raster_path):
    """Calculate normalize raster by dividing by the mean value.

    Parameters:
        base_raster_path_band (tuple): raster path/band tuple to calculate
            mean.
        target_normalized_raster_path (string): path to target normalized
            raster from base_raster_path_band.

    Returns:
        None.

    """
    value_sum = 0.0
    value_count = 0.0
    base_nodata = pygeoprocessing.get_raster_info(
        base_raster_path_band[0])['nodata'][base_raster_path_band[1]-1]
    for _, raster_block in pygeoprocessing.iterblocks(
            base_raster_path_band):
        valid_block = raster_block[~numpy.isclose(raster_block, base_nodata)]
        value_sum = numpy.sum(valid_block)
        value_count += valid_block.size

    value_mean = value_sum
    if value_count > 0.0:
        value_mean /= value_count

    def normalize_raster_op(array):
        """Divide values by mean."""
        valid_mask = ~numpy.isclose(array, base_nodata)
        result = numpy.empty(valid_mask.shape, dtype=numpy.float32)
        result[:] = base_nodata
        result[valid_mask] = array[valid_mask]
        if value_mean != 0:
            result[valid_mask] /= value_mean
        return result

    pygeoprocessing.raster_calculator(
        [base_raster_path_band], normalize_raster_op,
        target_normalized_raster_path, gdal.GDT_Float32, base_nodata)


def map_load_function(
        lucode_to_parameters, load_type, subsurface_proportion_type=None):
    """Function generator to map arbitrary nutrient type to surface load.

    Parameters:
        lucode_to_parameters (dict): maps an integer code to a dictionary
            that contains at least the key `load_type` which maps to a
            float.
        load_type (string): either 'n' or 'p', used for indexing headers
        subsurface_proportion_type (string): if None no subsurface transfer
            is mapped.  Otherwise indexed from lucode_to_parameters

    Returns:
        map_load (function(lucode_array)): a function that can be passed to
            vectorize_datasets.
    """
    def map_load(lucode_array):
        """Convert unit load to total load & handle nodata."""
        result = numpy.empty(lucode_array.shape)
        result[:] = nodata_load
        for lucode in numpy.unique(lucode_array):
            if lucode != nodata_landuse:
                if subsurface_proportion_type is not None:
                    if lucode not in lucode_to_parameters:
                        raise KeyError(
                            "The LULC code %s can\'t be found in the"
                            " biophysical table." % lucode)
                    result[lucode_array == lucode] = (
                        lucode_to_parameters[lucode][load_type] *
                        (1 - lucode_to_parameters[lucode]
                         [subsurface_proportion_type]) *
                        cell_area_ha)
                else:
                    result[lucode_array == lucode] = (
                        lucode_to_parameters[lucode][load_type] *
                        cell_area_ha)
        return result
    return map_load


def calculate_load(
        lulc_raster_path, lucode_to_parameters, nutrient_id,
        subsurface_proportion_type, target_load_raster):
    """Calculate load raster.

    Returns:
        None.

    """
    lulc_raster_info = pygeoprocessing.get_raster_info(lulc_raster_path)
    nodata_landuse = lulc_raster_info['nodata'][0]
    cell_area_ha = abs(numpy.prod(lulc_raster_info['pixel_size'])) * 0.0001

    def map_load_function(load_type, subsurface_proportion_type=None):
        """Function generator to map arbitrary nutrient type to surface load.

        Parameters:
            load_type (string): either 'n' or 'p', used for indexing headers
            subsurface_proportion_type (string): if None no subsurface transfer
                is mapped.  Otherwise indexed from lucode_to_parameters

        Returns:
            map_load (function(lucode_array)): a function that can be passed to
                vectorize_datasets.
        """
        def map_load(lucode_array):
            """Convert unit load to total load & handle nodata."""
            result = numpy.empty(lucode_array.shape)
            result[:] = _TARGET_NODATA
            for lucode in numpy.unique(lucode_array):
                if lucode != nodata_landuse:
                    if subsurface_proportion_type is not None:
                        if lucode not in lucode_to_parameters:
                            raise KeyError(
                                "The LULC code %s can\'t be found in the"
                                " biophysical table." % lucode)
                        result[lucode_array == lucode] = (
                            lucode_to_parameters[lucode][load_type] *
                            (1 - lucode_to_parameters[lucode]
                             [subsurface_proportion_type]) *
                            cell_area_ha)
                    else:
                        result[lucode_array == lucode] = (
                            lucode_to_parameters[lucode][load_type] *
                            cell_area_ha)
            return result
        return map_load

    pygeoprocessing.raster_calculator(
        [(lulc_raster_path, 1)], map_load_function(
            'load_%s' % nutrient_id, subsurface_proportion_type),
        target_load_raster, gdal.GDT_Float32, _TARGET_NODATA)


def multiply_rasters(raster_path_list, target_nodata, target_result_path):
    """Multiply the rasters in `raster_path_list`.

    Parameters:
        raster_path_list (list): list of single band raster paths.
        target_nodata (float): desired target nodata value.
        target_result_path (string): path to float 32 target raster
            multiplied where all rasters are not nodata.

    Returns:
        None.

    """
    def mult_op(*array_nodata_list):
        """Multiply non-nodata stacks."""
        result = numpy.empty(array_nodata_list[0].shape)
        result[:] = target_nodata
        valid_mask = ~numpy.isclose(
            array_nodata_list[0], array_nodata_list[1])
        for array, nodata in zip(*[iter(array_nodata_list[2:])]*2):
            valid_mask &= ~numpy.isclose(array, nodata)
        result[valid_mask] = array_nodata_list[0][valid_mask]
        for array in array_nodata_list[2::2]:
            result[valid_mask] *= array[valid_mask]
        return result

    # make a list of (raster_path_band, nodata) tuples, then flatten it
    path_nodata_list = list(itertools.chain(*[
        ((path, 1),
         (pygeoprocessing.get_raster_info(path)['nodata'][0], 'raw'))
        for path in raster_path_list]))
    pygeoprocessing.raster_calculator(
        path_nodata_list, mult_op, target_result_path,
        gdal.GDT_Float32, target_nodata)


def map_subsurface_load(
        lulc_raster_path, lucode_to_parameters, load_id,
        subsurface_proportion_type, target_sub_load_path):
    """Calculate subsurface load from landcover raster.

    Parameters:
        lulc_raster_path (string): path to landcover raster.
        lucode_to_parameters (dict): maps landcover codes to a dictionary that
            can be indexed by `load_id` and then by
            `subsurface_proportion_type`.
        load_id (string): either 'load_n' or 'load_p', used to look up in
            `lucode_to_parameters`.
        subsurface_proportion_type (string): if None no subsurface transfer
            is mapped.  Otherwise indexed from lucode_to_parameters.
        target_sub_load_path (string): path to target raster.

    Returns:
        None.

    """
    lulc_raster_info = pygeoprocessing.get_raster_info(lulc_raster_path)
    nodata_landuse = lulc_raster_info['nodata'][0]
    cell_area_ha = abs(numpy.prod(lulc_raster_info['pixel_size'])) * 0.0001

    keys = sorted(numpy.array(lucode_to_parameters.keys()))
    surface_values = numpy.array(
        [lucode_to_parameters[x][load_id] for x in keys])
    if subsurface_proportion_type is not None:
        subsurface_values = numpy.array(
            [lucode_to_parameters[x][subsurface_proportion_type]
             for x in keys])

    def map_subsurface_load_op(lucode_array):
        """Convert unit load to total load & handle nodata."""
        # If we don't have subsurface, just return 0.0.
        if subsurface_proportion_type is None:
            return numpy.where(
                lucode_array != nodata_landuse, 0, _TARGET_NODATA)

        valid_mask = lucode_array != nodata_landuse
        result = numpy.empty(valid_mask.shape, dtype=numpy.float32)
        result[:] = _TARGET_NODATA
        index = numpy.digitize(
            lucode_array[valid_mask].ravel(), keys, right=True)
        result[valid_mask] = (
            surface_values[index] * subsurface_values[index] *
            cell_area_ha)
        return result

    pygeoprocessing.raster_calculator(
        [(lulc_raster_path, 1)], map_subsurface_load_op, target_sub_load_path,
        gdal.GDT_Float32, _TARGET_NODATA)


def calculate_ret_eff(
        lulc_raster_path, stream_path, lucode_to_parameters, eff_id,
        target_eff_path):
    """Make retention efficiency raster from landcover.

    Parameters:
        lulc_raster_path (string): path to landcover raster.
        stream_path (string) path to stream layer 0, no stream 1 stream.
        lucode_to_parameters (dict) mapping of landcover code to a dictionary
            that contains the key in `eff_id`
        eff_id (string): the id in the lookup table with values to map
            landcover to efficiency.
        target_eff_path (string): target raster that contains the mapping of
            landcover codes to retention efficiency values except where there
            is a stream in which case the retention efficiency is 0.

    Returns:
        None.

    """
    keys = sorted(numpy.array(lucode_to_parameters.keys()))
    values = numpy.array(
        [lucode_to_parameters[x][eff_id] for x in keys])

    nodata_landuse = pygeoprocessing.get_raster_info(
        lulc_raster_path)['nodata'][0]
    nodata_stream = pygeoprocessing.get_raster_info(stream_path)['nodata'][0]

    def map_eff_op(lucode_array, stream_array):
        """Map efficiency from LULC and handle nodata/streams."""
        valid_mask = (
            (lucode_array != nodata_landuse) &
            (stream_array != nodata_stream))
        result = numpy.empty(valid_mask.shape, dtype=numpy.float32)
        result[:] = _TARGET_NODATA
        index = numpy.digitize(
            lucode_array[valid_mask].ravel(), keys, right=True)
        result[valid_mask] = (
            values[index] * (1 - stream_array[valid_mask]))
        return result

    pygeoprocessing.raster_calculator(
        ((lulc_raster_path, 1), (stream_path, 1)), map_eff_op,
        target_eff_path, gdal.GDT_Float32, _TARGET_NODATA)
