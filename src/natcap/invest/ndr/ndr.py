"""InVEST Nutrient Delivery Ratio (NDR) module."""
from __future__ import absolute_import
import pickle
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
from . import ndr_core

LOGGER = logging.getLogger('natcap.invest.ndr.ndr')

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
    'surface_load_n_path': 'surface_load_n.tif',
    'surface_load_p_path': 'surface_load_p.tif',
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
    'dist_to_channel_path': 'dist_to_channel.tif',
    }

_CACHE_BASE_FILES = {
    'filled_dem_path': 'filled_dem.tif',
    'aligned_dem_path': 'aligned_dem.tif',
    'slope_path': 'slope.tif',
    'aligned_lulc_path': 'aligned_lulc.tif',
    'aligned_runoff_proxy_path': 'aligned_runoff_proxy.tif',
    'runoff_mean_pickle_path': 'runoff_mean.pickle',
    'surface_load_n_pickle_path': 'surface_load_n.pickle',
    'surface_load_p_pickle_path': 'surface_load_p.pickle',
    'subsurface_load_n_pickle_path': 'subsurface_load_n.pickle',
    'subsurface_load_p_pickle_path': 'subsurface_load_p.pickle',
    'export_n_pickle_path': 'export_n.pickle',
    'export_p_pickle_path': 'export_p.pickle',
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

    # these are used for aggregation in the last step
    field_pickle_map = {}
    field_header_order_list = []

    create_vector_task = task_graph.add_task(
        func=create_vector_copy,
        args=(args['watersheds_path'], f_reg['watershed_results_ndr_path']),
        target_path_list=[f_reg['watershed_results_ndr_path']],
        task_name='create target vector')

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
        func=_normalize_raster,
        args=((f_reg['aligned_runoff_proxy_path'], 1),
              f_reg['runoff_proxy_index_path']),
        target_path_list=[f_reg['runoff_proxy_index_path']],
        dependent_task_list=[align_raster_task],
        task_name='runoff proxy mean')

    s_task = task_graph.add_task(
        func=pygeoprocessing.routing.flow_accumulation_mfd,
        args=((f_reg['flow_direction_path'], 1), f_reg['s_accumulation_path']),
        kwargs={
            'weight_raster_path_band': (f_reg['thresholded_slope_path'], 1)},
        target_path_list=[f_reg['s_accumulation_path']],
        dependent_task_list=[flow_dir_task, threshold_slope_task],
        task_name='route s')

    s_bar_task = task_graph.add_task(
        func=s_bar_calculate,
        args=(f_reg['s_accumulation_path'], f_reg['flow_accumulation_path'],
              f_reg['s_bar_path']),
        target_path_list=[f_reg['s_bar_path']],
        dependent_task_list=[s_task, flow_accum_task],
        task_name='calculate s bar')

    d_up_task = task_graph.add_task(
        func=d_up_calculation,
        args=(f_reg['s_bar_path'], f_reg['flow_accumulation_path'],
              f_reg['d_up_path']),
        target_path_list=[f_reg['d_up_path']],
        dependent_task_list=[s_bar_task, flow_accum_task],
        task_name='d up')

    s_inv_task = task_graph.add_task(
        func=invert_raster_values,
        args=(f_reg['thresholded_slope_path'], f_reg['s_factor_inverse_path']),
        target_path_list=[f_reg['s_factor_inverse_path']],
        dependent_task_list=[threshold_slope_task],
        task_name='s inv')

    d_dn_task = task_graph.add_task(
        func=pygeoprocessing.routing.distance_to_channel_mfd,
        args=(
            (f_reg['flow_direction_path'], 1), (f_reg['stream_path'], 1),
            f_reg['d_dn_path']),
        kwargs={'weight_raster_path_band': (
            f_reg['s_factor_inverse_path'], 1)},
        dependent_task_list=[stream_extraction_task, s_inv_task],
        target_path_list=[f_reg['d_dn_path']],
        task_name='d dn')

    dist_to_channel_task = task_graph.add_task(
        func=pygeoprocessing.routing.distance_to_channel_mfd,
        args=(
            (f_reg['flow_direction_path'], 1), (f_reg['stream_path'], 1),
            f_reg['dist_to_channel_path']),
        dependent_task_list=[stream_extraction_task],
        target_path_list=[f_reg['dist_to_channel_path']],
        task_name='dist to channel')

    ic_task = task_graph.add_task(
        func=calculate_ic,
        args=(
            f_reg['d_up_path'], f_reg['d_dn_path'], f_reg['ic_factor_path']),
        target_path_list=[f_reg['ic_factor_path']],
        dependent_task_list=[d_dn_task, d_up_task],
        task_name='calc ic')

    for nutrient in nutrients_to_process:
        load_path = f_reg['load_%s_path' % nutrient]
        modified_load_path = f_reg['modified_load_%s_path' % nutrient]
        # Perrine says that 'n' is the only case where we could consider a
        # prop subsurface component.  So there's a special case for that.
        if nutrient == 'n':
            subsurface_proportion_type = 'proportion_subsurface_n'
        else:
            subsurface_proportion_type = None
        load_task = task_graph.add_task(
            func=_calculate_load,
            args=(
                f_reg['aligned_lulc_path'], lucode_to_parameters,
                'load_%s' % nutrient, load_path),
            dependent_task_list=[align_raster_task],
            target_path_list=[load_path],
            task_name='%s load' % nutrient)

        modified_load_task = task_graph.add_task(
            func=_multiply_rasters,
            args=([load_path, f_reg['runoff_proxy_index_path']],
                  _TARGET_NODATA, modified_load_path),
            target_path_list=[modified_load_path],
            dependent_task_list=[load_task, runoff_proxy_index_task],
            task_name='modified load %s' % nutrient)

        surface_load_path = f_reg['surface_load_%s_path' % nutrient]
        surface_load_task = task_graph.add_task(
            func=_map_surface_load,
            args=(modified_load_path, f_reg['aligned_lulc_path'],
                  lucode_to_parameters, subsurface_proportion_type,
                  surface_load_path),
            target_path_list=[surface_load_path],
            dependent_task_list=[modified_load_task, align_raster_task],
            task_name='map surface load %s' % nutrient)

        subsurface_load_path = f_reg['sub_load_%s_path' % nutrient]
        subsurface_load_task = task_graph.add_task(
            func=_map_subsurface_load,
            args=(modified_load_path, f_reg['aligned_lulc_path'],
                  lucode_to_parameters,
                  subsurface_proportion_type, subsurface_load_path),
            target_path_list=[subsurface_load_path],
            dependent_task_list=[modified_load_task, align_raster_task],
            task_name='map subsurface load %s' % nutrient)

        eff_path = f_reg['eff_%s_path' % nutrient]
        eff_task = task_graph.add_task(
            func=_map_lulc_to_val_mask_stream,
            args=(
                f_reg['aligned_lulc_path'], f_reg['stream_path'],
                lucode_to_parameters, 'eff_%s' % nutrient, eff_path),
            target_path_list=[eff_path],
            dependent_task_list=[align_raster_task, stream_extraction_task],
            task_name='ret eff %s' % nutrient)

        crit_len_path = f_reg['crit_len_%s_path' % nutrient]
        crit_len_task = task_graph.add_task(
            func=_map_lulc_to_val_mask_stream,
            args=(
                f_reg['aligned_lulc_path'], f_reg['stream_path'],
                lucode_to_parameters, 'crit_len_%s' % nutrient, crit_len_path),
            target_path_list=[crit_len_path],
            dependent_task_list=[align_raster_task, stream_extraction_task],
            task_name='ret eff %s' % nutrient)

        effective_retention_path = (
            f_reg['effective_retention_%s_path' % nutrient])
        ndr_eff_task = task_graph.add_task(
            func=ndr_core.ndr_eff_calculation,
            args=(
                f_reg['flow_direction_path'], f_reg['stream_path'], eff_path,
                crit_len_path, effective_retention_path),
            target_path_list=[effective_retention_path],
            dependent_task_list=[
                stream_extraction_task, eff_task, crit_len_task],
            task_name='eff ret %s' % nutrient)

        ndr_path = f_reg['ndr_%s_path' % nutrient]
        ndr_task = task_graph.add_task(
            func=_calculate_ndr,
            args=(
                effective_retention_path, f_reg['ic_factor_path'],
                float(args['k_param']), ndr_path),
            target_path_list=[ndr_path],
            dependent_task_list=[ndr_eff_task, ic_task],
            task_name='calc ndr %s' % nutrient)

        sub_ndr_path = f_reg['sub_ndr_%s_path' % nutrient]
        sub_ndr_task = task_graph.add_task(
            func=_calculate_sub_ndr,
            args=(
                float(args['subsurface_eff_%s' % nutrient]),
                float(args['subsurface_critical_length_%s' % nutrient]),
                f_reg['dist_to_channel_path'], sub_ndr_path),
            target_path_list=[sub_ndr_path],
            dependent_task_list=[dist_to_channel_task],
            task_name='sub ndr %s' % nutrient)

        export_path = f_reg['%s_export_path' % nutrient]
        calculate_export_task = task_graph.add_task(
            func=_calculate_export,
            args=(
                surface_load_path, ndr_path, subsurface_load_path,
                sub_ndr_path, export_path),
            target_path_list=[export_path],
            dependent_task_list=[
                load_task, ndr_task, surface_load_task, subsurface_load_task,
                sub_ndr_task],
            task_name='export %s' % nutrient)

        aggregate_export_task = task_graph.add_task(
            func=_aggregate_and_pickle_total,
            args=(
                (export_path, 1), f_reg['watershed_results_ndr_path'],
                f_reg['export_%s_pickle_path' % nutrient]),
            target_path_list=[f_reg['export_%s_pickle_path' % nutrient]],
            dependent_task_list=[calculate_export_task],
            task_name='aggregate %s export' % nutrient)

        aggregate_surface_load_task = task_graph.add_task(
            func=_aggregate_and_pickle_total,
            args=(
                (surface_load_path, 1), f_reg['watershed_results_ndr_path'],
                f_reg['surface_load_%s_pickle_path' % nutrient]),
            target_path_list=[f_reg['surface_load_%s_pickle_path' % nutrient]],
            dependent_task_list=[surface_load_task, create_vector_task],
            task_name='aggregate %s surface load' % nutrient)

        aggregate_subsurface_load_task = task_graph.add_task(
            func=_aggregate_and_pickle_total,
            args=(
                (subsurface_load_path, 1), f_reg['watershed_results_ndr_path'],
                f_reg['subsurface_load_%s_pickle_path' % nutrient]),
            target_path_list=[
                f_reg['subsurface_load_%s_pickle_path' % nutrient]],
            dependent_task_list=[subsurface_load_task, create_vector_task],
            task_name='aggregate %s subsurface load' % nutrient)

        field_pickle_map['surf_%s_ld' % nutrient] = (
            f_reg['surface_load_%s_pickle_path' % nutrient])
        field_pickle_map['sub_%s_ld' % nutrient] = (
            f_reg['subsurface_load_%s_pickle_path' % nutrient])
        field_pickle_map['%s_exp_tot' % nutrient] = (
            f_reg['export_%s_pickle_path' % nutrient])
        field_header_order_list = (
            [x % nutrient for x in [
                'surf_%s_ld', 'sub_%s_ld', '%s_exp_tot']] +
            field_header_order_list)

    task_graph.close()
    task_graph.join()

    LOGGER.info('Writing summaries to output shapefile')
    _add_fields_to_shapefile(
        field_pickle_map, field_header_order_list,
        f_reg['watershed_results_ndr_path'])

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
        field_pickle_map, field_header_order, target_vector_path):
    """Add fields and values to an OGR layer open for writing.

    Parameters:
        field_pickle_map (dict): maps field name to a pickle file that is a
            result of pygeoprocessing.zonal_stats with FIDs that match
            `target_vector_path`.
        field_header_order (list of string): a list of field headers in the
            order to appear in the output table.
        target_vector_path (string): path to target vector file.

    Returns:
        None.

    """
    target_vector = gdal.OpenEx(
        target_vector_path, gdal.OF_VECTOR | gdal.GA_Update)
    target_layer = target_vector.GetLayer()
    field_summaries = {}
    for field_name in field_header_order:
        field_def = ogr.FieldDefn(field_name, ogr.OFTReal)
        field_def.SetWidth(24)
        field_def.SetPrecision(11)
        target_layer.CreateField(field_def)
        with open(field_pickle_map[field_name]) as pickle_file:
            field_summaries[field_name] = pickle.load(pickle_file)

    for feature in target_layer:
        fid = feature.GetFID()
        for field_name in field_header_order:
            feature.SetField(
                field_name, float(field_summaries[field_name][fid]['sum']))
        # Save back to datasource
        target_layer.SetFeature(feature)
    target_layer = None
    target_vector = None


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
                    raster = gdal.OpenEx(args[key], gdal.OF_RASTER)
                    if raster is None:
                        validation_error_list.append(
                            ([key], 'not a raster'))
                    del raster
                elif key_type == 'vector':
                    vector = gdal.OpenEx(args[key], gdal.OF_VECTOR)
                    if vector is None:
                        validation_error_list.append(
                            ([key], 'not a vector'))
                    del vector

    return validation_error_list


def _normalize_raster(base_raster_path_band, target_normalized_raster_path):
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
        value_sum += numpy.sum(valid_block)
        value_count += valid_block.size

    value_mean = value_sum
    if value_count > 0.0:
        value_mean /= value_count

    def _normalize_raster_op(array):
        """Divide values by mean."""
        valid_mask = ~numpy.isclose(array, base_nodata)
        result = numpy.empty(valid_mask.shape, dtype=numpy.float32)
        result[:] = base_nodata
        result[valid_mask] = array[valid_mask]
        if value_mean != 0:
            result[valid_mask] /= value_mean
        return result

    pygeoprocessing.raster_calculator(
        [base_raster_path_band], _normalize_raster_op,
        target_normalized_raster_path, gdal.GDT_Float32, base_nodata)


def _calculate_load(
        lulc_raster_path, lucode_to_parameters, load_type,
        target_load_raster):
    """Calculate load raster by mapping landcover and multiplying by area.

    Parameters:
        lulc_raster_path (string): path to integer landcover raster.
        lucode_to_parameters (dict): a mapping of landcover IDs to a
            dictionary indexed by the value of `load_{load_type}` that
            represents a per-area nutrient load.
        load_type (string): represent nutrient to map, either 'load_n' or
            'load_p'.
        target_load_raster (string): path to target raster that will have
            total load per pixel.

    Returns:
        None.

    """
    lulc_raster_info = pygeoprocessing.get_raster_info(lulc_raster_path)
    nodata_landuse = lulc_raster_info['nodata'][0]
    cell_area_ha = abs(numpy.prod(lulc_raster_info['pixel_size'])) * 0.0001

    def _map_load_op(lucode_array):
        """Convert unit load to total load & handle nodata."""
        result = numpy.empty(lucode_array.shape)
        result[:] = _TARGET_NODATA
        for lucode in numpy.unique(lucode_array):
            if lucode != nodata_landuse:
                    result[lucode_array == lucode] = (
                        lucode_to_parameters[lucode][load_type] *
                        cell_area_ha)
        return result

    pygeoprocessing.raster_calculator(
        [(lulc_raster_path, 1)], _map_load_op, target_load_raster,
        gdal.GDT_Float32, _TARGET_NODATA)


def _multiply_rasters(raster_path_list, target_nodata, target_result_path):
    """Multiply the rasters in `raster_path_list`.

    Parameters:
        raster_path_list (list): list of single band raster paths.
        target_nodata (float): desired target nodata value.
        target_result_path (string): path to float 32 target raster
            multiplied where all rasters are not nodata.

    Returns:
        None.

    """
    def _mult_op(*array_nodata_list):
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
        path_nodata_list, _mult_op, target_result_path,
        gdal.GDT_Float32, target_nodata)


def _map_surface_load(
        modified_load_path, lulc_raster_path, lucode_to_parameters,
        subsurface_proportion_type, target_surface_load_path):
    """Calculate surface load from landcover raster.

    Parameters:
        modified_load_path (string): path to modified load raster with units
            of kg/pixel.
        lulc_raster_path (string): path to landcover raster.
        lucode_to_parameters (dict): maps landcover codes to a dictionary that
            can be indexed by `subsurface_proportion_type`.
        subsurface_proportion_type (string): if None no subsurface transfer
            is mapped.  Otherwise indexed from lucode_to_parameters.
        target_surface_load_path (string): path to target raster.

    Returns:
        None.

    """
    lulc_raster_info = pygeoprocessing.get_raster_info(lulc_raster_path)
    nodata_landuse = lulc_raster_info['nodata'][0]

    keys = sorted(numpy.array(lucode_to_parameters.keys()))
    if subsurface_proportion_type is not None:
        subsurface_values = numpy.array(
            [lucode_to_parameters[x][subsurface_proportion_type]
             for x in keys])

    def _map_surface_load_op(lucode_array, modified_load_array):
        """Convert unit load to total load & handle nodata."""
        # If we don't have subsurface, just return 0.0.
        if subsurface_proportion_type is None:
            return numpy.where(
                lucode_array != nodata_landuse, modified_load_array,
                _TARGET_NODATA)

        valid_mask = lucode_array != nodata_landuse
        result = numpy.empty(valid_mask.shape, dtype=numpy.float32)
        result[:] = _TARGET_NODATA
        index = numpy.digitize(
            lucode_array[valid_mask].ravel(), keys, right=True)
        result[valid_mask] = (
            modified_load_array[valid_mask] * (1 - subsurface_values[index]))
        return result

    pygeoprocessing.raster_calculator(
        [(lulc_raster_path, 1), (modified_load_path, 1)],
        _map_surface_load_op, target_surface_load_path, gdal.GDT_Float32,
        _TARGET_NODATA)


def _map_subsurface_load(
        modified_load_path, lulc_raster_path, lucode_to_parameters,
        subsurface_proportion_type, target_sub_load_path):
    """Calculate subsurface load from landcover raster.

    Parameters:
        modified_load_path (string): path to modified load raster.
        lulc_raster_path (string): path to landcover raster.
        lucode_to_parameters (dict): maps landcover codes to a dictionary that
            can be indexed by by `subsurface_proportion_type`.
        subsurface_proportion_type (string): if None no subsurface transfer
            is mapped.  Otherwise indexed from lucode_to_parameters.
        target_sub_load_path (string): path to target raster.

    Returns:
        None.

    """
    lulc_raster_info = pygeoprocessing.get_raster_info(lulc_raster_path)
    nodata_landuse = lulc_raster_info['nodata'][0]

    keys = sorted(numpy.array(lucode_to_parameters.keys()))
    if subsurface_proportion_type is not None:
        subsurface_permeance_values = numpy.array(
            [lucode_to_parameters[x][subsurface_proportion_type]
             for x in keys])

    def _map_subsurface_load_op(lucode_array, modified_load_array):
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
            modified_load_array[valid_mask] *
            subsurface_permeance_values[index])
        return result

    pygeoprocessing.raster_calculator(
        [(lulc_raster_path, 1), (modified_load_path, 1)],
        _map_subsurface_load_op, target_sub_load_path, gdal.GDT_Float32,
        _TARGET_NODATA)


def _map_lulc_to_val_mask_stream(
        lulc_raster_path, stream_path, lucode_to_parameters, map_id,
        target_eff_path):
    """Make retention efficiency raster from landcover.

    Parameters:
        lulc_raster_path (string): path to landcover raster.
        stream_path (string) path to stream layer 0, no stream 1 stream.
        lucode_to_parameters (dict) mapping of landcover code to a dictionary
            that contains the key in `map_id`
        map_id (string): the id in the lookup table with values to map
            landcover to efficiency.
        target_eff_path (string): target raster that contains the mapping of
            landcover codes to retention efficiency values except where there
            is a stream in which case the retention efficiency is 0.

    Returns:
        None.

    """
    keys = sorted(numpy.array(lucode_to_parameters.keys()))
    values = numpy.array(
        [lucode_to_parameters[x][map_id] for x in keys])

    nodata_landuse = pygeoprocessing.get_raster_info(
        lulc_raster_path)['nodata'][0]
    nodata_stream = pygeoprocessing.get_raster_info(stream_path)['nodata'][0]

    def _map_eff_op(lucode_array, stream_array):
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
        ((lulc_raster_path, 1), (stream_path, 1)), _map_eff_op,
        target_eff_path, gdal.GDT_Float32, _TARGET_NODATA)


def s_bar_calculate(
        s_accumulation_path, flow_accumulation_path, target_s_bar_path):
    """Calculate bar op which is s/flow."""
    s_nodata = pygeoprocessing.get_raster_info(
        s_accumulation_path)['nodata'][0]
    flow_nodata = pygeoprocessing.get_raster_info(
        flow_accumulation_path)['nodata'][0]

    def _bar_op(s_accumulation, flow_accumulation):
        """Calculate bar operation of s_accum / flow_accum."""
        valid_mask = (
            (s_accumulation != s_nodata) &
            (flow_accumulation != flow_nodata))
        result = numpy.empty(valid_mask.shape, dtype=numpy.float32)
        result[:] = _TARGET_NODATA
        result[valid_mask] = (
            s_accumulation[valid_mask] / flow_accumulation[valid_mask])
        return result

    pygeoprocessing.raster_calculator(
        ((s_accumulation_path, 1), (flow_accumulation_path, 1)), _bar_op,
        target_s_bar_path, gdal.GDT_Float32, _TARGET_NODATA)


def d_up_calculation(s_bar_path, flow_accum_path, target_d_up_path):
    """Calculate d_up = s_bar * sqrt(upstream area)."""
    s_bar_info = pygeoprocessing.get_raster_info(s_bar_path)
    s_bar_nodata = s_bar_info['nodata'][0]
    flow_accum_nodata = pygeoprocessing.get_raster_info(
        flow_accum_path)['nodata'][0]
    cell_area_m2 = abs(numpy.prod(s_bar_info['pixel_size']))

    def _d_up_op(s_bar, flow_accumulation):
        """Calculate d_up index."""
        valid_mask = (
            (s_bar != s_bar_nodata) &
            (flow_accumulation != flow_accum_nodata))
        result = numpy.empty(valid_mask.shape, dtype=numpy.float32)
        result[:] = _TARGET_NODATA
        result[valid_mask] = (
            s_bar[valid_mask] * numpy.sqrt(
                flow_accumulation[valid_mask] * cell_area_m2))
        return result

    pygeoprocessing.raster_calculator(
        [(s_bar_path, 1), (flow_accum_path, 1)], _d_up_op,
        target_d_up_path, gdal.GDT_Float32, _TARGET_NODATA)


def invert_raster_values(base_raster_path, target_raster_path):
    """Invert (1/x) the values in `base`.

    Parameters:
        base_raster_path (string): path to floating point raster.
        target_raster_path (string): path to created output raster whose
            values are 1/x of base.

    Returns:
        None.

    """
    base_nodata = pygeoprocessing.get_raster_info(
        base_raster_path)['nodata'][0]

    def _inverse_op(base_val):
        """Calculate inverse of S factor."""
        result = numpy.empty(base_val.shape, dtype=numpy.float32)
        result[:] = _TARGET_NODATA
        valid_mask = ~numpy.isclose(base_val, base_nodata)
        zero_mask = base_val == 0.0
        result[valid_mask & ~zero_mask] = (
            1.0 / base_val[valid_mask & ~zero_mask])
        result[zero_mask] = 0.0
        return result

    pygeoprocessing.raster_calculator(
        ((base_raster_path, 1),), _inverse_op,
        target_raster_path, gdal.GDT_Float32, _TARGET_NODATA)


def calculate_ic(d_up_path, d_dn_path, target_ic_path):
    """Calculate IC as log_10(d_up/d_dn)."""
    ic_nodata = float(numpy.finfo(numpy.float32).min)
    d_up_nodata = pygeoprocessing.get_raster_info(d_up_path)['nodata'][0]
    d_dn_nodata = pygeoprocessing.get_raster_info(d_dn_path)['nodata'][0]

    def _ic_op(d_up, d_dn):
        """Calculate IC0."""
        valid_mask = (
            (d_up != d_up_nodata) & (d_dn != d_dn_nodata) & (d_up != 0) &
            (d_dn != 0))
        result = numpy.empty(valid_mask.shape, dtype=numpy.float32)
        result[:] = ic_nodata
        result[valid_mask] = numpy.log10(d_up[valid_mask] / d_dn[valid_mask])
        return result

    pygeoprocessing.raster_calculator(
        [(d_up_path, 1), (d_dn_path, 1)], _ic_op,
        target_ic_path, gdal.GDT_Float32, ic_nodata)


def _calculate_ndr(
        effective_retention_path, ic_factor_path, k_param, target_ndr_path):
    """Calculate NDR as a function of Equation 4 in the user's guide."""
    ic_factor_raster = gdal.OpenEx(ic_factor_path, gdal.OF_RASTER)
    ic_factor_band = ic_factor_raster.GetRasterBand(1)
    ic_min, ic_max, _, _ = ic_factor_band.GetStatistics(0, 1)
    ic_factor_band = None
    ic_factor_raster = None
    ic_0_param = (ic_min + ic_max) / 2.0
    effective_retention_nodata = pygeoprocessing.get_raster_info(
        effective_retention_path)['nodata'][0]
    ic_nodata = pygeoprocessing.get_raster_info(ic_factor_path)['nodata'][0]

    def _calculate_ndr_op(effective_retention_array, ic_array):
        """Calculate NDR."""
        valid_mask = (
            (effective_retention_array != effective_retention_nodata) &
            (ic_array != ic_nodata))
        result = numpy.empty(valid_mask.shape, dtype=numpy.float32)
        result[:] = _TARGET_NODATA
        result[valid_mask] = (
            (1.0 - effective_retention_array[valid_mask]) /
            (1.0 + numpy.exp(
                (ic_0_param - ic_array[valid_mask]) / k_param)))
        return result

    pygeoprocessing.raster_calculator(
        [(effective_retention_path, 1), (ic_factor_path, 1)],
        _calculate_ndr_op, target_ndr_path, gdal.GDT_Float32, _TARGET_NODATA)


def _calculate_sub_ndr(
        eff_sub, crit_len_sub, dist_to_channel_path, target_sub_ndr_path):
    """Calculate subsurface: subndr = eff_sub(1-e^(-5*l/crit_len)."""
    dist_to_channel_nodata = pygeoprocessing.get_raster_info(
        dist_to_channel_path)['nodata'][0]

    def _sub_ndr_op(dist_to_channel_array):
        """Calculate subsurface NDR."""
        valid_mask = ~numpy.isclose(
            dist_to_channel_array, dist_to_channel_nodata)
        result = numpy.empty(valid_mask.shape, dtype=numpy.float32)
        result[:] = _TARGET_NODATA
        result[valid_mask] = 1.0 - eff_sub * (
            1-numpy.exp(-5*dist_to_channel_array[valid_mask]/crit_len_sub))
        return result

    pygeoprocessing.raster_calculator(
        [(dist_to_channel_path, 1)], _sub_ndr_op, target_sub_ndr_path,
        gdal.GDT_Float32, _TARGET_NODATA)


def _calculate_export(
        surface_load_path, ndr_path, subsurface_load_path,
        subsurface_ndr_path, target_export_path):
    """Calculate export."""
    load_nodata = pygeoprocessing.get_raster_info(
        surface_load_path)['nodata'][0]
    subsurface_load_nodata = pygeoprocessing.get_raster_info(
        subsurface_load_path)['nodata'][0]
    ndr_nodata = pygeoprocessing.get_raster_info(
        ndr_path)['nodata'][0]
    sub_ndr_nodata = pygeoprocessing.get_raster_info(
        subsurface_ndr_path)['nodata'][0]

    def _calculate_export_op(
            modified_load_array, ndr_array, modified_sub_load_array,
            sub_ndr_array):
        """Combine NDR and subsurface NDR."""
        valid_mask = (
            ~numpy.isclose(modified_load_array, load_nodata) &
            ~numpy.isclose(ndr_array, ndr_nodata) &
            ~numpy.isclose(modified_sub_load_array, subsurface_load_nodata) &
            ~numpy.isclose(sub_ndr_array, sub_ndr_nodata))
        result = numpy.empty(valid_mask.shape, dtype=numpy.float32)
        result[:] = _TARGET_NODATA
        result[valid_mask] = (
            modified_load_array[valid_mask] * ndr_array[valid_mask] +
            modified_sub_load_array[valid_mask] *
            sub_ndr_array[valid_mask])
        return result

    pygeoprocessing.raster_calculator(
        [(surface_load_path, 1), (ndr_path, 1),
         (subsurface_load_path, 1), (subsurface_ndr_path, 1)],
        _calculate_export_op, target_export_path, gdal.GDT_Float32,
        _TARGET_NODATA)


def _aggregate_and_pickle_total(
        base_raster_path_band, aggregate_vector_path, target_pickle_path):
    """Aggregate base raster path to vector path FIDs and pickle result.

    Parameters:
        base_raster_path_band (tuple): raster/path band to aggregate over.
        aggregate_vector_path (string): path to vector to use geometry to
            aggregate over.
        target_pickle_path (string): path to a file that will contain the
            result of a pygeoprocessing.zonal_statistics call over
            base_raster_path_band from aggregate_vector_path.

    Returns:
        None.

    """
    result = pygeoprocessing.zonal_statistics(
        base_raster_path_band, aggregate_vector_path,
        working_dir=os.path.dirname(target_pickle_path))

    with open(target_pickle_path, 'w') as target_pickle_file:
        pickle.dump(result, target_pickle_file)


def create_vector_copy(base_vector_path, target_vector_path):
    """Create a copy of base vector."""
    if os.path.isfile(target_vector_path):
        os.remove(target_vector_path)
    base_vector = gdal.OpenEx(base_vector_path, gdal.OF_VECTOR)
    driver = gdal.GetDriverByName('ESRI Shapefile')
    target_vector = driver.CreateCopy(
        target_vector_path, base_vector)
    target_vector = None  # seemingly uncessary but gdal seems to like it.
