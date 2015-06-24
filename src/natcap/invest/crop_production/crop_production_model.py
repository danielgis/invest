'''
The Crop Production Model module contains functions for running the model
'''

import os
import logging

import pygeoprocessing.geoprocessing as pygeo

LOGGER = logging.getLogger('CROP_PRODUCTION')
logging.basicConfig(format='%(asctime)s %(name)-15s %(levelname)-8s \
    %(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')


def create_observed_yield_maps(vars_dict):
    '''
    Creates yield maps from observed data and stores in the temporary yield
        folder

    Args:
        vars_dict (dict): descr

    Example Args::

        vars_dict = {
            # ...

            'lulc_map_uri': 'path/to/lulc_map_uri',
            'crop_lookup_dict': {
                'code': 'crop_name',
                ...
            },
            'observed_yields_maps_dict': {
                'crop': 'path/to/crop_climate_bin_map',
                ...
            },
            'tmp_observed_dir': '/path/to/tmp_observed_dir',
        }

    Outputs:

    .
    |-- tmp
        |-- observed
            |-- yield
                |-- [crop]_yield_map (*.tif)
    '''
    def lulc_mask_over_yield_map():
        pass

    # Make yield folder

    lulc_map_uri = vars_dict['lulc_map_uri']
    # AOI from lulc_map
    lulc_map_aoi = pygeo.get_bounding_box(lulc_map_uri)  # Will this work?

    for crop in vars_dict['observed_yields_maps_dict'].keys():
        observed_yield_map_uri = vars_dict['observed_yields_maps_dict'][crop]

        # Transform lulc_map AOI into observed_yield_map projection
        pygeo.reproject_dataset_uri()

        # Clip observed_yields_map over AOI: clipped_observed_yield_map_uri

        # Reproject and resample clipped_observed_yield_map to match lulc_map

        # Perform masked operation to produce observed_yields_map
        dataset_uri_list = [lulc_map_uri, clipped_observed_yield_map_uri]
        dataset_pixel_op = lulc_mask_over_yield_map
        dataset_out_uri = os.path.join(
            tmp_observed_dir,
            'yield/' + crop + '_yield_map.tif')
        gdal_dtype = pygeo.get_datatype_from_uri(observed_yield_map_uri)
        nodata = -1
        cell_size = pygeo.get_cell_size_from_uri(lulc_map_uri)

        pygeo.vectorize_datasets(
            dataset_uri_list,
            dataset_pixel_op,
            dataset_out_uri,
            gdal_dtype,
            nodata,
            cell_size,
            "union")

    # For crop in observed_yields_maps_dict:
    #    

    # For Each Crop, Clip Corresponding Observed Crop Yield Map over AOI and Reproject
    #   Save to temporary directory

    # Create Crop Production Maps by Multiplying Yield by Cell Size Area
    #   Output: Crop Production Maps

    # Find Total Production for Given Crop by Summing Cells in Crop Production Maps
    #   Output: Crop Production Dictionary?

    pass


def create_percentile_yield_maps(vars_dict):
    '''
    About

    var_name (type): desc

    Example Args::

        vars_dict = {
            ...

            '': '',

            ...
        }

    Example Returns::

        vars_dict = {
            ...

            '': '',

            ...
        }
    '''
    # Get List of Crops in LULC Crop Map

    # Create Raster of Climate Bin Indices

    # For Each Yield Column in Percentile Yield Table:

    # For Each Crop, Create Crop Yield Map over AOI
        # Output: Crop Yield Maps

    # Create Crop Production Maps by Multiplying Yield by Cell Size Area
        # Output: Crop Production Maps
        # If 'create_crop_production_maps' selected, save to output folder

    # Find Total Production for Given Crop by Summing Cells in Crop Production Maps
        # Output: Crop Production Dictionary?

    # Generate Yield Results

    pass


def create_regression_yield_maps(vars_dict):
    '''
    About

    var_name (type): desc

    Example Args::

        vars_dict = {
            ...

            '': '',

            ...
        }

    Example Returns::

        vars_dict = {
            ...

            '': '',

            ...
        }
    '''
    pass


def create_production_maps(vars_dict):
    '''
    About

    var_name (type): desc

    Example Args::

        vars_dict = {
            ...

            '': '',

            ...
        }

    Example Returns::

        vars_dict = {
            ...

            '': '',

            ...
        }
    '''
    pass


def create_economic_returns_map(vars_dict):
    '''
    About

    var_name (type): desc

    Example Args::

        vars_dict = {
            ...

            '': '',

            ...
        }

    Example Returns::

        vars_dict = {
            ...

            '': '',

            ...
        }
    '''
    pass


def create_results_table(vars_dict):
    '''
    About

    var_name (type): desc

    Example Args::

        vars_dict = {
            ...

            '': '',

            ...
        }

    Example Returns::

        vars_dict = {
            ...

            '': '',

            ...
        }
    '''
    pass


# Raster Utils Wrapper for High-level Functions
def clip_raster_over_aoi():
    '''
    About

    var_name (type): desc

    Example Args::

        vars_dict = {
            ...

            '': '',

            ...
        }

    Example Returns::

        vars_dict = {
            ...

            '': '',

            ...
        }
    '''
    # Check that datasets are in same projection
    pygeo.assert_datasets_in_same_projection()

    # Reproject AOI onto Raster, find bounding box
    pygeo.reproject_dataset_uri()

    # Clip raster around bounding box
    pygeo.clip_dataset_uri()

    # Reproject Clipped Raster onto AOI
    pygeo.reproject_dataset_uri()

    pygeo.align_dataset_list()

    # Save/return clipped raster
    pass


def sum_cells_in_raster():
    '''
    About

    var_name (type): desc

    Example Args::

        vars_dict = {
            ...

            '': '',

            ...
        }

    Example Returns::

        vars_dict = {
            ...

            '': '',

            ...
        }
    '''
    # Option 1: use sum or mean in GDAL's raster stats
    # Option 2: extract numpy array using GDAL
    # Option 3: use Raster_Utils functionality memmap for large rasters

    pass


def element_wise_operation():
    '''
    About

    var_name (type): desc

    Example Args::

        vars_dict = {
            ...

            '': '',

            ...
        }

    Example Returns::

        vars_dict = {
            ...

            '': '',

            ...
        }
    '''
    pygeo.vectorize_datasets(
        dataset_uri_list,
        dataset_pixel_op,
        dataset_out_uri,
        datatype_out,
        nodata_out,
        pixel_size_out,
        bounding_box_mode)
    pass
