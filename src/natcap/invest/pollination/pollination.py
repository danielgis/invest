"""Pollinator service model for InVEST."""

import os
import logging
import re
import shutil
import collections

import numpy
from osgeo import osr
from osgeo import gdal
import pygeoprocessing

from natcap.invest import fileio as fileio

logging.basicConfig(format='%(asctime)s %(name)-18s %(levelname)-8s \
     %(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

LOGGER = logging.getLogger('natcap.invest.pollination.pollination')

def execute(args):
    """
    Execute the pollination model from the topmost, user-accessible level.

    Parameters:
        args['workspace_dir'] (string): a URI to the workspace folder.  Not
            required to exist on disk.  Additional folders will be created
            inside of this folder.  If there are any file name collisions, this
            model will overwrite those files.
        args['landuse_uri'] (string): a URI to a GDAL raster on disk that
            represents the landcover map.
        args['do_valuation'] - A boolean.  Indicates whether valuation should be
            performed.  This applies to all scenarios.
        args['landuse_attributes_uri'] (string): a URI to a CSV on disk.  See
            the model's documentation for details on the structure of this
            table.
        args['do_valuation'] (boolean): Indicates whether the model should
            include valuation
        args['half_saturation'] (float): a number between 0 and 1 indicating the
            half-saturation constant. See the pollination documentation for
            more information.
        args['wild_pollination_proportion'] (float): a number between 0 and 1
            indicating the proportion of all pollinators that are wild.
            See the pollination documentation for more information.
        args['guilds_uri'] (string): a URI to a CSV on disk.  See the model's
            documentation for details on the structure of this table.
        args['ag_classes'] (string): (Optional) a space-separated list of land
            cover classes that are to be considered as agricultural.  If this
            input is not provided, all land cover classes are considered to
            be agricultural.
        args['farms_shapefile'] (string): (Optional) shapefile containing points
            representing data collection points on the landscape.
        args['results_suffix'] (string): inserted into the URI of each file
            created by this model, right before the file extension.

    Example Args Dictionary::

        {
            'workspace_dir': 'path/to/workspace_dir',
            'landuse_cur_uri': 'path/to/raster',
            'landuse_attributes_uri': 'path/to/csv',
            'do_valuation': 'example',
            'half_saturation': 'example',
            'wild_pollination_proportion': 'example',
            'guilds_uri': 'path/to/csv',
            'ag_classes': 'example',
            'results_suffix': 'example',

        }

    The following args dictionary entries are optional, and
    will affect the behavior of the model if provided:
        1. ag_classes
        2. results_suffix

    If args['do_valuation'] is set to True, the following args dictionary
    entries are also required:
        1. half_saturation
        2. wild_pollination_proportion

    This function has no return value, though it does save a number of
    rasters to disk.  See the user's guide for details.

    """

    #append a _ to the suffix if it's not empty and doens't already have one
    try:
        file_suffix = args['results_suffix']
        if file_suffix != "" and not file_suffix.startswith('_'):
            file_suffix = '_' + file_suffix
    except KeyError:
        file_suffix = ''

    workspace = args['workspace_dir']

    # Check to see if each of the workspace folders exists.  If not, create the
    # folder in the filesystem.
    inter_dir = os.path.join(workspace, 'intermediate')
    out_dir = os.path.join(workspace, 'output')
    pygeoprocessing.create_directories([inter_dir, out_dir])

    LOGGER.info('Starting pollination model')

    # Open a Table Handler for the land use attributes table and a different
    # table handler for the Guilds table.
    LOGGER.info('Opening landuse attributes table')

    land_attribute_table = pygeoprocessing.get_lookup_from_table(
        args['landuse_attributes_uri'], 'lulc')
    attribute_table_fields = land_attribute_table.itervalues().next().keys()
    nesting_fields = [
        f[2:] for f in attribute_table_fields if re.match('^n_', f)]
    floral_fields = [
        f[2:] for f in attribute_table_fields if re.match('^f_', f)]
    fields_to_check = [(nesting_fields, 'nesting'), (floral_fields, 'floral')]
    for field_list, field_type in fields_to_check:
        if len(field_list) == 0:
            raise ValueError(
                'LULC attribute table must have '
                ' %s fields but none were found.' % field_type)

    LOGGER.info('Opening guilds table')
    guilds_table = pygeoprocessing.get_lookup_from_table(
        args['guilds_uri'], 'species')

    # Convert agricultural classes (a space-separated list of ints) into a
    # list of ints.  If the user has not provided a string list of ints,
    # then use an empty list instead.
    LOGGER.info('Processing agricultural classes')
    try:
        # This approach will create a list with only ints, even if the user
        # has accidentally entered additional spaces.  Any other incorrect
        # input will throw a ValueError exception.
        user_ag_list = args['ag_classes'].split(' ')
        ag_class_list = [int(r) for r in user_ag_list if r != '']
    except KeyError:
        # If the 'ag_classes' key is not present in the args dictionary,
        # use an empty list in its stead.
        ag_class_list = []

    # Defined which rasters need to be created at the global level (at the
    # top level of the model dictionary).  the global_rasters list has this
    # structure:
    #   (model_args key, raster_uri base, folder to be saved to)
    global_rasters = [
        ('foraging_total', 'frm_tot', out_dir),
        ('foraging_average', 'frm_avg', out_dir),
        ('farm_value_sum', 'frm_val_sum', inter_dir),
        ('service_value_sum', 'sup_val_sum', out_dir),
        ('abundance_total', 'sup_tot', out_dir),
        ('ag_mask', 'agmap', inter_dir)]

    # loop through the global rasters provided and actually create the uris,
    # saving them to the model args dictionary.
    LOGGER.info('Creating top-level raster URI paths')
    global_raster_uris = {}
    for raster_id, basename, directory in global_rasters:
        raster_uri = os.path.join(
            directory, '%s%s.tif' % (basename, file_suffix))
        global_raster_uris[raster_id] = raster_uri

    # Fetch a list of all species from the guilds table.
    species_list = guilds_table.keys()

    # Make new rasters for each species.  In this list of tuples, the first
    # value of each tuple is the args dictionary key, and the second value
    # of each tuple is the raster prefix.  The third value is the folder in
    # which the created raster should exist.
    species_rasters = [
        ('nesting', 'hn', inter_dir),
        ('floral', 'hf', inter_dir),
        ('species_abundance', 'sup', inter_dir),
        ('farm_abundance', 'frm', inter_dir),
        ('farm_value', 'frm_val', inter_dir),
        ('value_abundance_ratio', 'val_sup_ratio', inter_dir),
        ('value_abundance_ratio_blur', 'val_sup_ratio_blur', inter_dir),
        ('service_value', 'sup_val', inter_dir)]

    # Loop through each species and define the necessary raster URIs, as
    # defined by the species_rasters list.
    LOGGER.info('Creating species-specific raster URIs')
    species_raster_uris = collections.defaultdict(dict)
    for species in species_list:
        if not isinstance(species, unicode):
            # Casting to UTF-8 so that LOGGER won't crash if species is UTF-8
            species = unicode(species, "utf-8")
        LOGGER.info('Creating rasters for %s', species)
        for group, prefix, folder in species_rasters:
            raster_uri = os.path.join(
                folder, '%s_%s%s.tif' % (prefix, species, file_suffix))
            species_raster_uris[species][group] = raster_uri

    create_ag_mask(
        args['landuse_uri'], global_raster_uris['ag_mask'], args['ag_classes'])

    # Loop through all species and perform the necessary calculations.
    for species, species_dict in guilds_table.iteritems():
        LOGGER.debug(species)
        LOGGER.debug(species_dict)

        # Calculate species abundance.  This represents the relative index of
        # how much of a species we can expect to find across the landscape given
        # the floral and nesting patterns (based on land cover) and the
        # specified use of these resources (defined in the guild_dict).
        LOGGER.info('Calculating %s abundance on the landscape', species)
        LOGGER.debug(land_attribute_table)
        calculate_species_abundance_index(
            args['landuse_uri'], land_attribute_table, species_dict,
            nesting_fields, floral_fields, uris={
                'nesting': species_dict['nesting'],
                'floral': species_dict['floral'],
                'species_abundance': species_dict['species_abundance'],
                'temp': args['paths']['temp']
            })

        # Add the newly-calculated abundance to the abundance_sum raster.
        LOGGER.info('Adding %s species abundance to the total', species)
        add_two_rasters(
            args['abundance_total'], species_dict['species_abundance'],
            args['abundance_total'])

        # Calculate the farm abundance.  This takes the species abundance and
        # calculates roughly how much of a species we can expect to find on farm
        # pixels.
        LOGGER.info('Calculating %s abundance on farms ("foraging")', species)
        calculate_farm_abundance(
            species_dict['species_abundance'],
            args['ag_mask'], guild_dict['alpha'], species_dict['farm_abundance'])

        # Add the newly calculated farm abundance raster to the total.
        LOGGER.info('Adding %s foraging abundance raster to total', species)
        add_two_rasters(
            species_dict['farm_abundance'], args['foraging_average'],
            args['foraging_average'])

        if args['do_valuation'] == True:
            LOGGER.info('Starting species-specific valuation for %s', species)

            # Apply the half-saturation yield function from the documentation
            # and write it to its raster
            LOGGER.info('Calculating crop yield due to %s', species)
            calculate_yield(
                species_dict['farm_abundance'], species_dict['farm_value'],
                args['half_saturation'], args['wild_pollination_proportion'],
                -1.0)

            # Add the new farm_value_matrix to the farm value sum matrix.
            LOGGER.info(
                'Adding crop yield due to %s to the crop yield total', species)
            add_two_rasters(
                args['farm_value_sum'], species_dict['farm_value'],
                args['farm_value_sum'])

            LOGGER.info('Calculating service value for %s', species)
            guild_dict = args['guilds'].get_table_row('species', species)

            calculate_service(
                rasters={
                    'farm_value': species_dict['farm_value'],
                    'farm_abundance': species_dict['farm_abundance'],
                    'species_abundance': species_dict['species_abundance'],
                    'ag_mask': args['ag_mask']
                },
                nodata=-1.0,
                alpha=float(guild_dict['alpha']),
                part_wild=args['wild_pollination_proportion'],
                out_uris={
                    'species_value': species_dict['value_abundance_ratio'],
                    'species_value_blurred':\
                        species_dict['value_abundance_ratio_blur'],
                    'service_value': species_dict['service_value'],
                    'temp': args['paths']['temp']
                })

            # Add the new service value to the service value sum matrix
            LOGGER.info(
                'Adding the %s service value raster to the sum', species)
            add_two_rasters(
                args['service_value_sum'], species_dict['service_value'],
                args['service_value_sum'])

    # Calculate the average foraging index based on the total
    # Divide the total pollination foraging index by the number of pollinators
    # to get the mean pollinator foraging index and save that to its raster.
    num_species = float(len(args['species'].values()))

    # Calculate the mean foraging values per species.
    divide_raster(
        args['foraging_average'], num_species, args['foraging_average'])

    # Calculate the mean pollinator supply (pollinator abundance) by taking the
    # abundance_total_matrix and dividing it by the number of pollinators.
    divide_raster(args['abundance_total'], num_species, args['abundance_total'])


def calculate_species_abundance_index(
        landuse_uri, nesting_suitability_table, floral_resources_table,
        species_abundance_uri):
    """Calculate pollinator abundance on the landscape.  The calculated
    pollinator abundance raster will be created at `species_abundance_uri`

    Parameters:
        landuse_uri (string): a path to a GDAL dataset of the LULC
        nesting_suitability_table(dict): maps lucodes from `landuse_uri` to


        lu_attr - a TableHandler
        guild - a dictionary containing information about the pollinator.  All
            entries are required:
            'alpha' - the typical foraging distance in m
            'species_weight' - the relative weight
            resource_n - One entry for each nesting field for each fieldname
                denoted in nesting_fields.  This value must be either 0 and 1,
                indicating whether the pollinator uses this nesting resource for
                nesting sites.
            resource_f - One entry for each floral field denoted in
                floral_fields.  This value must be between 0 and 1,
                representing the liklihood that this species will forage during
                this season.
        nesting_fields - a list of string fieldnames.  Used to extract nesting
            fields from the guild dictionary, so fieldnames here must exist in
            guild.
        floral_fields - a list of string fieldnames.  Used to extract floral
            season fields from the guild dictionary.  Fieldnames here must exist
            in guild.
        uris - a dictionary with these entries:
            'nesting' - a URI to where the nesting raster will be saved.
            'floral' - a URI to where the floral resource raster will be saved.
            'species_abundance' - a URI to where the species abundance raster
                will be saved.
            'temp' - a URI to a folder where temp files will be saved

        Returns nothing."""
    nodata = -1.0


    floral_raster_temp_uri = pygeoprocessing.geoprocessing.temporary_filename()
    map_attribute(
        landuse_uri, lu_attr, guild, floral_fields, floral_raster_temp_uri, sum)
    map_attribute(landuse_uri, lu_attr, guild, nesting_fields, uris['nesting'], max)

    # Now that the per-pixel nesting and floral resources have been
    # calculated, the floral resources still need to factor in
    # neighborhoods.
    lulc_ds = gdal.Open(landuse_uri)
    pixel_size = abs(lulc_ds.GetGeoTransform()[1])
    lulc_ds = None
    expected_distance = guild['alpha'] / pixel_size
    kernel_uri = pygeoprocessing.geoprocessing.temporary_filename()
    make_exponential_decay_kernel_uri(expected_distance, kernel_uri)

    # Fetch the floral resources raster and matrix from the args dictionary
    # apply an exponential convolution filter and save the floral resources raster to the
    # dataset.
    pygeoprocessing.geoprocessing.convolve_2d_uri(
        floral_raster_temp_uri, kernel_uri, uris['floral'])
    os.remove(kernel_uri)

    # Calculate the pollinator abundance index (using Math! to simplify the
    # equation in the documentation.  We're still waiting on Taylor
    # Rickett's reply to see if this is correct.
    # Once the pollination supply has been calculated, we add it to the
    # total abundance matrix.
    try:
        species_weight = float(guild['species_weight'])
    except KeyError:
        # If the species_weight field has not been provided, assume that all
        # species weights should be equal (1.0).
        species_weight = 1.0

    pygeoprocessing.geoprocessing.vectorize_datasets(
        [uris['nesting'], uris['floral']],
        lambda x, y: numpy.where(x != nodata, x * y * species_weight, nodata),
        dataset_out_uri=uris['species_abundance'],
        datatype_out=gdal.GDT_Float32, nodata_out=nodata,
        pixel_size_out=pixel_size, bounding_box_mode='intersection',
        vectorize_op=False)


def calculate_farm_abundance(species_abundance, ag_mask, alpha, uri):
    """Calculate the farm abundance raster.  The final farm abundance raster
    will be saved to uri.

        species_abundance - a URI to a GDAL dataset of species abundance.
        ag_mask - a uri to a GDAL dataset of values where ag pixels are 1
            and non-ag pixels are 0.
        alpha - the typical foraging distance of the current pollinator.
        uri - the output URI for the farm_abundance raster.

        Returns nothing."""

    farm_abundance_temp_uri = pygeoprocessing.geoprocessing.temporary_filename()
    species_abundance_uri = species_abundance
    species_abundance = gdal.Open(species_abundance_uri)

    pixel_size = abs(species_abundance.GetGeoTransform()[1])
    expected_distance = alpha / pixel_size

    kernel_uri = pygeoprocessing.geoprocessing.temporary_filename()
    make_exponential_decay_kernel_uri(expected_distance, kernel_uri)

    pygeoprocessing.geoprocessing.convolve_2d_uri(
        species_abundance_uri, kernel_uri, farm_abundance_temp_uri)
    os.remove(kernel_uri)

    nodata = species_abundance.GetRasterBand(1).GetNoDataValue()

    # Mask the farm abundance raster according to whether the pixel is
    # agricultural.  If the pixel is agricultural, the value is preserved.
    # Otherwise, the value is set to nodata.
    pygeoprocessing.geoprocessing.vectorize_datasets(
        dataset_uri_list=[farm_abundance_temp_uri, ag_mask],
        dataset_pixel_op=lambda x, y: numpy.where(y == 1.0, x, nodata),
        dataset_out_uri=uri,
        datatype_out=gdal.GDT_Float32,
        nodata_out=nodata,
        pixel_size_out=pixel_size,
        bounding_box_mode='intersection',
        vectorize_op=False)

def create_ag_mask(landuse_uri, out_uri, ag_classes):
    """Reclassify the landuse raster into a raster demarcating the agricultural
    state of a given pixel.  The reclassed ag raster will be saved to uri.

    Parameters:

        landuse_uri (string): a path to the landcover raster.
        out_uri (string): a path to the ouput mask raster
        ag_classes (string): an integer list of landuse classes that could be
            found in `landuse_uri` that are agricultural.

    Returns:
        None"""

    lulc_nodata = pygeoprocessing.get_nodata_from_uri(landuse_uri)
    ag_mask_nodata = 2
    pixel_size_out = pygeoprocessing.get_cell_size_from_uri(landuse_uri)
    def mask_ag_op(lulc_array):
        """Masks lulc array by ag class and preserves nodata"""
        ag_mask = numpy.in1d(lulc_array.flatten(), ag_classes).reshape(
            lulc_array.shape)
        result = numpy.zeros(lulc_array.shape, dtype=numpy.int8)
        result[ag_mask] = 1
        result[lulc_array == lulc_nodata] = ag_mask_nodata
        return result

    pygeoprocessing.vectorize_datasets(
        [landuse_uri], mask_ag_op, out_uri, gdal.GDT_Byte, ag_mask_nodata,
        pixel_size_out, "intersection", vectorize_op=False)


def add_two_rasters(raster_1, raster_2, out_uri):
    """Add two rasters where pixels in raster_1 are not nodata.  Pixels are
        considered to have a nodata value iff the pixel value in raster_1 is
        nodata.  Raster_2's pixel value is not checked for nodata.

        raster_1 - a uri to a GDAL dataset
        raster_2 - a uri to a GDAL dataset
        out_uri - the uri at which to save the resulting raster.

        Returns nothing."""

    # If the user wants us to write the output raster to the URI of one of the
    # input files, create a temporary directory and save the output file to the
    # temp folder.
    temp_dir = None
    if out_uri in [raster_1, raster_2]:
        old_out_uri = out_uri
        temp_dir = True
        out_uri = pygeoprocessing.geoprocessing.temporary_filename()

    nodata = pygeoprocessing.geoprocessing.get_nodata_from_uri(raster_1)
    min_pixel_size = min([
        pygeoprocessing.geoprocessing.get_cell_size_from_uri(x) for x in
        [raster_1, raster_2]])

    pygeoprocessing.geoprocessing.vectorize_datasets(
        dataset_uri_list=[raster_1, raster_2],
        dataset_pixel_op=lambda x, y: numpy.where(
            y != nodata, numpy.add(x, y), nodata),
        dataset_out_uri=out_uri,
        datatype_out=gdal.GDT_Float32,
        nodata_out=nodata,
        pixel_size_out=min_pixel_size,
        bounding_box_mode='intersection',
        vectorize_op=False)

    # If we saved the output file to a temp folder, remove the file that we're
    # trying to avoid and save the temp file to the old file's location.
    if temp_dir != None:
        os.remove(old_out_uri)
        shutil.move(out_uri, old_out_uri)


def calculate_service(rasters, nodata, alpha, part_wild, out_uris):
    """Calculate the service raster.  The finished raster will be saved to
    out_uris['service_value'].

        rasters - a dictionary with these entries:
            'farm_value' - a GDAL dataset.
            'farm_abundance' - a GDAL dataset.
            'species_abundance' - a GDAL dataset.
            'ag_mask' - a GDAL dataset.  Values are either nodata, 0 (if not an
                    ag pixel) or 1 (if an ag pixel).
        nodata - the nodata value for output rasters.
        alpha - the expected distance
        part_wild - a number between 0 and 1 representing the proportion of all
            pollination that is done by wild pollinators.
        out_uris - a dictionary with these entries:
            'species_value' - a URI.  The raster created at this URI will
                represent the part of the farm's value that is attributed to the
                current species.
            'species_value_blurred' - the raster created at this URI
                will be a copy of the species_value raster that has had a
                exponential convolution filter applied to it.
            'service_value' - a URI.  The raster created at this URI will be the
                calculated service value raster.
            'temp' - a folder in which to store temp files.

        Returns nothing."""

    # Open the species foraging matrix and then divide
    # the yield matrix by the foraging matrix for this pollinator.
    min_pixel_size = min([
        pygeoprocessing.geoprocessing.get_cell_size_from_uri(x) for x in
        [rasters['farm_value'], rasters['farm_abundance']]])

    pygeoprocessing.geoprocessing.vectorize_datasets(
        dataset_uri_list=[rasters['farm_value'], rasters['farm_abundance']],
        dataset_pixel_op=lambda x, y: numpy.where(
            x != nodata, x / y, nodata),
        dataset_out_uri=out_uris['species_value'],
        datatype_out=gdal.GDT_Float32,
        nodata_out=nodata,
        pixel_size_out=min_pixel_size,
        bounding_box_mode='intersection',
        vectorize_op=False)

    expected_distance = alpha / min_pixel_size
    kernel_uri = pygeoprocessing.geoprocessing.temporary_filename()
    make_exponential_decay_kernel_uri(expected_distance, kernel_uri)
    pygeoprocessing.geoprocessing.convolve_2d_uri(
        out_uris['species_value'], kernel_uri, out_uris['species_value_blurred'])
    os.remove(kernel_uri)

    # Vectorize the ps_vectorized function
    temp_service_uri = pygeoprocessing.geoprocessing.temporary_filename()

    pygeoprocessing.geoprocessing.vectorize_datasets(
        [rasters['species_abundance'], out_uris['species_value_blurred']],
        lambda x, y: numpy.where(
            x != nodata, part_wild * x * y, nodata),
        dataset_out_uri=temp_service_uri,
        datatype_out=gdal.GDT_Float32,
        nodata_out=nodata,
        pixel_size_out=min_pixel_size,
        bounding_box_mode='intersection',
        vectorize_op=False)

    # Set all agricultural pixels to 0.  This is according to issue 761.
    pygeoprocessing.geoprocessing.vectorize_datasets(
        dataset_uri_list=[rasters['ag_mask'], temp_service_uri],
        dataset_pixel_op=lambda x, y: numpy.where(x == 0, 0.0, y),
        dataset_out_uri=out_uris['service_value'],
        datatype_out=gdal.GDT_Float32,
        nodata_out=nodata,
        pixel_size_out=min_pixel_size,
        bounding_box_mode='intersection',
        vectorize_op=False)


def calculate_yield(in_raster, out_uri, half_sat, wild_poll, out_nodata):
    """Calculate the yield raster.

        in_raster - a uri to a GDAL dataset
        out_uri -a uri for the output (yield) dataset
        half_sat - the half-saturation constant, a python int or float
        wild_poll - the proportion of crops that are pollinated by wild
            pollinators.  An int or float from 0 to 1.
        out_nodata - the nodata value for the output raster

        Returns nothing"""

    # Calculate the yield raster
    kappa_c = float(half_sat)
    nu_c = float(wild_poll)
    nu_c_invert = 1.0 - nu_c
    in_nodata = pygeoprocessing.geoprocessing.get_nodata_from_uri(in_raster)

    # This function is a vectorize-compatible implementation of the yield
    # function from the documentation.
    def calc_yield(frm_avg):
        """Calculate the yield for a farm pixel.  frm_avg is the average
        foraging score on the landscape on this pixel aross all pollinators.
        This function applies the 'expected yield' function from the
        documentation."""
        return numpy.where(
            frm_avg == in_nodata, out_nodata,
            nu_c_invert + (nu_c * frm_avg / (frm_avg + kappa_c)))

    # Apply the yield calculation to the foraging_average raster
    pygeoprocessing.geoprocessing.vectorize_datasets(
        dataset_uri_list=[in_raster],
        dataset_pixel_op=calc_yield,
        dataset_out_uri=out_uri,
        datatype_out=gdal.GDT_Float32,
        nodata_out=out_nodata,
        pixel_size_out=pygeoprocessing.geoprocessing.get_cell_size_from_uri(in_raster),
        bounding_box_mode='intersection',
        vectorize_op=False)


def divide_raster(raster, divisor, uri):
    """Divide all non-nodata values in raster_1 by divisor and save the output
        raster to uri.

        raster - a uri to a GDAL dataset
        divisor - the divisor (a python scalar)
        uri - the uri to which to save the output raster.

        Returns nothing."""

    temp_dir = None
    if raster == uri:
        old_out_uri = uri
        temp_dir = True
        uri = pygeoprocessing.geoprocessing.temporary_filename()

    nodata = pygeoprocessing.geoprocessing.get_nodata_from_uri(raster)
    pygeoprocessing.geoprocessing.vectorize_datasets(
        dataset_uri_list=[raster],
        dataset_pixel_op=lambda x: numpy.where(
            x == nodata, nodata, x / divisor),
        dataset_out_uri=uri,
        datatype_out=gdal.GDT_Float32,
        nodata_out=nodata,
        pixel_size_out=pygeoprocessing.geoprocessing.get_cell_size_from_uri(
            raster),
        bounding_box_mode='intersection',
        vectorize_op=False)

    raster = None
    if temp_dir != None:
        os.remove(old_out_uri)
        shutil.move(uri, old_out_uri)


def map_attribute(base_raster, lu_table_dict, guild_dict, resource_fields,
                  out_uri, list_op):
    """Make an intermediate raster where values are mapped from the base raster
        according to the mapping specified by key_field and value_field.

        base_raster - a URI to a GDAL dataset
        lu_table_dict - a subclass of fileio.AbstractTableHandler
        guild_dict - a python dictionary representing the guild row for this
            species.
        resource_fields - a python list of string resource fields
        out_uri - a uri for the output dataset
        list_op - a python callable that takes a list of numerical arguments
            and returns a python scalar.  Examples: sum; max

        returns nothing."""

    # Get the input raster's nodata value
    base_nodata = pygeoprocessing.geoprocessing.get_nodata_from_uri(base_raster)
    LOGGER.debug(resource_fields)
    LOGGER.debug(guild_dict)
    value_list = dict((r, guild_dict[r]) for r in resource_fields)

    reclass_rules = {base_nodata: -1}
    for lulc in lu_table_dict.keys():
        resource_values = [
            value_list[r] * lu_table_dict[lulc][r] for r in resource_fields]
        reclass_rules[lulc] = list_op(resource_values)

    pygeoprocessing.geoprocessing.reclassify_dataset_uri(
        base_raster, reclass_rules, out_uri, gdal.GDT_Float32, -1)


def make_exponential_decay_kernel_uri(expected_distance, kernel_uri):
    """Create a GDAL raster whose pixel values correspond to an exponentially
    weighted decay distance of the form exp(-distance / expected_distance)

    Paramters:
        expected_distance (int): expected distance to decay to 0 in pixels
        kernel_uri (string): path to a raster that will contain the gaussian
            decay kernel of size (expected_distance * 5  + 1)^2

    Returns:
        None
    """
    max_distance = expected_distance * 5
    kernel_size = int(numpy.round(max_distance * 2 + 1))

    driver = gdal.GetDriverByName('GTiff')
    kernel_dataset = driver.Create(
        kernel_uri.encode('utf-8'), kernel_size, kernel_size, 1, gdal.GDT_Float32,
        options=['BIGTIFF=IF_SAFER'])

    #Make some kind of geotransform, it doesn't matter what but
    #will make GIS libraries behave better if it's all defined
    kernel_dataset.SetGeoTransform([444720, 30, 0, 3751320, 0, -30])
    srs = osr.SpatialReference()
    srs.SetUTM(11, 1)
    srs.SetWellKnownGeogCS('NAD27')
    kernel_dataset.SetProjection(srs.ExportToWkt())

    kernel_band = kernel_dataset.GetRasterBand(1)
    kernel_band.SetNoDataValue(-9999)

    col_index = numpy.array(xrange(kernel_size))
    integration = 0.0
    for row_index in xrange(kernel_size):
        distance_kernel_row = numpy.sqrt(
            (row_index - max_distance) ** 2 + (
                col_index - max_distance) ** 2).reshape(1, kernel_size)
        kernel = numpy.where(
            distance_kernel_row > max_distance, 0.0, numpy.exp(
                -distance_kernel_row / expected_distance))
        integration += numpy.sum(kernel)
        kernel_band.WriteArray(kernel, xoff=0, yoff=row_index)

    for row_index in xrange(kernel_size):
        kernel_row = kernel_band.ReadAsArray(
            xoff=0, yoff=row_index, win_xsize=kernel_size, win_ysize=1)
        kernel_row /= integration
        kernel_band.WriteArray(kernel_row, 0, row_index)
