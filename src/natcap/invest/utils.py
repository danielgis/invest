"""InVEST specific code utils."""
import math
import os

import numpy
from osgeo import gdal
from osgeo import osr
from osgeo import ogr
import pygeoprocessing


def is_gdal_type(arg):
    """Test if input argument is a path to a gdal raster."""
    if (isinstance(arg, str) or
            isinstance(arg, unicode)) and os.path.exists(arg):
        raster = gdal.Open(arg)
        if raster is not None:
            return True
    return False


def is_ogr_type(arg):
    """Test if input argument is a path to an ogr vector."""
    if (isinstance(arg, str) or
            isinstance(arg, unicode)) and os.path.exists(arg):
        vector = ogr.Open(arg)
        if vector is not None:
            return True
    return False


def reproject_bounding_box(bounding_box, source_proj, out_proj):
    """Reproject bounding box to destination source reference system.

    Reprojects the bounding box to account for the the possibility of warping
    due to different projections.

    Parameters:
        bounding_box (list): format of
            [upper_left_x, upper_left_y, lower_right_x, lower_right_y]
            This convention comes from pygeoprocessing v0.3.0 and is expected
            to change in later versions.
        source_proj (osr.SpatialReference): spatial reference for bounding_box
            coordinates
        out_proj (osr.SpatialReference): spatial reference for desired
            transformed bounding_box coordinates

    Returns:
        Transformed `bounding_box` coordinates from `source_proj` to `out_proj`
    """
    # make points for each corner
    bounding_box_points = []
    bounding_box_points.append((bounding_box[0], bounding_box[1]))
    bounding_box_points.append((bounding_box[2], bounding_box[1]))
    bounding_box_points.append((bounding_box[2], bounding_box[3]))
    bounding_box_points.append((bounding_box[0], bounding_box[3]))

    # calculate intermediate points to account for warping when reprojecting
    bounding_box_points.append(
        ((bounding_box[0]+bounding_box[2])/2.0, bounding_box[1]))
    bounding_box_points.append(
        (bounding_box[2], (bounding_box[1]+bounding_box[3])/2.0))
    bounding_box_points.append(
        ((bounding_box[0]+bounding_box[2])/2.0, bounding_box[3]))
    bounding_box_points.append(
        (bounding_box[0], (bounding_box[1]+bounding_box[3])/2.0))

    transformer = osr.CoordinateTransformation(
        source_proj, out_proj)

    transformed_points = []
    for point in bounding_box_points:
        x_coord, y_coord, _ = transformer.TransformPoint(point[0], point[1])
        transformed_points.append((x_coord, y_coord))

    # find the biggest bounding box around the points, initialize to the
    # first point
    out_bounding_box = [
        transformed_points[0][0],
        transformed_points[0][1],
        transformed_points[0][0],
        transformed_points[0][1],
        ]

    union_functions = [min, max, max, min]
    for point in transformed_points:
        for i in range(2):
            # x compare
            out_bounding_box[i*2] = union_functions[i*2](
                point[0], out_bounding_box[i*2])
            # y compare
            out_bounding_box[i*2+1] = union_functions[i*2+1](
                point[1], out_bounding_box[i*2+1])
    return out_bounding_box


def calculate_args_bounding_box(args_dict):
    """Calculate the bounding boxes of paths to GIS data in `args_dict`.

    Parameters:
        args_dict (dict): a string key and any value pair dictionary. Where
            some of the values could be paths to GIS types on disk.

    Returns:
        bb_intersection, bb_union tuple that's either the lat/lng bounding
            intersection and union bounding boxes of the gis types referred to
            in args_dict.  If no GIS types are present, this is a (None, None)
            tuple.
    """

    def _merge_bounding_boxes(bb1, bb2, mode):
        """Merge two bounding boxes through union or intersection.

            Parameters:
                bb1 (list of float): bounding box of the form
                    [minx, maxy, maxx, miny] or None
                bb2 (list of float): bounding box of the form
                    [minx, maxy, maxx, miny] or None
                mode (string): either "union" or "intersection" indicating the
                    how to combine the two bounding boxes.

            Returns:
                either the intersection or union of bb1 and bb2 depending
                on mode.  If either bb1 or bb2 is None, the other is returned.
                If both are None, None is returned.
            """
        if bb1 is None:
            return bb2
        if bb2 is None:
            return bb1

        if mode == "union":
            comparison_ops = [min, max, max, min]
        if mode == "intersection":
            comparison_ops = [max, min, min, max]

        bb_out = [op(x, y) for op, x, y in zip(comparison_ops, bb1, bb2)]
        return bb_out

    def _merge_local_bounding_boxes(arg, bb_intersection=None, bb_union=None):
        """Allows us to recursively walk a potentially nested dictionary
        and merge the bounding boxes that might be found in the GIS
        types

        Args:
            arg (dict): contains string keys and pairs that might be files to
                gis types.  They can be any other type, including dicts.
            bb_intersection (list or None): if list, has the form
                [upper_left_x, upper_left_y, lower_right_x, lower_right_y],
                where coordinates are in lng, lat
            bb_union (list or None): if list, has the form
                [upper_left_x, upper_left_y, lower_right_x, lower_right_y],
                where coordinates are in lng, lat

        Returns:
            (intersection, union) bounding box tuples of all filepaths to GIS
            data types found in the dictionary and bb_intersection and
            bb_union inputs.  None, None if no arguments were GIS data types
            and input bounding boxes are None.
        """
        if isinstance(arg, dict):
            # if dict, grab the bb's for all the members in it
            for value in arg.itervalues():
                bb_intersection, bb_union = _merge_local_bounding_boxes(
                    value, bb_intersection, bb_union)
        elif isinstance(arg, list):
            # if list, grab the bb's for all the members in it
            for value in arg:
                bb_intersection, bb_union = _merge_local_bounding_boxes(
                    value, bb_intersection, bb_union)
        else:
            # singular value, test if GIS type, if not, don't update bb's
            # this is an undefined bounding box that gets returned when ogr
            # opens a table only
            local_bb = [0., 0., 0., 0.]
            if is_gdal_type(arg):
                local_bb = pygeoprocessing.get_bounding_box(arg)
                projection_wkt = (
                    pygeoprocessing.get_dataset_projection_wkt_uri(arg))
                spatial_ref = osr.SpatialReference()
                spatial_ref.ImportFromWkt(projection_wkt)
            elif is_ogr_type(arg):
                local_bb = pygeoprocessing.get_datasource_bounding_box(arg)
                spatial_ref = pygeoprocessing.get_spatial_ref_uri(arg)

            try:
                # means there's a GIS type with a well defined bounding box
                # create transform, and reproject local bounding box to lat/lng
                lat_lng_ref = osr.SpatialReference()
                lat_lng_ref.ImportFromEPSG(4326)  # EPSG 4326 is lat/lng
                to_lat_trans = osr.CoordinateTransformation(
                    spatial_ref, lat_lng_ref)
                for point_index in [0, 2]:
                    local_bb[point_index], local_bb[point_index + 1], _ = (
                        to_lat_trans.TransformPoint(
                            local_bb[point_index], local_bb[point_index + 1]))

                bb_intersection = _merge_bounding_boxes(
                    local_bb, bb_intersection, 'intersection')
                bb_union = _merge_bounding_boxes(
                    local_bb, bb_union, 'union')
            except Exception:
                # All kinds of exceptions from bad transforms or CSV files
                # or dbf files could get us to this point, just don't bother
                # with the local_bb at all
                pass

        return bb_intersection, bb_union

    return _merge_local_bounding_boxes(args_dict)


def make_suffix_string(args, suffix_key):
    """Make an InVEST appropriate suffix string.

    Creates an InVEST appropriate suffix string  given the args dictionary and
    suffix key.  In general, prepends an '_' when necessary and generates an
    empty string when necessary.

    Parameters:
        args (dict): the classic InVEST model parameter dictionary that is
            passed to `execute`.
        suffix_key (string): the key used to index the base suffix.

    Returns:
        If `suffix_key` is not in `args`, or `args['suffix_key']` is ""
            return "",
        If `args['suffix_key']` starts with '_' return `args['suffix_key']`
            else return '_'+`args['suffix_key']`
    """
    try:
        file_suffix = args[suffix_key]
        if file_suffix != "" and not file_suffix.startswith('_'):
            file_suffix = '_' + file_suffix
    except KeyError:
        file_suffix = ''

    return file_suffix


def exponential_decay_kernel_raster(expected_distance, kernel_filepath):
    """Create a raster-based exponential decay kernel.

    The raster created will be a tiled GeoTiff, with 256x256 memory blocks.

    Parameters:
        expected_distance (int or float): The distance (in pixels) of the
            kernel's radius, the distance at which the value of the decay
            function is equal to `1/e`.
        kernel_filepath (string): The path to the file on disk where this
            kernel should be stored.  If this file exists, it will be
            overwritten.

    Returns:
        None
    """
    max_distance = expected_distance * 5
    kernel_size = int(numpy.round(max_distance * 2 + 1))

    driver = gdal.GetDriverByName('GTiff')
    kernel_dataset = driver.Create(
        kernel_filepath.encode('utf-8'), kernel_size, kernel_size, 1,
        gdal.GDT_Float32, options=[
            'BIGTIFF=IF_SAFER', 'TILED=YES', 'BLOCKXSIZE=256',
            'BLOCKYSIZE=256'])

    # Make some kind of geotransform, it doesn't matter what but
    # will make GIS libraries behave better if it's all defined
    kernel_dataset.SetGeoTransform([444720, 30, 0, 3751320, 0, -30])
    srs = osr.SpatialReference()
    srs.SetUTM(11, 1)
    srs.SetWellKnownGeogCS('NAD27')
    kernel_dataset.SetProjection(srs.ExportToWkt())

    kernel_band = kernel_dataset.GetRasterBand(1)
    kernel_band.SetNoDataValue(-9999)

    cols_per_block, rows_per_block = kernel_band.GetBlockSize()

    n_cols = kernel_dataset.RasterXSize
    n_rows = kernel_dataset.RasterYSize

    n_col_blocks = int(math.ceil(n_cols / float(cols_per_block)))
    n_row_blocks = int(math.ceil(n_rows / float(rows_per_block)))

    integration = 0.0
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

            # Numpy creates index rasters as ints by default, which sometimes
            # creates problems on 32-bit builds when we try to add Int32
            # matrices to float64 matrices.
            row_indices, col_indices = numpy.indices((row_block_width,
                                                      col_block_width),
                                                     dtype=numpy.float)

            row_indices += numpy.float(row_offset - max_distance)
            col_indices += numpy.float(col_offset - max_distance)

            kernel_index_distances = numpy.hypot(
                row_indices, col_indices)
            kernel = numpy.where(
                kernel_index_distances > max_distance, 0.0,
                numpy.exp(-kernel_index_distances / expected_distance))
            integration += numpy.sum(kernel)

            kernel_band.WriteArray(kernel, xoff=col_offset,
                                   yoff=row_offset)

    # Need to flush the kernel's cache to disk before opening up a new Dataset
    # object in interblocks()
    kernel_dataset.FlushCache()

    for block_data, kernel_block in pygeoprocessing.iterblocks(
            kernel_filepath):
        kernel_block /= integration
        kernel_band.WriteArray(kernel_block, xoff=block_data['xoff'],
                               yoff=block_data['yoff'])


def build_file_registry(base_file_path_list, file_suffix):
    """Combine file suffixes with key names, base filenames, and directories.

    Parameters:
        base_file_tuple_list (list): a list of (dict, path) tuples where
            the dictionaries have a 'file_key': 'basefilename' pair, or
            'file_key': list of 'basefilename's.  'path'
            indicates the file directory path to prepend to the basefile name.
        file_suffix (string): a string to append to every filename, can be
            empty string

    Returns:
        dictionary of 'file_keys' from the dictionaries in
        `base_file_tuple_list` mapping to full file paths with suffixes or
        lists of file paths with suffixes depending on the original type of
        the 'basefilename' pair.

    Raises:
        ValueError if there are duplicate file keys or duplicate file paths.
    """
    all_paths = set()
    duplicate_keys = set()
    duplicate_paths = set()
    f_reg = {}

    def _build_path(base_filename, path):
        """Internal helper to avoid code duplication."""
        pre, post = os.path.splitext(base_filename)
        full_path = os.path.join(path, pre+file_suffix+post)

        # Check for duplicate keys or paths
        if full_path in all_paths:
            duplicate_paths.add(full_path)
        else:
            all_paths.add(full_path)
        return full_path

    for base_file_dict, path in base_file_path_list:
        for file_key, file_payload in base_file_dict.iteritems():
            # check for duplicate keys
            if file_key in f_reg:
                duplicate_keys.add(file_key)
            else:
                # handle the case whether it's a filename or a list of strings
                if isinstance(file_payload, basestring):
                    full_path = _build_path(file_payload, path)
                    f_reg[file_key] = full_path
                elif isinstance(file_payload, list):
                    f_reg[file_key] = []
                    for filename in file_payload:
                        full_path = _build_path(filename, path)
                        f_reg[file_key].append(full_path)

    if len(duplicate_paths) > 0 or len(duplicate_keys):
        raise ValueError(
            "Cannot consolidate because of duplicate paths or keys: "
            "duplicate_keys: %s duplicate_paths: %s" % (
                duplicate_keys, duplicate_paths))

    return f_reg
