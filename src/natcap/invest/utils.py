"""File and GIS utilities that are used in more than one place across the
project but don't fall into a formal GIS or ecosystem service classification"""

import os

from osgeo import gdal
from osgeo import ogr
from osgeo import osr
import pygeoprocessing


def is_gdal_type(arg):
    """tests if input argument is a path to a gdal raster"""
    if (isinstance(arg, str) or
            isinstance(arg, unicode)) and os.path.exists(arg):
        raster = gdal.Open(arg)
        if raster is not None:
            return True
    return False


def is_ogr_type(arg):
    """tests if input argument is a path to an ogr vector"""
    if (isinstance(arg, str) or
            isinstance(arg, unicode)) and os.path.exists(arg):
        vector = ogr.Open(arg)
        if vector is not None:
            return True
    return False


def reproject_bounding_box(bounding_box, source_proj, out_proj):
    """Reproject bounding box to destination source reference system accounting
    for the possibility of warping due to different projections.

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

    print transformed_points
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
    """Parse through an InVEST style args dict and calculate the bounding boxes
    of any paths that can be interpreted as GIS types.

    Parameters:
        args_dict (dict): a string key and any value pair dictionary. Where
            some of the values could be paths to GIS types on disk.

    Returns:
        bb_intersection, bb_union tuple that's either the lat/lng bounding
            intersection and union bounding boxes of the gis types referred to
            in args_dict.  If no GIS types are present, this is a (None, None)
            tuple."""

    def _merge_bounding_boxes(bb1, bb2, mode):
        """Helper function to merge two bounding boxes through union or
            intersection

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
                gis types.  They can be any other type, including dictionaries.
            bb_intersection (list or None): if list, has the form
                [upper_left_x, upper_left_y, lower_right_x, lower_right_y], where coordinates are in lng, lat
            bb_union (list or None): if list, has the form
                [upper_left_x, upper_left_y, lower_right_x, lower_right_y], where coordinates are in lng, lat

        Returns:
            (intersection, union) bounding box tuples of all filepaths to GIS
            data types found in the dictionary and bb_intersection and bb_union
            inputs.  None, None if no arguments were GIS data types and input
            bounding boxes are None."""

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
                projection_wkt = pygeoprocessing.get_dataset_projection_wkt_uri(
                    arg)
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
