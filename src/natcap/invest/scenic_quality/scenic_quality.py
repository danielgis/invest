import os
import sys
import math
import shutil

import numpy as np
import scipy.stats
import logging
import cProfile
import pstats
import h5py

from osgeo import gdal
from osgeo import ogr
from osgeo import osr

import pygeoprocessing.geoprocessing
from natcap.invest.scenic_quality import scenic_quality_core
import scenic_quality_cython_core

logging.basicConfig(format='%(asctime)s %(name)-20s %(levelname)-8s \
%(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

LOGGER = logging.getLogger('scenic_quality')

def old_reproject_dataset_uri(original_dataset_uri, *args, **kwargs):
    """A URI wrapper for reproject dataset that opens the original_dataset_uri
        before passing it to reproject_dataset.

       original_dataset_uri - a URI to a gdal Dataset on disk

       All other arguments to reproject_dataset are passed in.

       return - nothing"""

    original_dataset = gdal.Open(original_dataset_uri)
    reproject_dataset(original_dataset, *args, **kwargs)

    pygeoprocessing.geoprocessing.calculate_raster_stats_uri(original_dataset_uri)

def reproject_dataset_uri(original_dataset_uri, output_wkt, output_uri,
                      output_type = gdal.GDT_Float32):
    """A function to reproject and resample a GDAL dataset given an output pixel size
        and output reference and uri.

       original_dataset - a gdal Dataset to reproject
       pixel_spacing - output dataset pixel size in projected linear units (probably meters)
       output_wkt - output project in Well Known Text (the result of ds.GetProjection())
       output_uri - location on disk to dump the reprojected dataset
       output_type - gdal type of the output    

       return projected dataset"""

    original_dataset = gdal.Open(original_dataset_uri)

    original_sr = osr.SpatialReference()
    original_sr.ImportFromWkt(original_dataset.GetProjection())

    output_sr = osr.SpatialReference()
    output_sr.ImportFromWkt(output_wkt)

    vrt = gdal.AutoCreateWarpedVRT(original_dataset, None, output_wkt, gdal.GRA_Bilinear) 

    # Get the Geotransform vector
    geo_t = vrt.GetGeoTransform()
    x_size = vrt.RasterXSize # Raster xsize
    y_size = vrt.RasterYSize # Raster ysize

    # Work out the boundaries of the new dataset in the target projection
    

    gdal_driver = gdal.GetDriverByName('GTiff')
    # The size of the raster is given the new projection and pixel spacing
    # Using the values we calculated above. Also, setting it to store one band
    # and to use Float32 data type.

    output_dataset = gdal_driver.Create(output_uri, x_size, 
                              y_size, 1, output_type)

    # Set the nodata value
    original_band = original_dataset.GetRasterBand(1)
    out_nodata = original_band.GetNoDataValue()
    original_band.SetNoDataValue(out_nodata)

    # Set the geotransform
    output_dataset.SetGeoTransform(geo_t)
    output_dataset.SetProjection (output_sr.ExportToWkt())

    # Perform the projection/resampling 
    gdal.ReprojectImage(original_dataset, output_dataset,
                        original_sr.ExportToWkt(), output_sr.ExportToWkt(),
                        gdal.GRA_Bilinear)
    
    pygeoprocessing.geoprocessing.calculate_raster_stats_uri(output_uri)


def reclassify_quantile_dataset_uri( \
    dataset_uri, quantile_list, dataset_out_uri, datatype_out, nodata_out):

    nodata_ds = pygeoprocessing.geoprocessing.get_nodata_from_uri(dataset_uri)

    memory_file_uri = pygeoprocessing.geoprocessing.temporary_filename()
    memory_array = pygeoprocessing.geoprocessing.load_memory_mapped_array(dataset_uri, memory_file_uri)
    memory_array_flat = memory_array.reshape((-1,))

    quantile_breaks = [0]
    for quantile in quantile_list:
        quantile_breaks.append(scipy.stats.scoreatpercentile(
                memory_array_flat, quantile, (0.0, np.amax(memory_array_flat))))
        LOGGER.debug('quantile %f: %f', quantile, quantile_breaks[-1])

    def reclass(value):
        if value == nodata_ds:
            return nodata_out
        else:
            for new_value,quantile_break in enumerate(quantile_breaks):
                if value <= quantile_break:
                    return new_value
        raise ValueError, "Value was not within quantiles."

    cell_size = pygeoprocessing.geoprocessing.get_cell_size_from_uri(dataset_uri)

    pygeoprocessing.geoprocessing.vectorize_datasets([dataset_uri],
                                    reclass,
                                    dataset_out_uri,
                                    datatype_out,
                                    nodata_out,
                                    cell_size,
                                    "union",
                                    dataset_to_align_index=0)

    pygeoprocessing.geoprocessing.calculate_raster_stats_uri(dataset_out_uri)

def get_data_type_uri(ds_uri):
    raster_ds = gdal.Open(ds_uri)
    band = raster_ds.GetRasterBand(1)
    raster_data_type = band.DataType
    band = None
    raster_ds = None

    return raster_data_type

def compute_viewshed_uri(in_dem_uri, out_viewshed_uri, in_structure_uri, 
    curvature_correction, refr_coeff, args):
    """ Compute the viewshed as it is defined in ArcGIS where the inputs are:

        -in_dem_uri: URI to input surface raster
        -out_viewshed_uri: URI to the output raster
        -in_structure_uri: URI to a point shapefile that contains the location
        of the observers and the viewshed radius in (negative) meters
        -curvature_correction: flag for the curvature of the earth. Either
        FLAT_EARTH or CURVED_EARTH. Not used yet.
        -refraction: refraction index between 0 (max effect) and 1 (no effect).
        Default is 0.13."""

    # Extract cell size from input DEM
    cell_size = pygeoprocessing.geoprocessing.get_cell_size_from_uri(in_dem_uri)

    # Extract nodata
    nodata = pygeoprocessing.geoprocessing.get_nodata_from_uri(in_dem_uri)
    
    ## Build I and J arrays, and save them to disk
    rows, cols = pygeoprocessing.geoprocessing.get_row_col_from_uri(in_dem_uri)
    I, J = np.meshgrid(range(rows), range(cols), indexing = 'ij')
    # Base path uri
    base_uri = os.path.split(out_viewshed_uri)[0]
    I_uri = os.path.join(base_uri, 'I.tif')
    J_uri = os.path.join(base_uri, 'J.tif')
    pygeoprocessing.geoprocessing.new_raster_from_base_uri(in_dem_uri, I_uri, 'GTiff', \
        -32768., gdal.GDT_Float32, fill_value = -32768.)
    I_raster = gdal.Open(I_uri, gdal.GA_Update)
    I_band = I_raster.GetRasterBand(1)
    I_band.WriteArray(I)
    I_band = None
    I_raster = None
    pygeoprocessing.geoprocessing.new_raster_from_base_uri(in_dem_uri, J_uri, 'GTiff', \
        -32768., gdal.GDT_Float32, fill_value = -32768.)
    J_raster = gdal.Open(J_uri, gdal.GA_Update)
    J_band = J_raster.GetRasterBand(1)
    J_band.WriteArray(J)
    J_band = None
    J_raster = None
    # Extract the input raster geotransform
    GT = pygeoprocessing.geoprocessing.get_geotransform_uri(in_dem_uri)

    # Open the input URI and extract the numpy array
    input_nodata = pygeoprocessing.geoprocessing.get_nodata_from_uri(in_dem_uri)
    input_raster = gdal.Open(in_dem_uri)
    input_array = input_raster.GetRasterBand(1).ReadAsArray()
    input_array[input_array == input_nodata] = 0.
    input_raster = None

    # Create a raster from base before passing it to viewshed
    visibility_uri = out_viewshed_uri #pygeoprocessing.geoprocessing.temporary_filename()
    pygeoprocessing.geoprocessing.new_raster_from_base_uri(in_dem_uri, visibility_uri, 'GTiff', \
        0., gdal.GDT_Float64, fill_value = 0.)

    # Call the non-uri version of viewshed.
#    compute_viewshed(input_array, visibility_uri, in_structure_uri, \
#        cell_size, rows, cols, nodata, GT, I_uri, J_uri, \
#        curvature_correction, refr_coeff, args)

    cProfile.runctx( \
        'compute_viewshed(input_array, visibility_uri, in_structure_uri, \
        cell_size, rows, cols, nodata, GT, I_uri, J_uri, \
        curvature_correction, refr_coeff, args)' \
        , globals(), locals(), 'stats')
    p = pstats.Stats('stats')
    p.sort_stats("time").print_stats(20)
    p.sort_stats('cumulative').print_stats(20)

    os.remove(I_uri)
    os.remove(J_uri)

def compute_viewshed(input_array, visibility_uri, in_structure_uri, \
    cell_size, rows, cols, nodata, GT, I_uri, J_uri, curvature_correction, \
    refr_coeff, args):
    """ array-based function that computes the viewshed as is defined in ArcGIS
    """
    # default parameter values that are not passed to this function but that
    # scenic_quality_core.viewshed needs
    obs_elev = 1.0 # Observator's elevation in meters
    tgt_elev = 0.0  # Extra elevation applied to all the DEM
    max_dist = -1.0 # max. viewing distance(m). Distance is infinite if negative
    coefficient = 1.0 # Used to weight the importance of individual viewsheds
    height = 0.0 # Per viewpoint height offset--updated as we read file info

    # Setup valuation function
    a = args["a_coefficient"]
    b = args["b_coefficient"]
    c = args["c_coefficient"]
    d = args["d_coefficient"]

    valuation_function = None
    max_valuation_radius = args['max_valuation_radius']
    if "polynomial" in args["valuation_function"]:
        print("Polynomial")
        valuation_function = scenic_quality_cython_core.polynomial
    elif "logarithmic" in args['valuation_function']:
        print("logarithmic")
        valuation_function = scenic_quality_cython_core.logarithmic

    assert valuation_function is not None
    
    # Make sure the values don't become too small at max_valuation_radius:
#    test_function = valuation_function(0, 0, max_valuation_radius, 1)
#    edge_value = test_function(np.array([0]), \
#                                    np.array([1]), \
#                                    np.array([1]), \
#                                    np.array([0]))
#    edge_value = valuation_function(max_valuation_radius, 1)
#    message = "Valuation function can't be negative if evaluated at " + \
#    str(max_valuation_radius) + " meters (value is " + str(edge_value) + ")"
#    assert edge_value >= 0., message
        
    # Base path uri
    base_uri = os.path.split(visibility_uri)[0]

    # The model extracts each viewpoint from the shapefile
    point_list = []
    shapefile = ogr.Open(in_structure_uri)
    assert shapefile is not None
    layer = shapefile.GetLayer(0)
    assert layer is not None
    success, iGT = gdal.InvGeoTransform(GT)
    if not success:
        raise RuntimeError('Cannot invert transformation')

    feature_count = layer.GetFeatureCount()

    array_shape = (rows, cols)
    
    # Create arrays for the parameters used for the computation at each viewpoint
    max_distances = np.ones(feature_count).astype(int) * -1
    coefficients = np.ones(feature_count)
    obs_elevations = np.zeros(feature_count)
    tgt_elevations = np.zeros(feature_count)
    viewpoint_row = np.ones(feature_count).astype(int) * -1
    viewpoint_col = np.ones(feature_count).astype(int) * -1

    print('Number of viewpoints: ' + str(feature_count))
    for f in range(feature_count):
        feature = layer.GetFeature(f)
        field_count = feature.GetFieldCount()
        # Check for feature information (radius, coeff, height)
        for field in range(field_count):
            field_def = feature.GetFieldDefnRef(field)
            field_name = field_def.GetNameRef()
            if (field_name.upper() == 'RADIUS2') or \
                (field_name.upper() == 'RADIUS'):
                max_dist = abs(int(feature.GetField(field)))
                assert max_dist is not None, "max distance can't be None"
                if max_dist < args['max_valuation_radius']:
                    LOGGER.warning( \
                        'Valuation radius > maximum visibility distance: ' + \
                        '(' + str(args['max_valuation_radius']) + ' < ' + \
                        str(max_dist) + ')')
                    LOGGER.warning( \
                        'The valuation is performed beyond what is visible')
                max_distances[f] = int(max_dist/cell_size)
            if field_name.lower() == 'coeff':
                coefficients[f] = float(feature.GetField(field))
                assert coefficient is not None, "feature coeff can't be None"
            if field_name.lower() == 'offseta':
                obs_elevations[f] = float(feature.GetField(field))
                assert obs_elev is not None, "OFFSETA can't be None"
            if field_name.lower() == 'offsetb':
                tgt_elevations[f] = float(feature.GetField(field))
                assert tgt_elev is not None, "OFFSETB can't be None"
                
        geometry = feature.GetGeometryRef()
        assert geometry is not None
        message = 'geometry type is ' + str(geometry.GetGeometryName()) + \
        ' point is "POINT"'
        assert geometry.GetGeometryName() == 'POINT', message
        x = geometry.GetX()
        y = geometry.GetY()
        viewpoint_col[f] = int((iGT[0] + x*iGT[1] + y*iGT[2]))
        viewpoint_row[f] = int((iGT[3] + x*iGT[4] + y*iGT[5]))

        v_col, v_row = gdal.ApplyGeoTransform(iGT, x, y)

        assert abs(int(v_row) -viewpoint_row[f]) < 1e-8, \
            "coords don't agree: " + str([v_row, v_col]) + ' vs ' + \
            str([viewpoint_row[f], viewpoint_col[f]])

        assert abs(int(v_col) -viewpoint_col[f]) < 1e-8, \
            "coords don't agree: " + str([v_row, v_col]) + ' vs ' + \
            str([viewpoint_row[f], viewpoint_col[f]])

        assert viewpoint_row[f] >= 0, \
            'Viewpoint ' +str(f) +"'s row is negative:" + str(viewpoint_row[f])
        assert viewpoint_row[f] < input_array.shape[0], \
            'Viewpoint ' +str(f) +"'s row is too large:" + \
            str(viewpoint_row[f]) + ' raster size is ' + str(input_array.shape)
        assert viewpoint_col[f] >= 0, \
            'Viewpoint ' +str(f) +"'s col is negative:" + str(viewpoint_col[f])
        assert viewpoint_col[f] < input_array.shape[1], \
            'Viewpoint ' +str(f) +"'s col is too large:" + \
            str(viewpoint_col[f]) + ' raster size is ' + str(input_array.shape)


    layer = None
    shapefile = None

    arg_dist = np.argsort(max_distances)


    I_raster = gdal.Open(I_uri, gdal.GA_Update)
    I_band = I_raster.GetRasterBand(1)
    I_array = I_band.ReadAsArray()
    I_band = None
    I_raster = None
    J_raster = gdal.Open(J_uri, gdal.GA_Update)
    J_band = J_raster.GetRasterBand(1)
    J_array = J_band.ReadAsArray()
    J_band = None
    J_raster = None

    visibility_raster = gdal.Open(visibility_uri)
    visibility_band = visibility_raster.GetRasterBand(1)
    accum_visibility = visibility_band.ReadAsArray()
    accum_visibility[:] = 0.0
    visibility_band = None
    visibility_raster = None

    # Compute distances
    distances_array = np.copy(accum_visibility)

    last_dist = 0
    for dist in range(arg_dist.size-1, -1, -1):
        f = arg_dist[dist]
        print("Iteration " + str(dist) + ", processing viewpoint " + str(f))

        max_dist = max_distances[f]
        coefficient = coefficients[f]
        obs_elev = obs_elevations[f]
        tgt_elev = tgt_elevations[f]
        i = viewpoint_row[f]
        j = viewpoint_col[f]

        # Create a visibility map
        # Visibility convention: 1 visible, \
        # <0 is additional height to become visible
        visibility_map = np.zeros(input_array.shape)
        visibility_map[input_array == nodata] = 2. 

        # perimeter
        # 1- get perimeter cells
        perimeter_cells = \
            scenic_quality_core.get_perimeter_cells(input_array.shape, (i,j), max_dist)

        # angles
        # cell_angles + append the last element (2 PI) automatically
        angles = scenic_quality_core.cell_angles(perimeter_cells, (i,j))
        angles = np.append(angles, 2.0 * math.pi)
        
        # Viewshed information
        row_max = np.amax(perimeter_cells[0])
        row_min = np.amin(perimeter_cells[0])
        col_max = np.amax(perimeter_cells[1])    
        col_min = np.amin(perimeter_cells[1])
        
        # Shape of the viewshed
        viewshed_shape = (row_max-row_min + 1, col_max-col_min + 1)

        # Viewer's coordiantes relative to the viewshed 
        v = (i - row_min, j - col_min)

        # add_events, center_events, remove_events
        add_events, center_events, remove_events, I, J = \
            scenic_quality_cython_core.list_extreme_cell_angles(viewshed_shape, \
            v, max_dist)

        # arg_min, arg_center, arg_max
        arg_min = np.argsort(add_events).astype(np.int64)
        arg_max = np.argsort(remove_events).astype(np.int64)
        arg_center = np.argsort(center_events).astype(np.int64)

#        # Debug -- Testing the new memory efficient algorithm
#        scenic_quality_cython_core.memory_efficient_event_stream(viewshed_shape, \
#            v, max_dist)
#
#        return

        # I and J are relative to the viewshed_shape. Make them absolute
        I += row_min
        J += col_min

        # coord
        coord = np.array([I, J])

        # Compute distances
        distances = np.copy(I).astype(np.float64)
        distances_sq = np.copy(I).astype(np.float64)

        scenic_quality_cython_core.compute_distances( \
            i, j, cell_size, I, J, distances_sq, distances)

        # Computation of the visibility:
        # 1- get the height of the DEM w.r.t. the viewer's elevatoin (coord+elev)
        visibility = (input_array[(I, J)] - \
            input_array[i, j] - obs_elev).astype(np.float64)
        offset_visibility = visibility + tgt_elev
    
        # 2- Factor the effect of refraction in the elevation.
        # From the equation on the ArcGIS website:
        # http://resources.arcgis.com/en/help/main/10.1/index.html#//00q90000008v000000
        D_earth = 12740000. # Diameter of the earth in meters
        correction = distances_sq.astype(float) * \
            (refr_coeff - 1.) / D_earth
        visibility += correction
        offset_visibility += correction
        # 3- Divide the height by the distance to get a visibility score
        visibility /= distances
        offset_visibility /= distances


#        # Save viewshed function arguments in HDF5
#        # Inputs to save:
#        #   -add_events, center_events, remove_events
#        #   -arg_min, arg_center, arg_max
#        #   -distances, distances_sq
#        #   -visibility, offset_visibility
#        #
#        # Create the paths to the debug data
#        debug_uri = os.path.join(args['intermediate_dir'], \
#            'debug_data_' + str(f) + '.h5')
#
#        debug_data = h5py.File(debug_uri, 'w')
#
#
#        add_events_dataset = debug_data.create_dataset('add_events', 
#            add_events.shape, compression = 'gzip', fillvalue = -1)
#
#        center_events_dataset = debug_data.create_dataset('center_events', 
#            center_events.shape, compression = 'gzip', fillvalue = -1)
#        
#        remove_events_dataset = debug_data.create_dataset('remove_events', 
#            remove_events.shape, compression = 'gzip', fillvalue = -1)
#        
#
#        arg_min_dataset = debug_data.create_dataset('arg_min', 
#            arg_min.shape, compression = 'gzip', fillvalue = -1)
#
#        arg_center_dataset = debug_data.create_dataset('arg_center', 
#            arg_center.shape, compression = 'gzip', fillvalue = -1)
#
#        arg_max_dataset = debug_data.create_dataset('arg_max', 
#            arg_max.shape, compression = 'gzip', fillvalue = -1)
#
#
#        distances_dataset = debug_data.create_dataset('distances', 
#            distances.shape, compression = 'gzip', fillvalue = -1)
#
#        distances_sq_dataset = debug_data.create_dataset('distances_sq', 
#            distances_sq.shape, compression = 'gzip', fillvalue = -1)
#        
#
#        # Store data in the file
#        add_events_dataset[...] = add_events[...]
#        center_events_dataset[...] = center_events[...]
#        remove_events_dataset[...] = remove_events[...]
#        
#        arg_min_dataset[...] = arg_min[...]
#        arg_center_dataset[...] = arg_center[...]
#        arg_max_dataset[...] = arg_max[...]
#
#        distances_dataset[...] = distances[...]
#        distances_sq_dataset[...] = distances_sq[...]
#
#
#        # Close the files
#        debug_data.close()



        # compute viewshed
        visibility_array = scenic_quality_core.viewshed(
            input_array, cell_size, visibility_map, perimeter_cells, \
            (i,j), angles, v, viewshed_shape, row_min, col_min, \
            add_events, center_events, remove_events, I, J, \
            arg_min, arg_max, arg_center, \
            coord, distances_sq, distances, visibility, offset_visibility, \
            obs_elev, tgt_elev, max_dist, refr_coeff, alg_version='cython', 
            path = args['intermediate_dir'], index = dist)
        
        # apply valuation function
        valuation_function(a, b, c, d, \
            max_valuation_radius, i, j, cell_size, \
            coefficient , I, J, distances, \
            visibility_array, accum_visibility)

        # Write temporary result on disk
        tmp_uri = os.path.join(base_uri, 'viewshed_' + str(dist) + '_' + str(f) + '.tif')
        pygeoprocessing.geoprocessing.new_raster_from_base_uri(visibility_uri, tmp_uri, 'GTiff', \
            0., gdal.GDT_Float32, fill_value = 0)
        tmp_raster = gdal.Open(tmp_uri, gdal.GA_Update)
        tmp_band = tmp_raster.GetRasterBand(1)
        tmp_band.WriteArray(visibility_array)
        tmp_band = None
        tmp_raster = None
        pygeoprocessing.geoprocessing.calculate_raster_stats_uri(tmp_uri)
        visibility_raster = gdal.Open(visibility_uri, gdal.GA_Update)
        visibility_band = visibility_raster.GetRasterBand(1)
        visibility_band.WriteArray(accum_visibility)
        visibility_band = None
        visibility_raster = None
        pygeoprocessing.geoprocessing.calculate_raster_stats_uri(visibility_uri)
        
        last_dist = max_distances[f]

    # Write result on disk
    visibility_raster = gdal.Open(visibility_uri, gdal.GA_Update)
    visibility_band = visibility_raster.GetRasterBand(1)
    visibility_band.WriteArray(accum_visibility)
    visibility_band = None
    visibility_raster = None
    pygeoprocessing.geoprocessing.calculate_raster_stats_uri(visibility_uri)

def add_field_feature_set_uri(fs_uri, field_name, field_type):
    shapefile = ogr.Open(fs_uri, 1)
    layer = shapefile.GetLayer()
    new_field = ogr.FieldDefn(field_name, field_type)
    layer.CreateField(new_field)
    shapefile = None    

def add_id_feature_set_uri(fs_uri, id_name):
    shapefile = ogr.Open(fs_uri, 1)
    message = "Failed to open " + fs_uri + ": can't add new field."
    assert shapefile is not None, message
    layer = shapefile.GetLayer()
    new_field = ogr.FieldDefn(id_name, ogr.OFTInteger)
    layer.CreateField(new_field)

    for feature_id in xrange(layer.GetFeatureCount()):
        feature = layer.GetFeature(feature_id)
        feature.SetField(id_name, feature_id)
        layer.SetFeature(feature)
    shapefile = None

def set_field_by_op_feature_set_uri(fs_uri, value_field_name, op):
    shapefile = ogr.Open(fs_uri, 1)
    layer = shapefile.GetLayer()

    for feature_id in xrange(layer.GetFeatureCount()):
        feature = layer.GetFeature(feature_id)
        feature.SetField(value_field_name, op(feature))
        layer.SetFeature(feature)
    shapefile = None

def get_count_feature_set_uri(fs_uri):
    shapefile = ogr.Open(fs_uri)
    message = "Problem encoutered with " + fs_uri
    assert shapefile is not None, message
    layer = shapefile.GetLayer()
    assert layer is not None, message
    count = layer.GetFeatureCount()
    shapefile = None

    return count

def execute(args):
    """
    Entry point into the Scenic Quality Model

    Args:
        workspace_dir (string): Results will be saved to this location.
        aoi_uri (string): This AOI instructs the model where to clip
            the input data and the extent of analysis.
        cell_size (float): Length (in meters) of each side of the (square)
            cell.
        structure_uri (string): Shapefile indicating the locations of objects
            that contribute to negative scenic quality, such as aquaculture
            netpens or wave energy facilities.
        dem_uri (string): An elevation raster layer is required to conduct
            viewshed analysis.
        refraction (float): The earth curvature correction option corrects
            for the curvature of the earth and refraction of visible light in
            air.
        pop_uri (string): Raster used to determine population within the
            AOI's land-seascape where point features contributing to negative
            scenic quality are visible and not visible.
        overlap_uri (string): A shapefile used to determine the impact of
            objects on visual quality. This input must be a polygon and
            projected in meters.
        valuation_function (string): Indicates the functional form f(x) the
            model will use to value the visual impact for each viewpoint.
        a_coefficient (float): first coefficient for valuation function
        b_coefficient (float): second coefficient for valuation function
        c_coefficient (float): third coefficient for valuation function
        d_coefficient (float): fourth coefficient for valuation function
        max_valuation_radius (float): Radius beyond which the valuation is set
            to zero.

    Example Args Dictionary::

        {
            'workspace_dir': 'path/to/workspace_dir',
            'aoi_uri': 'path/to/shapefile',
            'cell_size': 500,
            'structure_uri': 'path/to/shapefile',
            'dem_uri': 'path/to/raster',
            'refraction': 0.13,
            'pop_uri': 'path/to/shapefile',
            'overlap_uri': 'path/to/shapefile',
            'valuation_function': 'polynomial: a + bx + cx^2 + dx^3',
            'a_coefficient': 1.0,
            'b_coefficient': 0.0,
            'c_coefficient': 0.0,
            'd_coefficient': 0.0,
            'max_valuation_radius': 8000.0,

        }

    """
    LOGGER.info("Start Scenic Quality Model")

    #validate input
    LOGGER.debug("Validating parameters.")
    dem_cell_size=pygeoprocessing.geoprocessing.get_cell_size_from_uri(args['dem_uri'])
    LOGGER.debug("DEM cell size: %f" % dem_cell_size)
    if "cell_size" in args:
        if args['cell_size'] < dem_cell_size:
            raise ValueError, "The cell size cannot be downsampled below %f" % dem_cell_size
    else:
        args['cell_size'] = dem_cell_size

    intermediate_dir = os.path.join(args['workspace_dir'], 'intermediate')
    if not os.path.isdir(intermediate_dir):
        os.makedirs(intermediate_dir)
    args['intermediate_dir'] = intermediate_dir

    output_dir = os.path.join(args['workspace_dir'], 'output')
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    args['output_dir'] = output_dir

    #local variables
    LOGGER.debug("Setting local variables.")
    z_factor=1
    curvature_correction=args['refraction']

    #intermediate files
    aoi_dem_uri=os.path.join(intermediate_dir,"aoi_dem.shp")
    aoi_pop_uri=os.path.join(intermediate_dir,"aoi_pop.shp")

    viewshed_dem_uri=os.path.join(intermediate_dir,"dem_vs.tif")
    viewshed_dem_reclass_uri=os.path.join(intermediate_dir,"dem_vs_re.tif")

    pop_clip_uri=os.path.join(intermediate_dir,"pop_clip.tif")
    pop_prj_uri=os.path.join(intermediate_dir,"pop_prj.tif")
    pop_vs_uri=os.path.join(intermediate_dir,"pop_vs.tif")

    viewshed_reclass_uri=os.path.join(intermediate_dir,"vshed_bool.tif")
    viewshed_polygon_uri=os.path.join(intermediate_dir,"vshed.shp")

    #outputs
    viewshed_uri=os.path.join(output_dir,"vshed.tif")
    viewshed_quality_uri=os.path.join(output_dir,"vshed_qual.tif")    
    pop_stats_uri=os.path.join(output_dir,"populationStats.html")
    overlap_uri=os.path.join(output_dir,"vp_overlap.shp")

    #determining best data type for viewshed
    features = get_count_feature_set_uri(args['structure_uri'])
    if features < 2 ** 16:
        viewshed_type = gdal.GDT_UInt16
        viewshed_nodata = (2 ** 16) - 1
    elif features < 2 ** 32:
        viewshed_type = gdal.GDT_UInt32
        viewshed_nodata = (2 ** 32) - 1
    else:
        raise ValueError, "Too many structures."
    
    #clip DEM by AOI and reclass
    LOGGER.info("Clipping DEM by AOI.")

    LOGGER.debug("Projecting AOI for DEM.")
    dem_wkt = pygeoprocessing.geoprocessing.get_dataset_projection_wkt_uri(args['dem_uri'])
    pygeoprocessing.geoprocessing.reproject_datasource_uri(args['aoi_uri'], dem_wkt, aoi_dem_uri)

    LOGGER.debug("Clipping DEM by projected AOI.")
    LOGGER.debug("DEM: %s, AOI: %s", args['dem_uri'], aoi_dem_uri)
    pygeoprocessing.geoprocessing.clip_dataset_uri(args['dem_uri'], aoi_dem_uri, viewshed_dem_uri, False)

    LOGGER.info("Reclassifying DEM to account for water at sea-level and resampling to specified cell size.")
    LOGGER.debug("Reclassifying DEM so negative values zero and resampling to save on computation.")

    nodata_dem = pygeoprocessing.geoprocessing.get_nodata_from_uri(args['dem_uri'])

    def no_zeros(value):
        if value == nodata_dem:
            return nodata_dem
        elif value < 0:
            return 0
        else:
            return value

    pygeoprocessing.geoprocessing.vectorize_datasets([viewshed_dem_uri],
                                    no_zeros,
                                    viewshed_dem_reclass_uri,
                                    get_data_type_uri(viewshed_dem_uri),
                                    nodata_dem,
                                    args["cell_size"],
                                    "union")

    #calculate viewshed
    LOGGER.info("Calculating viewshed.")
    compute_viewshed_uri(viewshed_dem_reclass_uri,
             viewshed_uri,
             args['structure_uri'],
             curvature_correction,
             args['refraction'],
             args)

    LOGGER.info("Ranking viewshed.")
    #rank viewshed
    quantile_list = [25,50,75,100]
    LOGGER.debug('reclassify input %s', viewshed_uri)
    LOGGER.debug('reclassify output %s', viewshed_quality_uri)
    reclassify_quantile_dataset_uri(viewshed_uri,
                                    quantile_list,
                                    viewshed_quality_uri,
                                    viewshed_type,
                                    viewshed_nodata)

    if "pop_uri" in args:
        #tabulate population impact
        LOGGER.info("Tabulating population impact.")
        LOGGER.debug("Tabulating unaffected population.")
        nodata_pop = pygeoprocessing.geoprocessing.get_nodata_from_uri(args["pop_uri"])
        LOGGER.debug("The no data value for the population raster is %s.", str(nodata_pop))
        nodata_viewshed = pygeoprocessing.geoprocessing.get_nodata_from_uri(viewshed_uri)
        LOGGER.debug("The no data value for the viewshed raster is %s.", str(nodata_viewshed))

        #clip population
        LOGGER.debug("Projecting AOI for population raster clip.")
        pop_wkt = pygeoprocessing.geoprocessing.get_dataset_projection_wkt_uri(args['pop_uri'])
        pygeoprocessing.geoprocessing.reproject_datasource_uri(args['aoi_uri'],
                                              pop_wkt,
                                              aoi_pop_uri)

        LOGGER.debug("Clipping population raster by projected AOI.")
        pygeoprocessing.geoprocessing.clip_dataset_uri(args['pop_uri'],
                                      aoi_pop_uri,
                                      pop_clip_uri,
                                      False)
        
        #reproject clipped population
        LOGGER.debug("Reprojecting clipped population raster.")
        vs_wkt = pygeoprocessing.geoprocessing.get_dataset_projection_wkt_uri(viewshed_uri)
        reproject_dataset_uri(pop_clip_uri,
                                           vs_wkt,
                                           pop_prj_uri,
                                           get_data_type_uri(pop_clip_uri))

        #align and resample population
        def copy(value1, value2):
            if value2 == nodata_viewshed:
                return nodata_pop
            else:
                return value1
        
        LOGGER.debug("Resampling and aligning population raster.")
        pygeoprocessing.geoprocessing.vectorize_datasets([pop_prj_uri, viewshed_uri],
                                       copy,
                                       pop_vs_uri,
                                       get_data_type_uri(pop_prj_uri),
                                       nodata_pop,
                                       args["cell_size"],
                                       "intersection",
                                       ["bilinear", "bilinear"],
                                       1)
        
        pop = gdal.Open(pop_vs_uri)
        messge = "Can't open " + pop_vs_uri
        assert pop is not None, message
        pop_band = pop.GetRasterBand(1)
        messge = "Can't extract band from " + pop_vs_uri
        assert pop_band is not None, message
        vs = gdal.Open(viewshed_uri)
        vs_band = vs.GetRasterBand(1)
        message = 'population and viewshed file sizes are incompatible'
        assert vs_band.XSize >= pop_band.XSize, message
        assert vs_band.YSize >= pop_band.YSize, message

        affected_pop = 0
        unaffected_pop = 0
        for row_index in range(pop_band.YSize):
            pop_row = pop_band.ReadAsArray(0, row_index, pop_band.XSize, 1)
            vs_row = vs_band.ReadAsArray(0, row_index, vs_band.XSize, 1).astype(np.float64)

            pop_row[pop_row == nodata_pop]=0.0
            vs_row[vs_row == nodata_viewshed]=-1

            affected_pop += np.sum(pop_row[vs_row > 0])
            unaffected_pop += np.sum(pop_row[vs_row == 0])

        pop_band = None
        pop = None
        vs_band = None
        vs = None

        table="""
        <html>
        <title>Marine InVEST</title>
        <center><H1>Scenic Quality Model</H1><H2>(Visual Impact from Objects)</H2></center>
        <br><br><HR><br>
        <H2>Population Statistics</H2>

        <table border="1", cellpadding="0">
        <tr><td align="center"><b>Number of Features Visible</b></td><td align="center"><b>Population (estimate)</b></td></tr>
        <tr><td align="center">None visible<br> (unaffected)</td><td align="center">%i</td>
        <tr><td align="center">1 or more<br>visible</td><td align="center">%i</td>
        </table>
        </html>
        """

        outfile = open(pop_stats_uri, 'w')
        outfile.write(table % (unaffected_pop, affected_pop))
        outfile.close()

    #perform overlap analysis
    LOGGER.info("Performing overlap analysis.")

    LOGGER.debug("Reclassifying viewshed")

    nodata_vs_bool = 0
    def non_zeros(value):
        if value == nodata_vs_bool:
            return nodata_vs_bool
        elif value > 0:
            return 1
        else:
            return nodata_vs_bool

    pygeoprocessing.geoprocessing.vectorize_datasets([viewshed_uri],
                                    non_zeros,
                                    viewshed_reclass_uri,
                                    gdal.GDT_Byte,
                                    nodata_vs_bool,
                                    args["cell_size"],
                                    "union")

    if "overlap_uri" in args:  
        LOGGER.debug("Copying overlap analysis features.")
        pygeoprocessing.geoprocessing.copy_datasource_uri(args["overlap_uri"], overlap_uri)

        LOGGER.debug("Adding id field to overlap features.")
        id_name = "investID"
        add_id_feature_set_uri(overlap_uri, id_name)

        LOGGER.debug("Add area field to overlap features.")
        area_name = "overlap"
        add_field_feature_set_uri(overlap_uri, area_name, ogr.OFTReal)
        
        LOGGER.debug("Count overlapping pixels per area.")
        values = pygeoprocessing.geoprocessing.aggregate_raster_values_uri(
            viewshed_reclass_uri, overlap_uri, id_name, ignore_nodata=True).total

        def calculate_percent(feature):
            if feature.GetFieldAsInteger(id_name) in values:
                return (values[feature.GetFieldAsInteger(id_name)] * \
                args["cell_size"]) / feature.GetGeometryRef().GetArea()
            else:
                return 0
            
        LOGGER.debug("Set area field values.")
        set_field_by_op_feature_set_uri(overlap_uri, area_name, calculate_percent)
