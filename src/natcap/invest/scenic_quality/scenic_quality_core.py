import math
import os
import numpy as np
import collections
import logging
from osgeo import gdal

import pygeoprocessing.geoprocessing
import scenic_quality_cython_core

logging.getLogger('pygeoprocessing.geoprocessing').setLevel(logging.WARNING)
logging.getLogger('raster_cython_utils').setLevel(logging.WARNING)

logging.basicConfig(format='%(asctime)s %(name)-15s %(levelname)-8s \
    %(message)s', level=logging.WARNING, datefmt='%m/%d/%Y %H:%M:%S ')
LOGGER = logging.getLogger('scenic_quality_core')

# Linked cells used for the active pixels
linked_cell_factory = collections.namedtuple('linked_cell', \
    ['previous', 'next', 'distance', 'visibility'])

# Links to the cells
cell_link_factory = collections.namedtuple('cell_link', \
    ['top', 'right', 'bottom', 'level', 'distance'])

def print_sweep_line(sweep_line):
    if sweep_line:
        pixel = sweep_line['closest']
        while pixel is not None:
            print( \
                'pixel', pixel['distance'], \
                'visibility', pixel['visibility'], 'next', \
                (pixel['next'] 
                    if pixel['next'] is None 
                    else pixel['next']['distance']))
            pixel = pixel['next']

def print_node(node):
    """Printing a node by displaying its 'distance' and 'next' fields"""
    print(str(None) if node is None else (node['distance'] + '-' + \
    (str(None) if node['next'] is None else node['next']['distance'])))


def update_visible_pixels(active_pixels, I, J, visibility_map):
    """Update the array of visible pixels from the active pixel's visibility
    
            Inputs:
                -active_pixels: a linked list of dictionaries containing the
                following fields:
                    -distance: distance between pixel center and viewpoint
                    -visibility: an elevation/distance ratio used by the
                    algorithm to determine what pixels are bostructed
                    -index: pixel index in the event stream, used to find the
                    pixel's coordinates 'i' and 'j'.
                    -next: points to the next pixel, or is None if at the end
                The linked list is implemented with a dictionary where the
                pixels distance is the key. The closest pixel is also
                referenced by the key 'closest'.
                -I: the array of pixel rows indexable by pixel['index']
                -J: the array of pixel columns indexable by pixel['index']
                -d: distance in meters
                -visibility_map: a python array the same size as the DEM
                with 1s for visible pixels and 0s otherwise. Viewpoint is
                always visible.
            
            Returns nothing"""
    # Update visibility and create a binary map of visible pixels
    # -Look at visibility from closer pixels out, keep highest visibility
    # -A pixel is not visible if its visibility <= highest visibility so far
    if not active_pixels:
        return

    pixel = active_pixels['closest']
    max_visibility = -1000000.
    while pixel is not None:
        # Pixel is visible
        if pixel['offset'] > max_visibility:
            #print('pixel is visible:', pixel['offset'], '>',max_visibility)
            visibility = 1
            visibility = pixel['offset'] - max_visibility
        # Pixel is not visible
        else:
            #print('pixel is not visible:', pixel['offset'], '<=',max_visibility)
            visibility = 0
            visibility = pixel['offset'] - max_visibility
        # Update max_visibility for this pixel
        if pixel['visibility'] > max_visibility:
            #print('Updating visibility', pixel['visibility'], '>', max_visibility)
	    max_visibility = pixel['visibility']
        # Update the visibility map for this pixel
        index = pixel['index']
        i = I[index]
        j = J[index]
        if visibility_map[i, j] == 0:
            visibility_map[i, j] = visibility
        pixel = pixel['next']


def find_active_pixel(sweep_line, distance):
    """Find an active pixel based on distance. Return None if can't be found"""
    if 'closest' in sweep_line:
        # Get information about first pixel in the list
        pixel = sweep_line[sweep_line['closest']['distance']] # won't change
        # Move on to next pixel if we're not done
        while (pixel is not None) and \
            (pixel['distance'] < distance):
            pixel = pixel['next']
        # We reached the end and didn't find anything
        if pixel is None:
            return None
        # We didn't reach the end: either pixel doesn't exist...
        if pixel['distance'] != distance:
            return None
        # ...or we found it
        else:
            return pixel
    else:
        return None


def remove_active_pixel(sweep_line, distance):
    """Remove a pixel based on distance. Do nothing if can't be found."""
    if 'closest' in sweep_line:
        # Get information about first pixel in the list
        previous = None
        pixel = sweep_line[sweep_line['closest']['distance']] # won't change
        # Move on to next pixel if we're not done
        while (pixel is not None) and \
            (pixel['distance'] < distance):
            previous = pixel
            pixel = pixel['next']
        # We reached the end and didn't find anything
        if pixel is None:
            return sweep_line
        # We didn't reach the end: either pixel doesn't exist:
        if pixel['distance'] != distance:
            return sweep_line
        # Or we found the value we want to delete
        # Make the previous element point to the next
        # We're at the beginning of the list: update the list's first element
        if previous is None:
            # No next pixel: we have to delete 'closest'
            if pixel['next'] is None:
                del sweep_line['closest']
            # Otherwise, update it
            else:
                sweep_line['closest'] = pixel['next']
        # We're not at the beginning of the list: only update previous
        else:
            previous['next'] = sweep_line[distance]['next']
        # Remove the value from the list
        del sweep_line[distance]
    return sweep_line


def add_active_pixel(sweep_line, index, distance, visibility, offset):
    """Add a pixel to the sweep line in O(n) using a linked_list of
    linked_cells."""
    #print('adding ' + str(distance) + ' to python list')
    #print_sweep_line(sweep_line)
    # Make sure we're not creating any duplicate
    message = 'Duplicate entry: the value ' + str(distance) + ' already exist'
    assert distance not in sweep_line, message
    new_pixel = \
    {'next':None, 'index':index, 'distance':distance, \
        'visibility':visibility, 'offset':offset}
    if 'closest' in sweep_line:
        # Get information about first pixel in the list
        previous = None
        pixel = sweep_line[sweep_line['closest']['distance']] # won't change
        # Move on to next pixel if we're not done
        while (pixel is not None) and \
            (pixel['distance'] < distance):
            previous = pixel
            pixel = pixel['next']
        # 1- Make the current pixel points to the next one
        new_pixel['next'] = pixel
        # 2- Insert the current pixel in the sweep line:
        sweep_line[distance] = new_pixel
        # 3- Make the preceding pixel point to the current one
        if previous is None:
            sweep_line['closest'] = new_pixel
        else:
            sweep_line[previous['distance']]['next'] = sweep_line[distance]
    else:
        sweep_line[distance] = new_pixel
        sweep_line['closest'] = new_pixel
    return sweep_line


def list_extreme_cell_angles(array_shape, viewpoint_coords, max_dist):
    """List the minimum and maximum angles spanned by each cell of a
        rectangular raster if scanned by a sweep line centered on
        viewpoint_coords.
    
        Inputs:
            -array_shape: a shape tuple (rows, cols) as is created from
                calling numpy.ndarray.shape()
            -viewpoint_coords: a 2-tuple of coordinates similar to array_shape
            where the sweep line originates
            -max_dist: maximum viewing distance
            
        returns a tuple (min, center, max, I, J) with min, center and max 
        Nx1 numpy arrays of each raster cell's minimum, center, and maximum 
        angles and coords as two Nx1 numpy arrays of row and column of the 
        coordinate of each point.
    """
    viewpoint = np.array(viewpoint_coords)

    pi = math.pi
    two_pi = 2. * pi
    rad_to_deg = 180.0 / pi
    deg_to_rad = 1.0 / rad_to_deg

    # Each row correspond to an angular sector (angle, min & max corners):
    # row 0: angle == 0         (min: lower left,  max: upper left)
    # row 1: 0 < angle < PI/2   (min: lower right, max: upper left)
    # row 2: angle == PI/2      (min: lower right, max: lower left)
    # row 3: PI/2 < angle < PI  (min: upper right, max: lower left)
    # row 4: angle == PI        (min: upper right, max: lower right)
    # row 5: PI < angle < 3PI/2 (min: upper left,  max: lower right)
    # row 6: angle == 3PI/2     (min: upper left,  max: upper right)
    # row 7: 3PI/2 < angle < 2PI(min: lower left,  max: upper right)
    extreme_cell_points = [ \
    {'min_angle':[0.5, -0.5], 'max_angle':[-0.5, -0.5]}, \
    {'min_angle':[0.5, 0.5], 'max_angle':[-0.5, -0.5]}, \
    {'min_angle':[0.5, 0.5], 'max_angle':[0.5, -0.5]}, \
    {'min_angle':[-0.5, 0.5], 'max_angle':[0.5, -0.5]}, \
    {'min_angle':[-0.5, 0.5], 'max_angle':[0.5, 0.5]}, \
    {'min_angle':[-0.5, -0.5], 'max_angle':[0.5, 0.5]}, \
    {'min_angle':[-0.5, -0.5], 'max_angle':[-0.5, 0.5]}, \
    {'min_angle':[0.5, -0.5], 'max_angle':[-0.5, 0.5]}]

    min_angles = []
    angles = []
    max_angles = []
    I = []
    J = []
    max_dist_sq = max_dist**2 if max_dist > 0 else 1000000000
    cell_count = array_shape[0]*array_shape[1]
    current_cell_id = 0
    discarded_cells = 0
    for row in range(array_shape[0]):
        for col in range(array_shape[1]):
            if (cell_count > 1000) and \
                (current_cell_id % (cell_count/1000)) == 0:
                progress = round(float(current_cell_id) / cell_count * 100.,1)
                print(str(progress) + '%')
            # Skip if cell is too far
            cell = np.array([row, col])
            viewpoint_to_cell = cell - viewpoint
            if np.sum(viewpoint_to_cell**2) > max_dist_sq:
                discarded_cells += 1
                continue
            # Skip if cell falls on the viewpoint
            if (row == viewpoint[0]) and (col == viewpoint[1]):
                discarded_cells += 1
                continue
            I.append(row)
            J.append(col)
            # Compute the angle of the cell center
            angle = np.arctan2(-viewpoint_to_cell[0], viewpoint_to_cell[1])
            angles.append((angle + two_pi) % two_pi)
            # find index in extreme_cell_points that corresponds to the current
            # angle to compute the offset from cell center
            sector = int(4. * angles[-1] / two_pi) * 2
            if np.amin(np.absolute(viewpoint_to_cell)) > 0:
                sector += 1
            min_corner_offset = \
                np.array(extreme_cell_points[sector]['min_angle'])
            max_corner_offset = \
                np.array(extreme_cell_points[sector]['max_angle'])
            # Use the offset to compute extreme angles
            min_corner = viewpoint_to_cell + min_corner_offset
            min_angle = np.arctan2(-min_corner[0], min_corner[1])
            min_angles.append((min_angle + two_pi) % two_pi) 

            max_corner = viewpoint_to_cell + max_corner_offset
            max_angle = np.arctan2(-max_corner[0], max_corner[1])
            max_angles.append((max_angle + two_pi) % two_pi)
            current_cell_id += 1
    # Create a tuple of ndarray angles before returning
    min_angles = np.array(min_angles)
    angles = np.array(angles)
    max_angles = np.array(max_angles)
    I = np.array(I)
    J = np.array(J)
    return (min_angles, angles, max_angles, I, J)


def get_perimeter_cells(array_shape, viewpoint, max_dist=-1):
    """Compute cells along the perimeter of an array.

        Inputs:
            -array_shape: tuple (row, col) as ndarray.shape containing the
            size of the array from which to compute the perimeter
            -viewpoint: tuple (row, col) indicating the position of the
            observer
            -max_dist: maximum distance in pixels from the center of the array.
            Negative values are ignored (same effect as infinite distance).
            
        Returns a tuple (rows, cols) of the cell rows and columns following
        the convention of numpy.where() where the first cell is immediately
        right to the viewpoint, and the others are enumerated clockwise."""
    message = 'viewpoint ' + str(viewpoint) + \
        ' outside of the raster boundaries ' + str(array_shape) + '.'
    assert viewpoint[0] >= 0, message 
    assert viewpoint[0] < array_shape[0], message
    assert viewpoint[1] >= 0, message
    assert viewpoint[1] < array_shape[1], message

    # Adjust max_dist to a very large value if negative
    if max_dist < 0:
        max_dist = 1000000000
    # Limit the perimeter cells to the intersection between the array's bounds
    # and the axis-aligned bounding box that extends around viewpoint out to
    # max_dist
    i_min = max(viewpoint[0] - max_dist, 0)
    i_max = min(viewpoint[0] + max_dist + 1, array_shape[0])
    j_min = max(viewpoint[1] - max_dist, 0)
    j_max = min(viewpoint[1] + max_dist + 1, array_shape[1])
    # list all perimeter cell center angles
    row_count = i_max - i_min 
    col_count = j_max - j_min
    # Create top row, except cell (0,0)
    rows = np.zeros(col_count - 1)
    cols = np.array(range(col_count-1, 0, -1))
    # Create left side, avoiding repeat from top row
    rows = np.concatenate((rows, np.array(range(row_count -1))))
    cols = np.concatenate((cols, np.zeros(row_count - 1)))
    # Create bottom row, avoiding repat from left side
    rows = np.concatenate((rows, np.ones(col_count - 1) * (row_count -1)))
    cols = np.concatenate((cols, np.array(range(col_count - 1))))
    # Create last part of the right side, avoiding repeat from bottom row
    rows = np.concatenate((rows, np.array(range(row_count - 1, 0, -1))))
    cols = np.concatenate((cols, np.ones(row_count - 1) * (col_count - 1)))
    # Add a offsets to rows & cols to be consistent with viewpoint
    rows += i_min
    cols += j_min
    # Roll the arrays so the first point's angle at (rows[0], cols[0]) is 0
    rows = np.roll(rows, viewpoint[0] - i_min).astype(np.int64)
    cols = np.roll(cols, viewpoint[0] - i_min).astype(np.int64)
    return (rows, cols)


def cell_angles(cell_coords, viewpoint):
    """Compute angles between cells and viewpoint where 0 angle is right of
    viewpoint.
        Inputs:
            -cell_coords: coordinate tuple (rows, cols) as numpy.where()
            from which to compute the angles
            -viewpoint: tuple (row, col) indicating the position of the
            observer. Each of row and col is an integer.
            
        Returns a sorted list of angles"""
    rows, cols = cell_coords
    # List the angles between each perimeter cell
    two_pi = 2.0 * math.pi
    r = np.array(range(rows.size))
    p = (rows[r] - viewpoint[0], cols[r] - viewpoint[1])
    angles = (np.arctan2(-p[0], p[1]) + two_pi) % two_pi
    return angles


def viewshed(input_array, cell_size, visibility_map, perimeter_cells, coordinates, \
    angles, v, viewshed_shape, row_min, col_min, \
    add_events, center_events, remove_events, I, J, \
    arg_min, arg_max, arg_center, \
    coord, distances_sq, distances, visibility, offset_visibility, \
    obs_elev=1.75, tgt_elev=0.0, max_dist=-1.0, \
    refraction_coeff=None, alg_version='cython', path = '', index = None):
    """Compute the viewshed for a single observer. 
        Inputs: 
            -input_array: a numpy array of terrain elevations
            -nodata: input_array's nodata value
            -coordinates: tuple (east, north) of coordinates of viewing
                position
            -obs_elev: observer elevation above the raster map.
            -tgt_elev: offset for target elevation above the ground. Applied to
                every point on the raster
            -max_dist: maximum visibility radius. By default infinity (-1), 
            -cell_size: cell size in meters (integer)
            -refraction_coeff: refraction coefficient (0.0-1.0), not used yet
            -alg_version: name of the algorithm to be used. Either 'cython'
                (default) or 'python'.
            -path: file path used to save data for debugging purposes

        Returns the visibility map for the DEM as a numpy array"""
    if alg_version is 'python':
        sweep_through_angles( \
            coordinates, \
            angles, add_events, center_events, remove_events,\
            I, J, distances, offset_visibility, visibility, \
            visibility_map)
    else:
        scenic_quality_cython_core.sweep_through_angles( \
            np.array(coordinates).astype(np.int64), \
            np.array([perimeter_cells[0], \
            perimeter_cells[1]]), angles, \
            add_events, center_events, remove_events, \
            arg_min, arg_center, arg_max, \
            coord, distances, offset_visibility, visibility, \
            visibility_map, path, index)

    # Set the viewpoint visible as a convention
    visibility_map[coordinates] = 1

    # Return the output to save computation
    return visibility_map

def active_pixel_index(O, P, E):
    return scenic_quality_cython_core._active_pixel_index(O, P, E)

def sweep_through_angles(viewpoint, angles, add_events, center_events, \
    remove_events, I, J, distances, offset_visibility, visibility, \
    visibility_map):
    """Update the active pixels as the algorithm consumes the sweep angles"""
    angle_count = len(angles)
    # 4- build event lists
    add_event_id = 0
    add_event_count = add_events.size
    center_event_id = 0
    center_event_count = center_events.size
    remove_event_id = 0
    remove_event_count = remove_events.size
    # 5- Sort event lists
    arg_min = np.argsort(add_events)
    arg_center = np.argsort(center_events)
    arg_max = np.argsort(remove_events)

    # Updating active cells
    active_line = {}
    # 1- add cells at angle 0
    #LOGGER.debug('Creating python event stream')
    # Create cell_center_events
    while (center_event_id < center_event_count) and \
        (center_events[arg_center[center_event_id]] < angles[1]):
        c = arg_center[center_event_id]
        d = distances[c]
        v = visibility[c]
        o = offset_visibility[c]
        active_line = add_active_pixel(active_line, c, d, v, o)
        center_event_id += 1
    # The sweep line is current, now compute pixel visibility
    update_visible_pixels(active_line, I, J, visibility_map)
    
    # 2- loop through line sweep angles:
    for a in range(angle_count-2):
    #   2.2- remove cells
        while (remove_event_id < remove_event_count) and \
            (remove_events[arg_max[remove_event_id]] <= angles[a+1]):
            c = arg_max[remove_event_id]
            d = distances[c]
            v = visibility[c]
            active_line = remove_active_pixel(active_line, d)
            remove_event_id += 1
        #   2.1- add cells
        while (add_event_id < add_event_count) and \
            (add_events[arg_min[add_event_id]] < angles[a+1]):
            c = arg_min[add_event_id]
            d = distances[c]
            v = visibility[c]
            o = offset_visibility[c]
            active_line = add_active_pixel(active_line, c, d, v, o)
            add_event_id += 1
        # The sweep line is current, now compute pixel visibility
        update_visible_pixels(active_line, I, J, visibility_map)


def execute(args):
    """Entry point for scenic quality core computation.
    
        Inputs:
        
        Returns
    """
    pass
