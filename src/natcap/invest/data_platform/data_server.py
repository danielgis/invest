"""Data server implementation"""

import os
import sys
import tempfile
import logging
import sqlite3
import zipfile
import glob
import hashlib
import pprint

from osgeo import gdal
from osgeo import ogr
import Pyro4
import pygeoprocessing

import natcap.invest.utils

logging.basicConfig(format='%(asctime)s %(name)-20s %(levelname)-8s \
%(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

LOGGER = logging.getLogger('natcap.invest.data_platform.data_server')

Pyro4.config.SERIALIZER = 'marshal'  # lets us pass null bytes in strings


class DataServer(object):
    """Pyro4 RPC server for advertising and serving data required to run InVEST
    """

    _STATIC_DATA_TYPES = [
        'dem',
        'erodibility',
        'erosivity',
        'lulc',
        'pawc',
        'precipitation',
        'rootdepth',
        'streams',
        'watersheds',
    ]

    _DATA_TABLE_NAME = 'data_available_table'

    _DATA_SERVER_SCHEMA = [
        ('gis_type', ' text'),
        ('data_type', ' text'),
        ('bounding_box', ' text'),
        ('path', ' text'),
        ('path_hash', ' text PRIMARY KEY'),
    ]

    def __init__(self, database_filepath, search_directory_list):
        """build a server w/ files in that directory"""

        self.database_filepath = database_filepath
        filepath_directory = os.path.dirname(self.database_filepath)
        if filepath_directory != '' and not os.path.exists(filepath_directory):
            os.mkdir(os.path.dirname(self.database_filepath))

        db_connection = sqlite3.connect(self.database_filepath)
        db_cursor = db_connection.cursor()
        db_cursor.execute(
            'CREATE TABLE IF NOT EXISTS %s (%s)' % (
                self._DATA_TABLE_NAME, ','.join(
                    ['%s %s' % (field_id, modifier) for
                     field_id, modifier in self._DATA_SERVER_SCHEMA])))
        db_connection.commit()
        db_connection.close()

        raster_paths = []
        vector_paths = []

        for root, dirs, files in os.walk(search_directory_list):
            if any(x in root for x in ['.svn']):
                continue
            # check if any dirs are GIS types and prune if so
            for dir_index in reversed(xrange(len(dirs))):
                dir_path = os.path.join(root, dirs[dir_index])
                if natcap.invest.utils.is_gdal_type(dir_path):
                    raster_paths.append(dir_path)
                elif natcap.invest.utils.is_ogr_type(dir_path):
                    vector_paths.append(dir_path)
                else:
                    continue
                # if we get here, the directory is either a raster or vector
                # so no need to pick up subdirs
                del dirs[dir_index]
            for filename in files:
                file_path = os.path.join(root, filename)
                if any(file_path.endswith(suffix) for suffix in [
                        '.xml', '.hdr', '.tfw', '.gfs', '.lyr', '.xls',
                        '.pdf', '.txt']):
                    continue
                if natcap.invest.utils.is_gdal_type(file_path):
                    raster_paths.append(file_path)
                elif natcap.invest.utils.is_ogr_type(file_path):
                    vector_paths.append(file_path)

        data_hash = {}
        for gis_type, paths in [
                ('raster', raster_paths), ('vector', vector_paths)]:
            for path in paths:
                try:
                    if gis_type == 'raster':
                        bounding_box = pygeoprocessing.get_bounding_box(path)
                    elif gis_type == 'vector':
                        bounding_box = (
                            pygeoprocessing.get_datasource_bounding_box(path))
                    else:
                        raise ValueError("Unknown gis_type: %s", gis_type)
                except AttributeError:
                    LOGGER.error("Can't calculate bounds of %s", path)

                if bounding_box is None:
                    continue

                # parse out the underscored prefix of the base directory name
                # by convention, we name this to be the GIS data type
                if os.path.isfile(path):
                    data_type = (
                        os.path.basename(os.path.dirname(path)).split('_')[0])
                elif os.path.isdir(path):
                    data_type = (
                        os.path.basename(path).split('_')[0])
                else:
                    raise ValueError("What the hell is this %s" % path)

                path_hash = hashlib.sha1(path).hexdigest()
                data_hash[path_hash] = {
                    'path': path,
                    'gis_type': gis_type,
                    'bounding_box': bounding_box,
                    'data_type': data_type,
                }

        #LOGGER.debug(pprint.pformat(data_hash))

        db_connection = sqlite3.connect(self.database_filepath)
        db_cursor = db_connection.cursor()

        for path_hash, data_table in data_hash.iteritems():
            position_format = ','.join(['?'] * (len(data_table)+1))
            field_names = data_table.keys()
            field_names += ['path_hash']
            data_table['path_hash'] = path_hash

            insert_command = (
                'INSERT OR REPLACE INTO %s ' % self._DATA_TABLE_NAME +
                '(%s) VALUES (%s)' % (','.join(field_names), position_format))
            ordered_table_data = [
                str(data_table[field_id]) for field_id in field_names]
            db_cursor.execute(insert_command, ordered_table_data)
        db_connection.commit()
        db_connection.close()


    @staticmethod
    def get_server_version():
        """Returns a server version string to the client"""
        return natcap.invest.__version__

    def get_data(self, data_id, bounding_box):
        """get data from bounding box"""

        #TODO: the following is stolen from logging
        try:
            from osgeo import ogr
            from osgeo import osr
        except ImportError:
            LOGGER.error("osgeo.ogr and osgeo.osr not installed, aborting")
            raise

        #ZIP and stream the result back
        workspace_dir = tempfile.mkdtemp()
        shapefile_filename = os.path.join(
            workspace_dir, 'model_run_summary.shp')
        driver = ogr.GetDriverByName('ESRI Shapefile')

        if os.path.isfile(shapefile_filename):
            os.remove(shapefile_filename)
        datasource = driver.CreateDataSource(shapefile_filename)

        lat_lng_ref = osr.SpatialReference()
        lat_lng_ref.ImportFromEPSG(4326)  # EPSG 4326 is lat/lng

        polygon_layer = datasource.CreateLayer(
            'model_run_summary', lat_lng_ref, ogr.wkbPolygon)

        polygon_layer.CreateField(ogr.FieldDefn('n_runs', ogr.OFTInteger))
        polygon_layer.CreateField(ogr.FieldDefn('model', ogr.OFTString))

        db_connection = sqlite3.connect(self.database_filepath)
        db_cursor = db_connection.cursor()
        selection_command = (
            "SELECT model_name, bounding_box_intersection, count(model_name) "
            "FROM %s " % self._MODEL_LOG_TABLE_NAME +
            "WHERE bounding_box_intersection not LIKE 'None' "
            "GROUP BY model_name, bounding_box_intersection;")
        db_cursor.execute(selection_command)

        for line in db_cursor:
            try:
                model_name, bounding_box_string, n_runs = line
                n_runs = int(n_runs)
                bounding_box = list(
                    [float(x) for x in bounding_box_string[1:-1].split(',')])
                ring = ogr.Geometry(ogr.wkbLinearRing)
                ring.AddPoint(bounding_box[0], bounding_box[3])
                ring.AddPoint(bounding_box[0], bounding_box[1])
                ring.AddPoint(bounding_box[2], bounding_box[1])
                ring.AddPoint(bounding_box[2], bounding_box[3])
                ring.AddPoint(bounding_box[0], bounding_box[3])
                poly = ogr.Geometry(ogr.wkbPolygon)
                poly.AddGeometry(ring)
                feature = ogr.Feature(polygon_layer.GetLayerDefn())
                feature.SetGeometry(poly)
                feature.SetField('n_runs', n_runs)
                feature.SetField('model', str(model_name))
                polygon_layer.CreateFeature(feature)
            except Exception as exception:
                LOGGER.warn(
                    'unable to create a bounding box for %s (%s)', line,
                    str(exception))

        datasource.SyncToDisk()

        model_run_summary_zip_name = os.path.join(
            workspace_dir, 'model_run_summary.zip')
        with zipfile.ZipFile(model_run_summary_zip_name, 'w') as myzip:
            for filename in glob.glob(
                    os.path.splitext(shapefile_filename)[0] + '.*'):
                myzip.write(filename, os.path.basename(filename))
        model_run_summary_binary = open(model_run_summary_zip_name, 'rb').read()
        return model_run_summary_binary


def launch_data_server(data_directory, hostname, port):
    """Function to start a remote procedure call server

    Parameters:
        database_filepath (string): local filepath to the sqlite database
        hostname (string): network interface to bind to
        port (int): TCP port to bind to

    Returns:
        never"""

    daemon = Pyro4.Daemon(hostname, port)
    uri = daemon.register(
        DataServer(data_directory),
        'natcap.invest.data_platform.data_server')
    LOGGER.info(
        "natcap.invest.data_platform.data_server ready. Object uri = %s", uri)
    daemon.requestLoop()


if __name__ == '__main__':
    # attempt to set up a lock to prevent multiple running servers
    try:
        import fcntl
        PID_FILE = 'program.pid'
        FP = open(PID_FILE, 'w')
        try:
            fcntl.lockf(FP, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except IOError:
            # another instance is running
            sys.exit(0)
    except ImportError:
        LOGGER.warn(
            "Cannot import fcntl, this may allow multiple instances to run.")

    if len(sys.argv) != 4:
        LOGGER.error(
            'Error, incorrect usage try:\n%s data_directory hostname port',
            sys.argv[0])
        sys.exit(-1)
    DIRECTORY_PATH = sys.argv[1]
    HOSTNAME = sys.argv[2]
    PORT = int(sys.argv[3])
    launch_data_server(DIRECTORY_PATH, HOSTNAME, PORT)
