"""Data server implementation"""

import os
import sys
import tempfile
import logging
import sqlite3
import zipfile
import hashlib
import shutil

import shapely
from osgeo import gdal
from osgeo import ogr
from osgeo import osr
import Pyro4
import pygeoprocessing

import natcap.invest.utils

logging.basicConfig(format='%(asctime)s %(name)-20s %(levelname)-8s \
%(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

LOGGER = logging.getLogger('natcap.invest.data_platform.data_server')

Pyro4.config.SERIALIZER = 'marshal'  # lets us pass null bytes in strings


def path_to_zip_string(path):
    """Recursively zips the files in the path and returns a binary string that
    is a zipfile of those data.

    Parameters:
        path (string): path in which to recursively gather files into zipfile

    Returns:
        a binary string which can be written to disk as a zipfile archive of
        the files and directories under `path`
    """

    tmp_fd, tmp_path = tempfile.mkstemp(suffix='.zip')
    abs_path = os.path.abspath(path)
    with zipfile.ZipFile(tmp_path, 'w') as out_zip_file:
        for root, _, files in os.walk(path):
            local_path = os.path.relpath(root, abs_path)
            for filename in files:
                out_zip_file.write(
                    os.path.join(root, filename),
                    os.path.join(local_path, filename))
    os.close(tmp_fd)
    tmp_file = open(tmp_path, 'rb')
    result = tmp_file.read()
    tmp_file.close()
    os.remove(tmp_path)
    return result


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

    def __init__(self, database_filepath):
        """Initialize GDAP server and its database if it doesn't exist.

        Parameters:
            database_filepath (string): path to either an existing SQLite
                database, or a path to the desired location for one.  If
                the file doesn't exist a new one is created with the
                predefined data schema.

        Returns:
            None."""

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

    def fetch_data_tile(self, bounding_box, data_id):
        """Return a binary raster or vector zipfile string clipped to the
        bounding box of the indicated id

        Parameters:
            bounding_box (list): lat/lng of bounding box of the form
                TODO: DEFINE form
            data_id (string): hash to index into the local database file

        Returns:
            binary string that can be saved as a zipfile that contains the
            clipped requested data
        """

        # Make a bounding box polygon for clipping
        bounding_box_dir = tempfile.mkdtemp()
        bounding_box_path = os.path.join(
            bounding_box_dir, 'bounding_box.geojson')
        driver = ogr.GetDriverByName('GeoJSON')
        vector = driver.CreateDataSource(bounding_box_path)
        lat_lng_projection = osr.SpatialReference()
        lat_lng_projection.ImportFromEPSG(4326)  # EPSG 4326 is WGS84 lat/lng
        polygon_layer = vector.CreateLayer(
            'bounding_box', lat_lng_projection, ogr.wkbPolygon)
        ring = ogr.Geometry(ogr.wkbLinearRing)
        x_min, y_max, x_max, y_min = bounding_box
        ring.AddPoint(x_min, y_max)
        ring.AddPoint(x_min, y_min)
        ring.AddPoint(x_max, y_min)
        ring.AddPoint(x_max, y_max)
        ring.AddPoint(x_min, y_max)
        poly = ogr.Geometry(ogr.wkbPolygon)
        poly.AddGeometry(ring)
        feature = ogr.Feature(polygon_layer.GetLayerDefn())
        feature.SetGeometry(poly)
        polygon_layer.CreateFeature(feature)
        polygon_layer = None
        ogr.DataSource.__swig_destroy__(vector)
        vector = None

        db_connection = sqlite3.connect(self.database_filepath)
        db_cursor = db_connection.cursor()
        selection_command = (
            "SELECT path "
            "FROM %s " % self._DATA_TABLE_NAME +
            " WHERE path_hash = '%s';" % data_id)
        db_cursor.execute(selection_command)
        path = db_cursor.fetchone()[0]
        db_connection.close()

        working_dir = tempfile.mkdtemp()
        out_raster_path = os.path.join(working_dir, os.path.basename(path))
        pygeoprocessing.clip_dataset_uri(
            path, bounding_box_path, out_raster_path, assert_projections=False,
            all_touched=True)

        result = path_to_zip_string(working_dir)
        # clean up intermediate result
        shutil.rmtree(working_dir)
        shutil.rmtree(bounding_box_dir)

        return result


    def add_search_directory(self, search_directory_list):
        """Recursively search through a list of directories and add any viable
        GIS data types to the database.

        Parameters:
            search_directory_list (list): A list of directory paths in which to
                recursively search through for GIS data types to index.

        Returns:
            None
        """
        raster_paths = []
        vector_paths = []

        LOGGER.info("scanning directory")
        for directory_path in search_directory_list:
            for root, dirs, files in os.walk(directory_path):
                # skip any .svn directories
                if any(x in root for x in ['.svn']):
                    continue

                # check to see if the root is a GIS type, skip any subdirs if
                # so
                if natcap.invest.utils.is_gdal_type(root):
                    raster_paths.append(root)
                    dirs.clear()
                    continue
                if natcap.invest.utils.is_ogr_type(root):
                    vector_paths.append(root)
                    dirs.clear()
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
                    # if we get here, the directory is either a raster or
                    # vector so no need to pick up subdirs
                    del dirs[dir_index]

                # test all the raw files
                for filename in files:
                    file_path = os.path.join(root, filename)
                    # ignore these commonly existing, but known to not be GIS
                    # type files
                    if any(file_path.endswith(suffix) for suffix in [
                            '.xml', '.hdr', '.tfw', '.gfs', '.lyr', '.xls',
                            '.pdf', '.txt', '.zip']):
                        continue
                    if natcap.invest.utils.is_gdal_type(file_path):
                        raster_paths.append(file_path)
                    elif natcap.invest.utils.is_ogr_type(file_path):
                        vector_paths.append(file_path)

        data_hash = {}
        LOGGER.info("calculating bounding boxes")
        for gis_type, paths in [
                ('raster', raster_paths), ('vector', vector_paths)]:
            for path in paths:
                # packing the path as an args dictionary so
                # calculate_args_bounding_box can convert to lat/lng
                bounding_box = natcap.invest.utils.calculate_args_bounding_box(
                    {'path': path})[0]
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
                    raise ValueError(
                        "This is neither a file nor directory.  WTF is it? "
                        " %s" % path)

                path_hash = hashlib.sha1(path).hexdigest()
                data_hash[path_hash] = {
                    'path': path,
                    'gis_type': gis_type,
                    'bounding_box': bounding_box,
                    'data_type': data_type,
                }

        db_connection = sqlite3.connect(self.database_filepath)
        db_cursor = db_connection.cursor()
        LOGGER.info("populating database")
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

    def get_data_coverage(self, bounding_box, data_type_list):
        """Returns a list of (data_id, datatype) pairs that can be retrieved
        with `fetch_data_tile`

        Parameters:
            bounding_box (list): in WSG84 projection:
                [x_min, y_max, x_max, y_min]
            data_type_list (list): a list of strings that can be found in
            self._STATIC_DATA_TYPES.  If empty list, ALL types are returned.

        Returns:
            List of (data_id, datatype) pairs where the data_id can be
            passed to `fetch_data_tile`
        """
        db_connection = sqlite3.connect(self.database_filepath)
        db_cursor = db_connection.cursor()
        selection_command = (
            "SELECT bounding_box, path_hash, data_type "
            "FROM %s " % self._DATA_TABLE_NAME +
            " OR ".join([
                "WHERE data_type = '%s'" % data_type for data_type in
                data_type_list]) + ';')
        db_result = db_cursor.execute(selection_command)

        bounding_coords = [
            (bounding_box[0], bounding_box[1]),
            (bounding_box[2], bounding_box[1]),
            (bounding_box[2], bounding_box[3]),
            (bounding_box[0], bounding_box[3]),
            (bounding_box[0], bounding_box[1]),
            ]
        bounding_poly = shapely.geometry.Polygon(bounding_coords)
        result_list = []
        for bounding_box_string, path_hash, data_type in db_result:
            local_bounding_box = [
                float(x) for x in bounding_box_string[1:-1].split(',')]
            local_bounding_coords = [
                (local_bounding_box[0], local_bounding_box[1]),
                (local_bounding_box[2], local_bounding_box[1]),
                (local_bounding_box[2], local_bounding_box[3]),
                (local_bounding_box[0], local_bounding_box[3]),
                (local_bounding_box[0], local_bounding_box[1]),
                ]
            local_bounding_poly = shapely.geometry.Polygon(
                local_bounding_coords)
            if local_bounding_poly.intersects(bounding_poly):
                result_list.append((path_hash, data_type))
        return result_list

    def get_data_preview(self):
        """Build an OpenLayers based HTML preview page that highlights the
        GIS data sources' bounding boxes on a global map.

        Parameters:
            none

        Returns:
            A binary string which can be interpreted as a zipfile
        get data from bounding box"""

        working_dir = tempfile.mkdtemp()

        shapefile_path = os.path.join(working_dir, 'coverage_preview.geojson')
        driver = ogr.GetDriverByName('GeoJSON')
        vector = driver.CreateDataSource(shapefile_path)

        lat_lng_projection = osr.SpatialReference()
        lat_lng_projection.ImportFromEPSG(4326)  # EPSG 4326 is WGS84 lat/lng
        mercator_projection = osr.SpatialReference()
        mercator_projection.ImportFromEPSG(3857)  # EPSG 3857 is Mercator
        lat_to_merc_transform = osr.CoordinateTransformation(
            lat_lng_projection, mercator_projection)
        polygon_layer = vector.CreateLayer(
            'coverage_preview', mercator_projection, ogr.wkbPolygon)

        polygon_layer.CreateField(ogr.FieldDefn('data_type', ogr.OFTString))
        polygon_layer.CreateField(ogr.FieldDefn('gis_type', ogr.OFTString))
        polygon_layer.CreateField(ogr.FieldDefn('path', ogr.OFTString))

        db_connection = sqlite3.connect(self.database_filepath)
        db_cursor = db_connection.cursor()
        selection_command = (
            "SELECT data_type, gis_type, bounding_box, path "
            "FROM %s " % self._DATA_TABLE_NAME + ";")
            #" WHERE data_type = '%s';" % data_type)
        LOGGER.debug(selection_command)
        db_cursor.execute(selection_command)

        for line in db_cursor:
            try:
                data_type, gis_type, bounding_box_string, path = line
                bounding_box = list(
                    [float(x) for x in bounding_box_string[1:-1].split(',')])
                ring = ogr.Geometry(ogr.wkbLinearRing)

                point_min = ogr.CreateGeometryFromWkt(
                    "POINT (%f %f)" % (bounding_box[0], bounding_box[1]))

                point_max = ogr.CreateGeometryFromWkt(
                    "POINT (%f %f)" % (bounding_box[2], bounding_box[3]))

                point_min.Transform(lat_to_merc_transform)
                point_max.Transform(lat_to_merc_transform)

                x_min = point_min.GetX()
                y_min = point_min.GetY()
                x_max = point_max.GetX()
                y_max = point_max.GetY()

                # check for wraparound
                if x_max < x_min:
                    x_max *= -1

                ring.AddPoint(x_min, y_max)
                ring.AddPoint(x_min, y_min)
                ring.AddPoint(x_max, y_min)
                ring.AddPoint(x_max, y_max)
                ring.AddPoint(x_min, y_max)

                poly = ogr.Geometry(ogr.wkbPolygon)
                poly.AddGeometry(ring)
                poly.Transform(lat_to_merc_transform)
                feature = ogr.Feature(polygon_layer.GetLayerDefn())
                feature.SetGeometry(poly)
                feature.SetField('data_type', str(data_type))
                feature.SetField('gis_type', str(gis_type))
                feature.SetField('path', str(path))
                polygon_layer.CreateFeature(feature)
            except Exception as exception:
                LOGGER.warn(
                    'unable to do this thing for %s (%s)', line,
                    str(exception))

        # close the vector so we can delete it later
        vector.SyncToDisk()
        polygon_layer = None
        ogr.DataSource.__swig_destroy__(vector)
        vector = None

        webpage_out = open(
            os.path.join(working_dir, 'coverage_preview.html'), 'w')
        webpage_out.write("""<!DOCTYPE html>
<html>
<head>
<style>

.h1, h1 {
    font-size: 36px;
}

.h3, h3 {
    font-size: 24px;
}

.h1, .h2, .h3, .h4, .h5, .h6, h1, h2, h3, h4, h5, h6 {
    font-family: inherit;
    font-weight: 500;
    line-height: 1.1;
    color: inherit;
}

body {
    font-family: "Helvetica Neue",Helvetica,Arial,sans-serif;
    font-size: 14px;
    line-height: 1.42857143;
    color: #333;
    background-color: #fff;
}

div {
    display: block;
}

.map {
    height: 500px;
    width: 100%;
    margin-bottom: 10px;
}
.topRight {
    position: absolute;
    top: 0px;
    right: 0px;
    width:80px; /* you can use % */
    height: auto;
    margin: 5px;
}

td.datatype {
    font-weight: bold;
}

td.path {
  font-family: monospace;
}

table {
  border-spacing: 5px;
}

.alert-success {
    color: #3c763d;
    background-color: #dff0d8;
    border-color: #d6e9c6;
        padding: 15px;
    margin-bottom: 20px;
    border: 1px solid transparent;
    border-radius: 4px;
}
</style>

<title>Vector layer example</title>
<script src="https://code.jquery.com/jquery-1.11.2.min.js"></script>
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.5/js/bootstrap.min.js"></script>
<script src="http://openlayers.org/en/v3.4.0/build/ol.js"></script>

</head>
<body>

  <h1>GDAP Data Coverage Overview</h1>
  <div class="row">
      <div class="map" id="map"></div>
  </div>
  <h3>Highlighted Layer Data:</h3>
  <div class="span4 offset4">
    <div id="info" class="alert alert-success">
      &nbsp;
    </div>
  </div>

  <img class="topRight" src="http://www.naturalcapitalproject.org/images/NatCap-Logo.jpg"/>

</div>
<script>
  // Vector layer from GeoJSON source
  var vectorLayer = new ol.layer.Vector({
      source: new ol.source.GeoJSON({
          url: 'coverage_preview.geojson',
          projection: 'EPSG:3857'
      }),
      style: new ol.style.Style({
        stroke: new ol.style.Stroke({
          color: 'blue',
          width: 0.1,
          lineDash: [3,3]
        }),
        fill: new ol.style.Fill({
          color: 'rgba(0, 0, 255, 0.1)'
        })
      })
  });

  // OSM tile layer
  var map = new ol.Map({
    layers: [
      new ol.layer.Tile({
        source: new ol.source.MapQuest({layer: 'osm'})
      }),
      vectorLayer
    ],
    target: 'map',
    view: new ol.View({
      center: [0, 0],
      zoom: 2.5
    })
  });

  var highlightStyleCache = {};
  var featureOverlay = new ol.FeatureOverlay({
    map: map,
    style: function(feature, resolution) {
      var text = resolution < 5000 ? feature.get('name') : '';
      if (!highlightStyleCache[text]) {
        highlightStyleCache[text] = [new ol.style.Style({
          stroke: new ol.style.Stroke({
            color: '#000',
            width: 1,
            lineDash: [5,5]
          }),
          fill: new ol.style.Fill({
            color: 'rgba(255,0,0,0.0)'
          }),
          text: new ol.style.Text({
            font: '12px Calibri,sans-serif',
            text: text,
            fill: new ol.style.Fill({
              color: '#000'
            }),
            stroke: new ol.style.Stroke({
              color: '#f0',
              width: 3
            })
          })
        })];
      }
      return highlightStyleCache[text];
    }
  });

  var displayFeatureInfo = function(pixel) {
    var highlight_features = [];
    var info_string = '<table class="feature_highlight"><tr><th>Datatype</th><th>Remote Path</th></tr>';
    var features_displayed = {};
    var keys = [];
    featureOverlay.getFeatures().clear();

    var feature = map.forEachFeatureAtPixel(pixel, function(feature, layer) {
      if (!(feature.get('path') in features_displayed)) {
        features_displayed[feature.get('path')] = feature.get('data_type');
        keys.push(feature.get('path'));
        highlight_features.push(feature);
      }
    });
    keys.sort();
    for (var i = 0; i < keys.length; i++) {
      info_string = info_string + '<tr><td class="datatype">' + features_displayed[keys[i]] + '</td><td class="path">' + keys[i] + ': ' + '</td></tr>';
    }

    if (keys.length === 0) {
      info_string = info_string + "<tr><td>n/a</td><td>n/a</td></tr>"
    }

    info_string = info_string + '</table>';
    var info = document.getElementById('info');
    info.innerHTML = info_string;

    for (var i = 0; i < highlight_features.length; i++) {
      featureOverlay.addFeature(highlight_features[i]);
    }
  };

  map.on('pointermove', function(evt) {
    if (evt.dragging) {
      return;
    }
    var pixel = map.getEventPixel(evt.originalEvent);
    displayFeatureInfo(pixel);
  });

  map.on('click', function(evt) {
    displayFeatureInfo(evt.pixel);
  });



</script>
</body>
</html>""")
        webpage_out.close()
        result = path_to_zip_string(working_dir)
        # clean up intermediate result
        shutil.rmtree(working_dir)
        return result


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
