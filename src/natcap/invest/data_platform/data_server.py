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
from osgeo import osr
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

        LOGGER.info("scanning directory")
        for root, dirs, files in os.walk(search_directory_list):
            if any(x in root for x in ['.svn']):
                continue
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
                # if we get here, the directory is either a raster or vector
                # so no need to pick up subdirs
                del dirs[dir_index]
            for filename in files:
                file_path = os.path.join(root, filename)
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

    def get_data_preview(self):
        """get data from bounding box"""

        #for data_type in self._STATIC_DATA_TYPES:
        #ZIP and stream the result back
        shapefile_filename = 'coverage_preview.geojson'
        driver = ogr.GetDriverByName('GeoJSON')

        if os.path.isfile(shapefile_filename):
            os.remove(shapefile_filename)
        datasource = driver.CreateDataSource(shapefile_filename)

        lat_lng_projection = osr.SpatialReference()
        lat_lng_projection.ImportFromEPSG(4326)
        out_projection = osr.SpatialReference()
        #out_projection.ImportFromEPSG(4326)  # EPSG 4326 is lat/lng
        out_projection.ImportFromEPSG(3857)  # EPSG 4326 is lat/lng
        transform = osr.CoordinateTransformation(lat_lng_projection, out_projection)
        polygon_layer = datasource.CreateLayer(
            'coverage_preview', out_projection, ogr.wkbPolygon)

        polygon_layer.CreateField(
            ogr.FieldDefn('data_type', ogr.OFTString))
        polygon_layer.CreateField(
            ogr.FieldDefn('gis_type', ogr.OFTString))
        polygon_layer.CreateField(
            ogr.FieldDefn('path', ogr.OFTString))

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

                point_min.Transform(transform)
                point_max.Transform(transform)

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
                poly.Transform(transform)
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
        datasource.SyncToDisk()
        webpage_out = open('overview.html', 'w')
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
