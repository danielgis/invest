"""test for gdap_server"""

import sys
import os
import zipfile

import pygeoprocessing
import natcap.invest.data_platform.gdap_server
import osr
import ogr


def main():
    """entry point"""
    print 'create server'
    data_server = natcap.invest.data_platform.gdap_server.DataServer('foo.db')

    print 'adding search directories'
    data_server.add_search_directory(sys.argv[1:])

    print 'get server preview'
    data_preview = data_server.get_data_preview()
    open('data_preview.zip', 'wb').write(data_preview)

    aoi_path = r"C:\Users\Rich\Documents\svn_repos\invest-sample-data\forest_carbon_edge_effect\forest_carbon_edge_demo_aoi.shp"
    aoi_bounding_box = pygeoprocessing.get_datasource_bounding_box(aoi_path)
    vector = ogr.Open(aoi_path)
    layer = vector.GetLayer()
    aoi_projection = layer.GetSpatialRef()
    lat_lng_projection = osr.SpatialReference()
    lat_lng_projection.ImportFromEPSG(4326)  # EPSG 4326 is WGS84 lat/lng
    lat_lng_bounding_box = natcap.invest.utils.reproject_bounding_box(
        aoi_bounding_box, aoi_projection, lat_lng_projection)

    data_coverage_list = data_server.get_data_coverage(
        lat_lng_bounding_box, ['dem'])
    print '%d files to fetch ' % len(data_coverage_list)
    for data_id, _ in data_coverage_list:
        print 'fetch %s' % data_id
        data_tile = data_server.fetch_data_tile(aoi_bounding_box, data_id)
        result_zip_uri = 'data_tile.zip'
        open(result_zip_uri, 'wb').write(data_tile)
        zipfile.ZipFile(result_zip_uri, 'r').extractall('.')
        os.remove(result_zip_uri)


if __name__ == '__main__':
    main()
