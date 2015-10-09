"""test for gdap_server"""

import sys

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
    print data_coverage_list
    #data_tile = data_server.fetch_data_tile(bounding_box, data_id)
    #open('data_tile.zip', 'wb').write(data_tile)


if __name__ == '__main__':
    main()
