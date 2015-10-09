"""test for gdap_server"""

import sys

import natcap.invest.data_platform.gdap_server


def main():
    """entry point"""
    print 'create server'
    data_server = natcap.invest.data_platform.gdap_server.DataServer('foo.db')

    print 'adding search directories'
    data_server.add_search_directory(sys.argv[1:])

    print 'get server preview'
    data_preview = data_server.get_data_preview()
    open('data_preview.zip', 'wb').write(data_preview)

    bounding_box = [
        94.99958345439518, 25.00041712672993, 100.00041678772851,
        19.999583793396596]
    data_id = "3d921d0fe6d31d2a76fa1e3923098bf6ee0f62da"
    data_tile = data_server.fetch_data_tile(bounding_box, data_id)
    open('data_tile.zip', 'wb').write(data_tile)


if __name__ == '__main__':
    main()
