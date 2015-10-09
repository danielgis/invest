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


if __name__ == '__main__':
    main()
