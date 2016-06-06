"""GDAP client implementation."""

import logging
import urllib

import Pyro4

logging.basicConfig(format='%(asctime)s %(name)-20s %(levelname)-8s \
%(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

LOGGER = logging.getLogger('natcap.invest.data_platform.data_server')

INVEST_USAGE_LOGGER_URL = (
    'http://data.naturalcapitalproject.org/server_registry/data_hosts/')

Pyro4.config.SERIALIZER = 'marshal'  # lets us pass null bytes in strings


def _get_data_server(path=None):
    """Returns a remote procedure call logging server from the
    https://bitbucket.org/natcap/natcap_model_logger project.

    Parameters:
        path (string): A Pyro4 compatible url for getting a Proxy object.

    """

    if path is None:
        path = urllib.urlopen(INVEST_USAGE_LOGGER_URL).read().rstrip()
    logging_server = Pyro4.Proxy(path)
    return logging_server
