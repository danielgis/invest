"""Profile generator sample."""
import logging

import natcap.invest.profile_generator

logging.basicConfig(format='%(asctime)s %(name)-20s %(levelname)-8s \
%(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

LOGGER = logging.getLogger('test_profile_generator')


def main():
    """Entry point."""
    args = {
        'workspace_dir': 'delete_profile_generator_workspace_jess_debug',
        'results_suffix': 'mysuffix',
        'bathymetry_path': r"E:\Dropbox\jess_profile_debug_data\profile_data_for_jess\Andros_dem.tif",
        'shore_height': 0.0,  # shore elevation on bathy layer
        'representative_point_vector_path': r"E:\Dropbox\jess_profile_debug_data\profile_data_for_jess\Andros_beach_profile_beach.shp",
        # stepsize is (close distance step, max close distance)
        # stepsize is (far distance step, far distance definition)
        'step_size': ((10, 500), (100, 2000)),
        'smoothing_sigma': 0.0,  # sigma of gaussian filter of bathy layer
        'offshore_profile_length': 2000,
        'onshore_profile_length': 500,
        # list of tuples of the form: (shapefile path, habitat name field)
        'habitat_vector_path_list': [],
    }
    natcap.invest.profile_generator.execute(args)

if __name__ == '__main__':
    main()
