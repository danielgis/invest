"""Profile generator test."""
import logging

import natcap.invest.profile_generator

logging.basicConfig(format='%(asctime)s %(name)-20s %(levelname)-8s \
%(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

LOGGER = logging.getLogger('test_profile_generator')


def main():
    """Entry point."""
    args = {
        'workspace_dir': 'profile_generator_workspace',
        'results_suffix': 'test',
        'bathymetry_path': r"claybark_dem",
        'shore_height': 0.0,
        'representative_point_vector_path': r"representative_profile_points.shp",
        'step_size': 20,
        'profile_length': 600,
        'habitat_vector_path_list': [
            (r"sample_claybark_hab_a.shp", 'name')],
    }
    natcap.invest.profile_generator.execute(args)

if __name__ == '__main__':
    main()
