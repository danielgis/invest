"""Profile generator test."""
import logging

import natcap.invest.profile_generator

logging.basicConfig(format='%(asctime)s %(name)-20s %(levelname)-8s \
%(message)s', level=logging.DEBUG, datefmt='%m/%d/%Y %H:%M:%S ')

LOGGER = logging.getLogger('test_profile_generator')


def main():
    """Entry point."""
    args = {
        'workspace_dir': 'delete_profile_generator_workspace',
        'results_suffix': 'test',
        'bathymetry_path': r"E:\repositories\bitbucket_repos\invest\data\invest-data\Base_Data\Marine\DEMs\claybark_dem",
        'shore_height': 0.0,
        'sample_point_vector_path': r"C:\Users\Daddy\Documents\representative_profile_points.shp",
        'feature_id_key': 'name',
        'step_size': 90,
        'habitat_vector_path_list': [],
    }
    natcap.invest.profile_generator.execute(args)

if __name__ == '__main__':
    main()
