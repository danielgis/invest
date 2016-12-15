"""Profile generator test."""
import natcap.invest.profile_generator


def main():
    """Entry point."""
    args = {
        'workspace_dir': 'delete_profile_generator_workspace',
        'results_suffix': 'test',
        'bathymetry_path': r"C:\Users\rpsharp\Documents\bitbucket_repos\invest\data\invest-data\Base_Data\Marine\DEMs\claybark_dem",
        'shore_height': 0.0,
        'sample_point_vector_path': r"C:\Users\rpsharp\Documents\sample_points.shp",
        'feature_id_key': 'name',
        'step_size': 90,
        'habitat_vector_path_list': [],
    }
    natcap.invest.profile_generator.execute(args)

if __name__ == '__main__':
    main()
