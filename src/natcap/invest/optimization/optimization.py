"""Optimization Prototype."""


def optimization_prototype():
    """Optimization Prototype."""
    base_dir = r'C:/Users/rpsharp/Documents/bitbucket_repos/invest/data/invest-data/Base_Data/Freshwater'

    # Optimizer that works by changing habitat

    # * specify marginal gain map generators (say whether to min/max)
    #   * specify marginal gain map changers (vector, data, and change)
    #       * specify constraints?

    # algorithm
    # * do a change on the highest marginal gain with the least constraint
    # * if requested, update marginal gain maps
    # * if no change possible given constraints, then quit

    # report
    # * map of changes?
    # * total "cost"

    # ex:
    #   * min(sed_export)
    #   * max(ag_production) ag_production <= 10000
    #   * max(carbon)
    #   * max(PUD)
    #   habitat change vector per map
    #       * vector -> what to change
    #           * (sed) grid cells [zones] -> lulc map -> convert to lulc code "5" (forest)
    #               * update marginal gain burn solution to LULC and recalculate?
    #           * (ag) grid cells [zones] -> lulc map -> convert to lulc code "12" (ag)
    #               * update marginal gain burn solution to LULC and recalculate?
    #           * (carbon) grid cells [zones] -> lulc map -> convert to lulc code "5" (forest)
    #               * update marginal gain burn solution to LULC and recalculate?
    #           * (PUD) vector [park_locations] -> park location map -> pick existance
    #               * select park, remove overlaps, and re-run PUD for each vector
    #
    # there are sed_export_marginal, ag_production_marginal, carbon_marginal
    #   don't worry how these are calculated, they can all be different
    #


    optimizer_args = {
    }

if __name__ == '__main__':
    optimization_prototype()


"""'input_list': [
            ('raster', 'dem', 'path...'),
            ('raster', 'erosivity', 'path...'),
            ('table', 'biophysical_table': os.path.join('biophysical_table.csv'),
        u'dem_path': u'C:/Users/rpsharp/Documents/bitbucket_repos/invest/data/invest-data/Base_Data/Freshwater/dem',
        u'drainage_path': u'',
        u'erodibility_path': u'C:/Users/rpsharp/Documents/bitbucket_repos/invest/data/invest-data/Base_Data/Freshwater/erodibility',
        u'erosivity_path': u'C:/Users/rpsharp/Documents/bitbucket_repos/invest/data/invest-data/Base_Data/Freshwater/erosivity',
        u'ic_0_param': u'0.5',
        u'k_param': u'2',
        u'lulc_path': u'C:/Users/rpsharp/Documents/bitbucket_repos/invest/data/invest-data/Base_Data/Freshwater/landuse_90',
        u'sdr_max': u'0.8',
        u'threshold_flow_accumulation': u'1000',
        u'watersheds_path': u'C:/Users/rpsharp/Documents/bitbucket_repos/invest/data/invest-data/Base_Data/Freshwater/watersheds.shp',
        u'workspace_dir': u'C:\\Users\\rpsharp/Documents/delete_sed_work_drainag5',
}"""