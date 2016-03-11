"""InVEST optimization framework."""

from enum import Enum


def optimize_example():
    """Optimize example."""
    InputType = Enum('InputType', 'raster vector table scalar')

    MetricType = Enum(
        'MetricType',
        'cost total_sed_export total_nut_export total_carbon_stocks')

    constraint_a = 'MetricType.cost <= 1000'
    score = 'total_sed_export + total_carbon_stocks'



    print MetricType

if __name__ == '__main__':
    optimize_example()
