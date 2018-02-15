"""setup.py module for natcap.invest

InVEST - Integrated Valuation of Ecosystem Services and Tradeoffs

Common functionality provided by setup.py:
    build_sphinx

For other commands, try `python setup.py --help-commands`
"""
from Cython.Build import cythonize
from setuptools.extension import Extension
from setuptools import setup
import numpy


# Read in requirements.txt and populate the python readme with the non-comment
# contents.
_REQUIREMENTS = [req for req in open('requirements.txt').readlines()
                 if not req.startswith('#') and len(req) > 0 and not
                 req.startswith('hg+')]
README = open('README_PYTHON.rst').read().format(
    requirements='\n'.join(['    ' + r for r in _REQUIREMENTS]))


setup(
    name='natcap.invest',
    description="InVEST Ecosystem Service models",
    long_description=README,
    maintainer='James Douglass',
    maintainer_email='jdouglass@stanford.edu',
    url='http://bitbucket.org/natcap/invest',
    namespace_packages=['natcap'],
    packages=[
        'natcap',
        'natcap.invest',
        'natcap.invest.coastal_blue_carbon',
        'natcap.invest.coastal_vulnerability',
        'natcap.invest.finfish_aquaculture',
        'natcap.invest.fisheries',
        'natcap.invest.habitat_risk_assessment',
        'natcap.invest.hydropower',
        'natcap.invest.ui',
        'natcap.invest.ndr',
        'natcap.invest.overlap_analysis',
        'natcap.invest.recreation',
        'natcap.invest.reporting',
        'natcap.invest.routing',
        'natcap.invest.scenario_generator',
        'natcap.invest.scenic_quality',
        'natcap.invest.seasonal_water_yield',
        'natcap.invest.wave_energy',
        'natcap.invest.wind_energy',
        'natcap.invest.pygeoprocessing_0_3_3',
        'natcap.invest.pygeoprocessing_0_3_3.routing',
        'natcap.invest.pygeoprocessing_0_3_3.dbfpy',
        'natcap.invest.pygeoprocessing_0_3_3.testing',
    ],
    package_dir={
        'natcap': 'src/natcap'
    },
    use_scm_version={'version_scheme': 'post-release',
                     'local_scheme': 'node-and-date'},
    include_package_data=True,
    install_requires=_REQUIREMENTS,
    setup_requires=['setuptools_scm', 'numpy', 'cython'],
    license='BSD',
    zip_safe=False,
    keywords='gis invest',
    classifiers=[
        'Intended Audience :: Developers',
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 2 :: Only',
        'License :: OSI Approved :: BSD License',
        'Topic :: Scientific/Engineering :: GIS'
    ],
    ext_modules=cythonize([
        Extension(
            name="natcap.invest.recreation.out_of_core_quadtree",
            sources=[
                'src/natcap/invest/recreation/out_of_core_quadtree.pyx'],
            include_dirs=[numpy.get_include()],
            language="c++"),
        Extension(
            name="scenic_quality_cython_core",
            sources=[
                'src/natcap/invest/scenic_quality/scenic_quality_cython_core.pyx'],
            include_dirs=[numpy.get_include()],
            language="c++"),
        Extension(
            name="ndr_core",
            sources=['src/natcap/invest/ndr/ndr_core.pyx'],
            include_dirs=[numpy.get_include()],
            language="c++"),
        Extension(
            name="seasonal_water_yield_core",
            sources=['src/natcap/invest/seasonal_water_yield/seasonal_water_yield_core.pyx'],
            include_dirs=[numpy.get_include()],
            language="c++"),
        Extension(
            name="natcap.invest.pygeoprocessing_0_3_3.geoprocessing_core",
            sources=[
                'src/natcap/invest/pygeoprocessing_0_3_3/geoprocessing_core.pyx'],
            include_dirs=[numpy.get_include()],
            language="c++"),
        Extension(
            name="natcap.invest.pygeoprocessing_0_3_3.routing.routing_core",
            sources=[
                'src/natcap/invest/pygeoprocessing_0_3_3/routing/routing_core.pyx'],
            include_dirs=[numpy.get_include()],
            language="c++")
    ]),
    entry_points={
        'console_scripts': [
            'invest = natcap.invest.cli:main'
        ],
    },
    extras_require={
        'ui': ('qtpy>1.3', 'qtawesome', 'faulthandler'),
    },
    package_data={
        'natcap.invest.reporting': [
            'reporting_data/*.js',
            'reporting_data/*.css',
        ],
        'natcap.invest.scenario_generator': [
            '*.js',
        ],
        'natcap.invest.wave_energy': [
            'wave_energy_scripts/*.sh',
            'wave_energy_scripts/*.txt'
        ],
    }
)
