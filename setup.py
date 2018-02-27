from setuptools import setup, find_packages
from codecs import open  # To use a consistent encoding
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'DESCRIPTION.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='prost',
    version = file("VERSION").readline().rstrip(),

    description='An application to quantify and annotate miRNA expression.',
    long_description=long_description,

    # The project's main homepage.
    url='https://github.com/uoregon-postlethwait/prost',

    # Author details
    author='Jason Sydes and Peter Batzel',
    author_email='sydes@uoregon.edu',

    classifiers=[
        'Development Status :: 5 - Production/Stable',

        # Who Prost! is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        # License (must match license field above)
        'License :: Free for non-commercial use',

        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
    ],

    keywords='bioinformatics microRNA miRNA miR',

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    # packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    packages=['prost'],

    # List run-time dependencies here. These will be installed by pip when your
    # project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/technical.html#install-requires-vs-requirements-files
    install_requires=['enum34<=1.0.4', 'docopt', 'xlsxwriter', 'psutil'],

    # List additional groups of dependencies here (e.g. development dependencies).
    # You can install these using the following syntax, for example:
    # $ pip install -e .[dev,test]
    #extras_require = {
    #    'dev': ['check-manifest'],
    #    'test': ['coverage'],
    #},

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    #package_data={
    #    'sample': ['package_data.dat'],
    #},

    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages.
    # see http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    # data_files=[('my_data', ['data/data_file'])],
    #data_files = [("", ["LICENSE.txt", "DESCRIPTION.rst"])],

    # Reference the License
    license = 'LICENSE.txt',

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'prost=prost.command_line:main_wrapper',
            'prost_excel_converter=prost.excel:command_line',
        ],
    },
)

# vim: softtabstop=4:shiftwidth=4:expandtab
