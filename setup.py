#!/usr/bin/env python

import os

try:
    from setuptools import setup
except:
    from distutils.core import setup

requirements = [
'numpy',
'scipy',
'astropy',
'future',
]


PACKAGE_NAME = 'eidos'
__version__ = '0.1.0'

setup(name = PACKAGE_NAME,
    version = __version__,
    description = 'Modelling primary beams of radio telescope antennae',
    author = 'Khan M. B Asad',
    author_email = 'khmbasad@gmail.com',
    url = 'https://github.com/ratt-ru/eidos',
    packages = [PACKAGE_NAME],
    install_requires = requirements,
    include_package_data = True,
    package_data = { "eidos/data" : ["meerkat_coeff_dict.npy"] },
    scripts = ['bin/' + j for j in os.listdir('bin')],
    license = ['GNU GPL v3'],
    classifiers = [
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Astronomy",
    ])
