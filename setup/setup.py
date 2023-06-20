
from distutils.core import setup
from distutils.extension import Extension
import sys

import numpy
from Cython.Build import cythonize

pth = sys.path[-1]
extensions = Extension('iso_transport', sources=[pth + '\example_cy.pyx'])
setup(ext_modules=cythonize(extensions))


