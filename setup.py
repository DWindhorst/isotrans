from setuptools import setup, Extension
import numpy

module = Extension('iso_transport', sources=['example_cy.pyx'])

setup(
    name='iso_trans',
    version='',
    packages=['iso_trans'],
    url='',
    license='',
    author='poudel-b',
    author_email='',
    ext_modules=[module],
    include_dirs=[numpy.get_include()]

)
