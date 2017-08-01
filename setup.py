#!/usr/bin/env python
from os.path import join
import sys

def copy_dir():
    dir_path = 'YOUR_DIR_HERE'
    base_dir = os.path.join('MODULE_DIR_HERE', dir_path)
    for (dirpath, dirnames, files) in os.walk(base_dir):
        for f in files:
            yield os.path.join(dirpath.split('/', 1)[1], f)



if __name__ == '__main__':
    from setuptools import setup

    try:
      from version import __version__
    except:
      __version__ = ''

    

    long_description = """
#####################################################################################
#                              Waveletspec V-1.0                                    #
#####################################################################################

A code to determine the atmospheric parameters of FGK stars from
echelle spectra. 

This code provide command line functions to load in, normalise and clean spectra from
a variety of instruments. This code also compares wavelet coefficients to a grid of
stellar models to determine the best fitting atmospherioc parameters.

Disclaimer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The spectral processing scripts are backended with iSpec scripts:

>> http://adsabs.harvard.edu/abs/2014ASInC..11...85B

I do not own their scripts and use them simply to pre-process spectra. I created a 
python 3 version of their scripts using 2to3:

>> (https://docs.python.org/2/library/2to3.html)

and changing a few module imports to their Python 3 counterparts. 


Trouble shooting
~~~~~~~~~~~~~~~~~~~~~~~~~~

1. If you are having gcc issues installing mlpy in Ubuntu, have a look at:
   https://stackoverflow.com/questions/11494179/gsl-error-while-installing-mlpy

   The basic summary is that the experimental gsl headers need installing:

   >> sudo apt-get install libgsl0-dev  

   Running the setup scripts again in a new terminal should install it.
# waveletspec
"""

    license = """	"""


    setup(name='waveletspec',
      version='1.0',
        author='Samuel Gill',
        author_email='s.gill@keele.ac.uk',
        license='GNU GPLv3',
        url='https://github.com/samgill844/waveletspec',
	packages=['waveletspec'],
    description="Wavelet analysis of echelle spectra",
    long_description = long_description,
    data_files=[ ('waveletspec/Grids', ['waveletspec/data/dummy_grid.fits'])],


    classifiers = ['Development Status :: 4 - Beta']
	)   


'''
if (sys.version_info > (3, 0)):
    data_files=[ ('data/Gridss', ['data/Grids/spectrum.fits']), ( ], 
else:
'''
