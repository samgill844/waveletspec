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
