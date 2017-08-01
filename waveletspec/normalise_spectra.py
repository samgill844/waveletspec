from __future__ import print_function
import os,sys

from astropy.io import fits
import waveletspec,os,sys
import waveletspec,os,sys
this_dir, this_filename = os.path.split(waveletspec.__file__)
print('ISPEC_1: {}/iSpec_v20160930_py2'.format(this_dir))
# First search for the iSpec tar (eithe Python2 or 3 version)
try:
    if (sys.version_info < (3, 0)):
        os.stat(this_dir+'/iSpec_v20160930_py2')
        print('Found python 2 version of iSpec')
        ispec_dir = this_dir+'/iSpec_v20160930_py2'
        sys.path.insert(0, os.path.abspath((this_dir+'/iSpec_v20160930_py2')))
    if (sys.version_info > (3, 0)):
        os.stat(this_dir+'/iSpec_v20160930_py3')
        ispec_dir = this_dir+'/iSpec_v20160930_py3'
        print('Found python 3 version of iSpec')
        sys.path.insert(0, os.path.abspath((this_dir+'/iSpec_v20160930_py3')))
except OSError as e:
    # export the tar
    print('I cant find the iSpec directory, searching for the tar files...')
    if (sys.version_info > (3, 0)):
        os.system('tar -xf {}/iSpec_v20160930_py3.tar.gz -C {}'.format(this_dir,this_dir))
        print('Exctracted to {}/iSpec_v20160930_py3'.format(this_dir))
        ispec_dir = this_dir+'/iSpec_v20160930_py3'
        sys.path.insert(0, os.path.abspath((this_dir+'/iSpec_v20160930_py3')))
        
    elif (sys.version_info < (3, 0)):
        os.system('tar -xf {}/iSpec_v20160930_py2.tar.gz -C {}'.format(this_dir,this_dir))
        print('Exctracted to {}/iSpec_v20160930_py2'.format(this_dir))
        sys.path.insert(0, os.path.abspath((this_dir+'/iSpec_v20160930_py2')))
        ispec_dir = this_dir+'/iSpec_v20160930_py2'
    else:
        print('I CANT FIND ISPEC')
import ispec
def normalise_spectra(data,resolution=55000):
	model = "Splines" # "Polynomy"
	degree = 2
	nknots = None # Automatic: 1 spline every 5 nm

	# Strategy: Filter first median values and secondly MAXIMUMs in order to find the continuum
	order='median+max'
	median_wave_range=0.05
	max_wave_range=1.0

	strong_lines = ispec.read_line_regions(ispec_dir + "/input/regions/strong_lines/absorption_lines.txt")
	star_continuum_model = ispec.fit_continuum(data, from_resolution=resolution, \
		                ignore=strong_lines, \
		                nknots=nknots, degree=degree, \
		                median_wave_range=median_wave_range, \
		                max_wave_range=max_wave_range, \
		                model=model, order=order, \
		                automatic_strong_line_detection=True, \
		                strong_line_probability=0.5, \
		                use_errors_for_fitting=True)

	data = ispec.normalize_spectrum(data, star_continuum_model, consider_continuum_errors=False)
	return data
