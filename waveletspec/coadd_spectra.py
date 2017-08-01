from __future__ import print_function
import os,sys
import numpy as np
from astropy.io import fits
import waveletspec,os,sys
import matplotlib.pyplot as plt
this_dir, this_filename = os.path.split(waveletspec.__file__)
import glob

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
        print('Found python 3 version of iSpec')
        ispec_dir = this_dir+'/iSpec_v20160930_py3'
        sys.path.insert(0, os.path.abspath((this_dir+'/iSpec_v20160930_py3')))
except OSError as e:
    # export the tar
    print('I cant find the iSpec directory, searching for the tar files...')
    if (sys.version_info > (3, 0)):
        os.system('tar -xf {}/iSpec_v20160930_py3.tar.gz -C {}'.format(this_dir,this_dir))
        print('Exctracted to {}/iSpec_v20160930_py3'.format(this_dir))
        sys.path.insert(0, os.path.abspath((this_dir+'/iSpec_v20160930_py3')))
        ispec_dir = this_dir+'/iSpec_v20160930_py3'
        
    elif (sys.version_info < (3, 0)):
        os.system('tar -xf {}/iSpec_v20160930_py2.tar.gz -C {}'.format(this_dir,this_dir))
        print('Exctracted to {}/iSpec_v20160930_py2'.format(this_dir))
        sys.path.insert(0, os.path.abspath((this_dir+'/iSpec_v20160930_py2')))
        ispec_dir = this_dir+'/iSpec_v20160930_py2'
    else:
        print('I CANT FIND ISPEC')
import ispec

from progress.bar import Bar


def determine_radial_velocity_with_template(data,velocity_step):
    # - Read synthetic template
    #template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/Atlas.Arcturus.372_926nm/template.txt.gz")
    #template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/Atlas.Sun.372_926nm/template.txt.gz")
    template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/NARVAL.Sun.370_1048nm/template.txt.gz")
    #template = ispec.read_spectrum(ispec_dir + "/input/spectra/templates/Synth.Sun.300_1100nm/template.txt.gz")

    models, ccf = ispec.cross_correlate_with_template(data, template, \
                            lower_velocity_limit=-200, upper_velocity_limit=200, \
                            velocity_step=velocity_step, fourier=False)

    # Number of models represent the number of components
    components = len(models)
    print(components)
    # First component:
    rv = np.round(models[0].mu(), 2) # km/s
    rv_err = np.round(models[0].emu(), 2) # km/s
    return rv, rv_err


def coadd_spectra(files=None,scale_method = 'nothing', combine_method='coadd', velocity_step=1.0):
	########################
	# First sort files out
	########################
	if files==None:
		files = glob.glob('COR*.fits')
	if len(files) == 1:
		print('Use the "load_spectra" command for a single file.')
		return
	if len(files) == 0:
		print('No files found.')
		return

	
	##############################
	# Now cycle through the files
	##############################
	spectra = []
	wavelengths = []
	
	bar = Bar('Correcting for RVs', max=len(files))
	for i in range(len(files)):
		data = ispec.read_spectrum(files[i]) 
		#data = ispec. resample_spectrum(data, np.linspace(627.5, 680.5,2**14)) # FOR INT SPECTRA



		rv, rv_err = determine_radial_velocity_with_template(data,velocity_step)
		print('{}: {} +- {}'.format(files[i],rv,rv_err))
		#if abs(rv)>100:
		#	bar.next()
		#	continue
		data = ispec.correct_velocity(data, rv)
		spectra.append(data['flux'])
		wavelengths.append(data['waveobs'])
		bar.next()
	bar.finish()

	if len(spectra)==0:
		print('RV correction produced no spectra to co-add.')
		return

	##################################################################
	# Now find the max of the mins and mins of the max for wavelength
	##################################################################
	mins,maxs =[],[]
	for i in range(len(wavelengths)):
		mins.append(min(wavelengths[i]))
		maxs.append(max(wavelengths[i]))
	minimum_wavelength = max(mins)+1
	maximum_wavelength = min(maxs)-1
	print('The minimum usable wavelength range is {} - {} nm'.format(minimum_wavelength,maximum_wavelength))


	#################################
	# Now re-interpolate the spectra
	#################################

	wavelength_range = np.linspace(minimum_wavelength,maximum_wavelength,2**18)
	bar = Bar('Re-interpolating spectra', max=len(wavelengths))
	for i in range(len(wavelengths)):
		spectra[i] = np.interp(wavelength_range,wavelengths[i], spectra[i])
		bar.next()
	bar.finish()
	del wavelengths # keep useage low
	spectra = np.array(spectra) # for ease later


	#####################################
	# Now use method to scale spectra
	#####################################
	if scale_method == 'median':
		medians = np.median(spectra,axis=1)
		for i in range(len(spectra)):
			spectra[i] = spectra[i]/medians[i]
		print('Scale method: median')
	elif scale_method == 'mean':
		means = np.mean(spectra,axis=1)
		for i in range(len(spectra)):
			spectra[i] = spectra[i]/means[i]
		print('Scale method: mean')
	elif scale_method == 'nothing':
		print('Scale method: nothing')
		pass
	else:
		print('Scale method no understood')
		return

	
	#####################################################
	# Now plot spectra to ensure the method is behaving
	#####################################################
	fig1 = plt.figure()
	ax1 = fig1.add_subplot(111)
	number_of_plots = len(spectra)
	colormap = plt.cm.nipy_spectral #I suggest to use nipy_spectral, Set1,Paired
	ax1.set_color_cycle([colormap(i) for i in np.linspace(0, 1,number_of_plots)])
	
	for i in range(len(spectra)):
		plt.plot(wavelength_range,spectra[i],label=files[i])
	plt.legend()
	plt.xlabel('Wavelength (nm)')
	plt.ylabel('Counts')
	#plt.xlim(683.9,739.5)
	plt.show(block=False)

	response = input('Are you happy? ')
	if (response.lower() == 'y') or (response.lower() == 'yes'):
		pass
	else:
		print('Breaking')
		return


	#####################################
	# Now use method to combine spectra
	#####################################
	if combine_method == 'median':
		spectrae = np.std(spectra,axis=0)
		spectra = np.median(spectra,axis=0)
		print('Combination method: median')
	elif combine_method == 'mean':
		spectrae = np.std(spectra,axis=0)
		spectra = np.mean(spectra,axis=0)
		print('Combination method: mean')
	elif combine_method == 'coadd':
		spectra = np.sum(spectra,axis=0)
		spectrae = [0.0]*len(spectra)
		print('Combination method: coadd')
		pass
	else:
		print('Combinaton method no understood')
		return
		
	plt.plot(wavelength_range,spectra,'r--',linewidth=1.5)
	response = input('Check the red line. Happy?  ')
	if (response.lower() == 'y') or (response.lower() == 'yes'):
		pass
	else:
		print('Breaking')
		return

	data = ispec.create_spectrum_structure(wavelength_range)
	data['waveobs'] = wavelength_range
	data['flux'] = spectra
	data['err'] = spectrae

	return data

	

	
