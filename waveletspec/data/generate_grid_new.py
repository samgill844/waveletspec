from __future__ import print_function
import sys
import cPickle as pickle
import numpy as np
import logging
import multiprocessing as mp
from multiprocessing import Pool
import os,sys
import shelve
import mlpy.wavelet as w
from scipy.interpolate import interp1d
from astropy.io import fits
import glob
LOG_LEVEL = "critical"
logger = logging.getLogger() # root logger, common for all
logger.setLevel(logging.getLevelName(LOG_LEVEL.upper()))
import numpy as np

ispec_dir= '/home/sam/iSpec_v20160930/' # this needs to be where the package is installed
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec


###############
# Documentation
###############
# Synthesise a grid for the wavelet analysis procedure
# This has only been tested as by being run directly and not invoked through the module.
#	e.g. $ python generate_grid.py
#	>> follow instructions.

############################
# Parameters for synthesis
############################

# All we need is coordinates i.e.
for code in ['sme']:
        model = ispec_dir + "/input/atmospheres/MARCS.GES/modeled_layers_pack.dump"
	modeled_layers_pack = ispec.load_modeled_layers_pack(model)
	modeled_layers, used_values_for_layers, proximity, spherical, Trange, Lrange, Mrange, nlayers = modeled_layers_pack  

	Trange = np.arange(4000,8000,250)
	Vrange = [0]
	Lrange = np.arange(3.5,5.1,0.25)
	Mrange = np.arange(-1,1.1,0.25)
	limb_darkening_coeff = 0.5 
	# We could just as easily use np.arange or np.linspace to sdo the same thing...

	# Print number of spectra and approx memory
	print('\nWith the range  you have selected there will be a total of '+str(len(Trange)*len(Mrange)*len(Lrange)*len(Vrange))+' spectra synthesised. ')
	print('I estimate this to be '+str(0.001*len(Trange)*len(Mrange)*len(Lrange)*len(Vrange))+' Gb of storage for the files and is also needed in RAM to collect the files into a FITS file to hold the grid. The estimated time will be around '+str(len(Trange)*len(Mrange)*len(Lrange)*len(Vrange)/60)+' Hours \n\nYU HAVE BEEN WARNED!')


	######################
	# Main synthesis code
	######################


	resolution = None
	wave_step = 0.001
	regions = None
	#wave_base = 515.0 # Magnesium triplet region
	#wave_top = 525.0

	wave_base = 449.9 #wasp
	wave_top = 700.9


	# Selected model amtosphere, linelist and solar abundances
	#model = ispec_dir + "/input/atmospheres/MARCS/modeled_layers_pack.dump"
	model = ispec_dir + "/input/atmospheres/MARCS.GES/modeled_layers_pack.dump"
	#model = ispec_dir + "/input/atmospheres/MARCS.APOGEE/modeled_layers_pack.dump"
	#model = ispec_dir + "/input/atmospheres/ATLAS9.APOGEE/modeled_layers_pack.dump"
	#model = ispec_dir + "/input/atmospheres/ATLAS9.Castelli/modeled_layers_pack.dump"
	#model = ispec_dir + "/input/atmospheres/ATLAS9.Kurucz/modeled_layers_pack.dump"
	#model = ispec_dir + "/input/atmospheres/ATLAS9.Kirby/modeled_layers_pack.dump"

	atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.300_1100nm/atomic_lines.tsv"
	#atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.1100_2400nm/atomic_lines.tsv"
	#atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv5_atom_hfs_iso.420_920nm/atomic_lines.tsv"
	#atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv5_atom_nohfs_noiso.420_920nm/atomic_lines.tsv"

	isotope_file = ispec_dir + "/input/isotopes/SPECTRUM.lst"

	# Load chemical information and linelist
	atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, wave_base=wave_base, wave_top=wave_top)
	atomic_linelist = atomic_linelist[atomic_linelist['theoretical_depth'] >= 0.01] # Select lines that have some minimal contribution in the sun

	isotopes = ispec.read_isotope_data(isotope_file)

	#solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"
	#solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2005/stdatom.dat"
	solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2009/stdatom.dat"
	#solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
	#solar_abundances_file = ispec_dir + "/input/abundances/Anders.1989/stdatom.dat"

	# Load model atmospheres
	modeled_layers_pack = ispec.load_modeled_layers_pack(model)
	# Load SPECTRUM abundances
	fixed_abundances = None # No fixed abundances
	solar_abundances = ispec.read_solar_abundances(solar_abundances_file)
	ref = np.linspace(450.0,650.0,131072)




	def synthesize_spectrum(PARAMETER):
	    curr_proc=mp.current_process() 
	    curr_proc.daemon=False 
	    turbo=False


	    teff = PARAMETER[0]
	    vsini = PARAMETER[1]
	    MH = PARAMETER[2]
	    logg = PARAMETER[3]

	    #microturbulence_vel = 0.89 + (4.16*(teff-5777)*(10**-4))+(9.25*10**(-8))*((teff-5777)**2) # Km/s Amanda Doyles Thesis

	    #macroturbulence = 3.21 + 2.33*(teff-5777)*(10**-3) + 2.00*(((teff-5777)**2)*10**(-6))-2.00*(logg-4.44) # Km/s Amanda vmic calibration

	     # Validate parameters
	    if not ispec.valid_atmosphere_target(modeled_layers_pack, teff, logg, MH):
		msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] \
		        fall out of theatmospheric models."
		print(msg)
#		continue

	    # Enhance alpha elements + CNO abundances following MARCS standard composition
	    alpha_enhancement, c_enhancement, n_enhancement, o_enhancement = ispec.determine_abundance_enchancements(MH)
	    abundances = ispec.enhance_solar_abundances(solar_abundances, alpha_enhancement, c_enhancement, n_enhancement, o_enhancement)

	    # Prepare atmosphere model
	    atmosphere_layers = ispec.interpolate_atmosphere_layers(modeled_layers_pack, teff, logg, MH, code=code)

	    # Synthesis
	    synth_spectrum = ispec.create_spectrum_structure(np.arange(wave_base, wave_top, wave_step))
	    synth_spectrum['flux'] = ispec.generate_spectrum(synth_spectrum['waveobs'], \
		    atmosphere_layers, teff, logg, MH, atomic_linelist, isotopes, abundances, \
		    fixed_abundances, microturbulence_vel = 0., \
		    macroturbulence=0., vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, \
		    R=resolution, regions=regions, verbose=1,
		    code=code)

	    f = interp1d(synth_spectrum['waveobs'], synth_spectrum['flux'], kind='linear') # Interpolate Data for grid
	    ref = np.linspace(450,700,2**18)
	    tmp = np.array([f(ref),teff,MH,logg,vsini])
	    np.save(str(teff)+'_'+str(vsini)+'_'+str(MH)+'_'+str(logg),tmp)


	PARAMETERS=[]
	for i in range(0,len(Trange)):
		for j in range(0,len(Vrange)):
			for k in range(0,len(Mrange)):
				for l in range(0,len(Lrange)):
					PARAMETERS.append([Trange[i],Vrange[j],Mrange[k],Lrange[l]])
			


	print('Beggining the synthesis: ')

	pool = Pool(processes=8)
	pool.map(synthesize_spectrum,PARAMETERS)
	pool.close()
	pool.join()



	files = glob.glob('*.npy')
	spectra =[] # np2d array for primary HDU
	waveobs = np.linspace(450,700,2**18) # wavelength scale for spectra - leave
	# table params
	teff=[]
	mh=[]
	logg=[]
	Vsin=[]

	print('\n\nCreating arrays')
	for i in range(0,len(files)): # populate spectra and table parameters
		data = np.load(files[i])
		spectra.append(data[0]) # populate spectra

		teff.append(data[1]) # temperature
		mh.append(data[2]) #not botherd over vsini, skip to mh
		logg.append(data[3]) # finish with logg
		Vsin.append(data[4]) # added vsini for quick grid
		print(str(i*100./len(files))+'%')

	# now should have 5 arrays ready for fits storage.

	# FITS section
	print('array building complete')
	print('building fits')
	hdu = fits.PrimaryHDU(spectra) # set primary hdu for spectra info

	col1 = fits.Column(name='Temperature', format='E', array=teff)
	col2 = fits.Column(name='MH', format='E', array=mh)
	col3 = fits.Column(name='Logg', format='E', array=logg)
	col4 = fits.Column(name='Vsini', format='E', array = Vsin)
	cols = fits.ColDefs([col1, col2,col3,col4])
	tbhdu = fits.BinTableHDU.from_columns(cols)

	hdulist = fits.HDUList([hdu,tbhdu])
	hdulist.writeto('sme.fits')
	print('Fits building complete')
#	for i in files:
#		os.remove(i)















