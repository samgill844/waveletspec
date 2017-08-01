from __future__ import print_function
import os,sys

from astropy.io import fits
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

import numpy as np
from astropy.io import fits
import emcee
import corner
from progress.bar import Bar
################################################################################################################
# SECTION 0 - IMPORTS 
################################################################################################################
import logging
import multiprocessing as mp
from multiprocessing import Pool
import mlpy.wavelet as w
from scipy.interpolate import interp1d
import glob
LOG_LEVEL = "critical"
logger = logging.getLogger() # root logger, common for all
logger.setLevel(logging.getLevelName(LOG_LEVEL.upper()))
from heapq import nsmallest
import time
import scipy.stats as s
from scipy.misc import logsumexp
import math
from scipy.stats import norm
from scipy.signal import fftconvolve
from multiprocessing import Process
from astropy.stats import mad_std
import matplotlib
from astropy.table import hstack
from matplotlib.patches import Rectangle
import subprocess
from astropy.table import Table, Column
import itertools
from progress.bar import Bar
from scipy import signal
import matplotlib.pyplot as plt
################################################################################################################
# SECTION 1 - FUNCTIONS CALL OPTIMISED FOR MULTICORE PROCESSING
################################################################################################################
bar_vals = {}
class Progressbar(Bar):
    suffix = '%(parameters)s %(timeleft)s  %(index)d / %(max)d'	
    @property
    def parameters(self):
	    t,m,l,v,j = np.mean(bar_vals['step'],axis=0)
	    dt,dm,dl,dv,dj = np.std(bar_vals['step'],axis=0)
	    return '{:.0f} +- {:.0f}    {:.2f} +- {:.2f}    {:.2f} +- {:.2f}     {:.2f} +- {:.2f}   '.format(t,dt,m,dm,l,dl,v,dv)

    @property
    def timeleft(self):
        return bar_vals['timeleft']




def broaden_all_fast(wavelength,flux,vmac=2,vsini=2,resolution=55000):
    #####################################
    # convert to a uniform velocity space
    #####################################
    tic_b_ = time.time()
    lnwave_ = np.linspace(np.log(wavelength[0]),np.log(wavelength[-1]),len(wavelength))
    flnw_ = np.interp(lnwave_,np.log(wavelength),flux) # Flux on uniform log wave 
    c_kms =  299792.458
    velocity_step = (np.exp(lnwave_[1]-lnwave_[0])-1)*c_kms

    # define the kernals
    #vmac_kernel,vsini_kernel,resolution_kernel = [],[],[]
    #print('\t- Logspace: {}'.format(time.time()-	tic_b_))
    #tic_b = time.time()


    kernals_to_use=[]

    if (vmac <=0) and (vsini <= 0) and (resolution <= 0):
        #print('No kernals used')
        return wavelength,flux

	#########################
	# Construct vmac kernel
	#########################
    if vmac > 0:
        # mu represent angles that divide the star into equal area annuli,
        # ordered from disk center (mu=1) to the limb (mu=0).
        # But since we don't have intensity profiles at various viewing (mu) angles
        # at this point, we just take a middle point:
        m = 0.5
        # Calc projected simga for radial and tangential velocity distributions.
        sigma = vmac/np.sqrt(2.0) / velocity_step
        sigr = sigma * m
        sigt = sigma * np.sqrt(1.0 - m**2.)
        # Figure out how many points to use in macroturbulence kernel
        nmk = max(min(round(sigma*10), (len(flux)-3)/2), 3)
        # Construct radial macroturbulence kernel w/ sigma of mu*vmac/sqrt(2)
        if sigr > 0:
            xarg = (np.arange(2*nmk+1)-nmk) / sigr   # exponential arg
            #mrkern = np.exp(max((-0.5*(xarg**2)),-20.0))
            mrkern = np.exp(-0.5*(xarg**2))
            mrkern = mrkern/mrkern.sum()
        else:
            mrkern = np.zeros(2*nmk+1)
            mrkern[nmk] = 1.0    #delta function

        # Construct tangential kernel w/ sigma of sqrt(1-mu**2)*vmac/sqrt(2.)
        if sigt > 0:
            xarg = (np.arange(2*nmk+1)-nmk) /sigt
            mtkern = np.exp(-0.5*(xarg**2))
            mtkern = mtkern/mtkern.sum()
        else:
            mtkern = np.zeros(2*nmk+1)
            mtkern[nmk] = 1.0
        ## Sum the radial and tangential components, weighted by surface area
        area_r = 0.5
        area_t = 0.5
        vmac_kernel = area_r*mrkern + area_t*mtkern
        kernals_to_use.append(vmac_kernel)
    else:
	    pass

    #print('\t- Vmac: {}'.format(time.time()-tic_b))
    #tic_b = time.time()
    #########################
    # Construct vsini kernel
    #########################

    if vsini > 0:
        epsilon=0.6
        #wave_ = np.log(wavelength)
        #velo_ = np.linspace(wave_[0],wave_[-1],len(wave_))
        #flux_ = np.interp(lnwave_,wave_,flux)
        #dvelo = velo_[1]-velo_[0]
        dvelo = velocity_step/c_kms
        vrot = vsini/(299792458.*1e-3) # in Km/s
        #-- compute the convolution kernel and normalise it
        # PFLM/SG - changed calculation of kernel to use an odd number of pixels so that 
        # there is at least one pixel in the kernel - more reliable for low vrot
        n = int(1+2*vrot/dvelo) # always have at least one pixel in the kernel !
        #n = int(3+2*vrot/dvelo) # tried at least 3 elements - it fails
        velo_k = np.arange(0,n)*dvelo
        velo_k -= velo_k[-1]/2.
        y = 1 - (velo_k/vrot)**2 # transformation of velocity
        vsini_kernel = (2*(1-epsilon)*np.sqrt(y)+np.pi*epsilon/2.*y)/(np.pi*vrot*(1-epsilon/3.0))  # the kernel
        vsini_kernel /= vsini_kernel.sum()
        kernals_to_use.append(vsini_kernel)
    else:
	    pass
	
    #vsini_kernel = [1]
    #print('\t- Vsini: {}'.format(time.time()-tic_b))
    #tic_b = time.time()
    ##############################
    # Construct resolution kernel
    ##############################
    if resolution > 0:
	    #lnwave_ = np.linspace(np.log(wavelength[0]),np.log(wavelength[-1]),len(wavelength))
	    #flnw_ = np.interp(lnwave_,np.log(wavelength),flux_spec) # Flux on uniform log wave 
	    #dvelo = (np.exp(lnwave_[1]-lnwave_[0])-1)*c_kms
	    v_fwhm =c_kms/resolution # in km/s
	    n = int(1+4*v_fwhm/velocity_step) # always have at least one pixel in the kernel !
	    velo_k = np.linspace(-n*velocity_step,n*velocity_step,2*n + 1)
	    v_sigma = v_fwhm/(2*np.sqrt(2*np.log(2)))
	    resolution_kernel =  np.exp(-0.5*velo_k**2/(v_sigma**2)) # the kernel
	    resolution_kernel /= resolution_kernel.sum()
	    kernals_to_use.append(resolution_kernel)

    else:
	    pass
    #print('\t- Resolution: {}'.format(time.time()-tic_b))
    #tic_b = time.time()

    #OLD WAY
    '''
    Interpolation: 0.0085289478302
	    - Logspace: 0.00404000282288
	    - Vmac: 0.000119924545288
	    - Vsini: 6.89029693604e-05
	    - Resolution: 7.10487365723e-05
	    - kernel combining: 0.000430107116699
	    - Actual broadening: 0.0108051300049
	    - Interpolating: 0.00414109230042
    Broadening: 0.0197830200195
    Re-interplating: 0.00380778312683
    Dicting: 0.0040590763092
    Total: 0.0362620353699


    '''
    #print('vmac kernel')
    #print(vmac_kernel)

    #print('vsini kernel')
    #print(vsini_kernel)

    #print('resolution kernel')
    #print(resolution_kernel)

    ##########################
    # Create the master kernel
    ##########################
    #print('Number of kernels: {}'.format(len(kernals_to_use)))
    #print(kernals_to_use)

    kernals_to_use = sorted(kernals_to_use, key=len)
    kernals_to_use = kernals_to_use[::-1] # list of kernels - longest first

    master_kernel = kernals_to_use[0]
    for i in range(1,len(kernals_to_use)):
        master_kernel = signal.convolve(master_kernel,kernals_to_use[i],mode='same')

    #master_kernel = signal.convolve(vmac_kernel,vsini_kernel,mode='same')
    #master_kernel = fftconvolve(vsini_kernel,vmac_kernel,mode='same')
    #print(master_kernel)
    #master_kernel = signal.convolve(master_kernel,resolution_kernel,mode='same')
    #master_kernel = fftconvolve(master_kernel,resolution_kernel,mode='same')




    #print(master_kernel)
    #print('\t- kernel combining: {}'.format(time.time()-tic_b))
    #tic_b = time.time()
    ###############################
    # Now do the final convolution
    ###############################
    #flnw_conv = 1 - signal.convolve(1-flnw_, master_kernel, mode='same') # Fastest
    flnw_conv = 1 - fftconvolve(1-flnw_, master_kernel, mode='same') # Fastest
    #print('\t- Actual broadening: {}'.format(time.time()-tic_b))
    tic_b = time.time()
    flux_conv = np.interp(wavelength,np.exp(lnwave_),flnw_conv)
    #print('\t- Interpolating: {}'.format(time.time()-tic_b))
    return wavelength, flux_conv

def broaden_all_fastn(wavelength,flux,vmac=2,vsini=2,resolution=55000):
    #####################################
    # convert to a uniform velocity space
    #####################################
    tic_b_ = time.time()
    lnwave_ = np.linspace(np.log(wavelength[0]),np.log(wavelength[-1]),len(wavelength))
    flnw_ = np.interp(lnwave_,np.log(wavelength),flux) # Flux on uniform log wave 
    c_kms =  299792.458
    velocity_step = (np.exp(lnwave_[1]-lnwave_[0])-1)*c_kms

    # define the kernals
    vmac_kernel,vsini_kernel,resolution_kernel = [],[],[]
    #print('\t- Logspace: {}'.format(time.time()-	tic_b_))
    tic_b = time.time()


    kernals_to_use=[]

    if (vmac <=0) and (vsini <= 0) and (resolution <= 0):
	    print('No kernals used')
	    return wavelength,flux

	#########################
	# Construct vmac kernel
	#########################
    if vmac > 0:
	    # mu represent angles that divide the star into equal area annuli,
	    # ordered from disk center (mu=1) to the limb (mu=0).
	    # But since we don't have intensity profiles at various viewing (mu) angles
	    # at this point, we just take a middle point:
	    m = 0.5
	    # Calc projected simga for radial and tangential velocity distributions.
	    sigma = vmac/np.sqrt(2.0) / velocity_step
	    sigr = sigma * m
	    sigt = sigma * np.sqrt(1.0 - m**2.)
	    # Figure out how many points to use in macroturbulence kernel
	    nmk = max(min(round(sigma*10), (len(flux)-3)/2), 3)
	    # Construct radial macroturbulence kernel w/ sigma of mu*vmac/sqrt(2)
	    if sigr > 0:
	        xarg = (np.arange(2*nmk+1)-nmk) / sigr   # exponential arg
	        #mrkern = np.exp(max((-0.5*(xarg**2)),-20.0))
	        mrkern = np.exp(-0.5*(xarg**2))
	        mrkern = mrkern/mrkern.sum()
	    else:
	        mrkern = np.zeros(2*nmk+1)
	        mrkern[nmk] = 1.0    #delta function

	    # Construct tangential kernel w/ sigma of sqrt(1-mu**2)*vmac/sqrt(2.)
	    if sigt > 0:
	        xarg = (np.arange(2*nmk+1)-nmk) /sigt
	        mtkern = np.exp(-0.5*(xarg**2))
	        mtkern = mtkern/mtkern.sum()
	    else:
	        mtkern = np.zeros(2*nmk+1)
	        mtkern[nmk] = 1.0

	    ## Sum the radial and tangential components, weighted by surface area
	    area_r = 0.5
	    area_t = 0.5
	    vmac_kernel = area_r*mrkern + area_t*mtkern
	    kernals_to_use.append(vmac_kernel)

    else:
	    pass

    #print('\t- Vmac: {}'.format(time.time()-tic_b))
    tic_b = time.time()
    #########################
    # Construct vsini kernel
    #########################

    if vsini > 0:
        epsilon=0.6
        #wave_ = np.log(wavelength)
        #velo_ = np.linspace(wave_[0],wave_[-1],len(wave_))
        #flux_ = np.interp(lnwave_,wave_,flux)
        #dvelo = velo_[1]-velo_[0]
        dvelo = velocity_step/c_kms
        vrot = vsini/(299792458.*1e-3) # in Km/s
        #-- compute the convolution kernel and normalise it
        # PFLM/SG - changed calculation of kernel to use an odd number of pixels so that 
        # there is at least one pixel in the kernel - more reliable for low vrot
        n = int(1+2*vrot/dvelo) # always have at least one pixel in the kernel !
        #n = int(3+2*vrot/dvelo) # tried at least 3 elements - it fails
        velo_k = np.arange(0,n)*dvelo
        velo_k -= velo_k[-1]/2.
        y = 1 - (velo_k/vrot)**2 # transformation of velocity
        vsini_kernel = (2*(1-epsilon)*np.sqrt(y)+np.pi*epsilon/2.*y)/(np.pi*vrot*(1-epsilon/3.0))  # the kernel
        vsini_kernel /= vsini_kernel.sum()
        kernals_to_use.append(vsini_kernel)
    else:
	    pass

    #vsini_kernel = [1]
    #print('\t- Vsini: {}'.format(time.time()-tic_b))
    tic_b = time.time()
    ##############################
    # Construct resolution kernel
    ##############################
    if resolution > 0:
        #lnwave_ = np.linspace(np.log(wavelength[0]),np.log(wavelength[-1]),len(wavelength))
        #flnw_ = np.interp(lnwave_,np.log(wavelength),flux_spec) # Flux on uniform log wave 
        #dvelo = (np.exp(lnwave_[1]-lnwave_[0])-1)*c_kms
        v_fwhm =c_kms/resolution # in km/s
        n = int(1+4*v_fwhm/velocity_step) # always have at least one pixel in the kernel !
        velo_k = np.linspace(-n*velocity_step,n*velocity_step,2*n + 1)
        v_sigma = v_fwhm/(2*np.sqrt(2*np.log(2)))
        resolution_kernel =  np.exp(-0.5*velo_k**2/(v_sigma**2)) # the kernel
        resolution_kernel /= resolution_kernel.sum()
        kernals_to_use.append(resolution_kernel)
    else:
	    pass
    #print('\t- Resolution: {}'.format(time.time()-tic_b))
    tic_b = time.time()

    #OLD WAY
    '''
    Interpolation: 0.0085289478302
	    - Logspace: 0.00404000282288
	    - Vmac: 0.000119924545288
	    - Vsini: 6.89029693604e-05
	    - Resolution: 7.10487365723e-05
	    - kernel combining: 0.000430107116699
	    - Actual broadening: 0.0108051300049
	    - Interpolating: 0.00414109230042
    Broadening: 0.0197830200195
    Re-interplating: 0.00380778312683
    Dicting: 0.0040590763092
    Total: 0.0362620353699


    '''
    #print('vmac kernel')
    #print(vmac_kernel)

    #print('vsini kernel')
    #print(vsini_kernel)

    #print('resolution kernel')
    #print(resolution_kernel)

    ##########################
    # Create the master kernel
    ##########################
    #print('Number of kernels: {}'.format(len(kernals_to_use)))
    #print(kernals_to_use)

    kernals_to_use = sorted(kernals_to_use, key=len)
    kernals_to_use = kernals_to_use[::-1] # list of kernels - longest first

    master_kernel = kernals_to_use[0]
    for i in range(1,len(kernals_to_use)):
	    master_kernel = signal.convolve(master_kernel,kernals_to_use[i],mode='same')
	
    #master_kernel = signal.convolve(vmac_kernel,vsini_kernel,mode='same')
    #master_kernel = fftconvolve(vsini_kernel,vmac_kernel,mode='same')
    #print(master_kernel)
    #master_kernel = signal.convolve(master_kernel,resolution_kernel,mode='same')
    #master_kernel = fftconvolve(master_kernel,resolution_kernel,mode='same')




    #print(master_kernel)
    #print('\t- kernel combining: {}'.format(time.time()-tic_b))
    tic_b = time.time()
    ###############################
    # Now do the final convolution
    ###############################
    #flnw_conv = 1 - signal.convolve(1-flnw_, master_kernel, mode='same') # Fastest
    flnw_conv = 1 - fftconvolve(1-flnw_, master_kernel, mode='same') # Fastest
    #print('\t- Actual broadening: {}'.format(time.time()-tic_b))
    tic_b = time.time()
    flux_conv = np.interp(wavelength,np.exp(lnwave_),flnw_conv)
    #print('\t- Interpolating: {}'.format(time.time()-tic_b))
    return wavelength, flux_conv
	


def cut_spectrum_from_segments(data,regions_name):
    segments = ispec.read_segment_regions(os.path.join(this_dir, "bad_regions", regions_name))
    wfilter = ispec.create_wavelength_filter(data, regions=segments)
    return data[~wfilter]

def estimate_turbulence_from_doyle(teff,logg):
	microturbulence_vel =  0.89 + (4.16*(teff-5777)*(10**-4))+(9.25*10**(-8))*((teff-5777)**2) # Km/s Amanda Doyles Thesis
	macroturbulence =  3.21 + 2.33*(teff-5777)*(10**-3) + 2.00*(((teff-5777)**2)*10**(-6))-2.00*(logg-4.44) # Km/s Amanda vmic calibration

	#print('microturbulence_vel =  {}, macroturbulence =  {}'.format(microturbulence_vel,macroturbulence))
	return microturbulence_vel, macroturbulence

def wavelet_analysis(data, method = 'single_walker+emcee',wavemin=452,wavemax=648,n_use=17,n_low=4,n_high=14,Nweight=10,T=None,M=None,L=None,V=None,Prior_values='flat',draws=1000,dwt_type=bytes('d', 'utf-8'),dwt_number=4,save_data = True,grid_name='spectrum',resolution=55000,chain_file = 'chain.fits',threads=1,silent=False, cut_regions=None,email_notification = True):
    #####################
    # Check within range
    #####################
    if wavemin<min(data['waveobs']) or wavemax>max(data['waveobs']):
	    print('The standard range (452nm-648nm) is out of range of your data.\nYour current range is '+str(min(data['waveobs']))+'nm - '+str(max(data['waveobs']))+'nm\nPlease specify a range with the wavemin and wavemax argument')
	    return

    ###################
    # Set up priors
    ###################
    if Prior_values=='flat':
	    if not silent:
		    print('Using Flat Priors.\n')
	    Prior_values = np.array([False]*10)


    ###################################################
    # Set up log likelihood functions from each method
    ###################################################

    try:
	    rtc = grid_name+'.fits'
	    # path_to_grid = raw_input('Please specify path to grid: ')
	    path_to_grid = os.path.join(this_dir, "Grids", rtc)
	    h = fits.open(path_to_grid)
	    hdulist  =fits.open(path_to_grid)
	    if not silent:
		    print('Using grid from:')
		    print(path_to_grid)
    except:
	    print('Grid not found. change in wavelet_analysis file.')
	    print(path_to_grid)
	    return data

    def lnlike(theta,home,temperature,metalicity,loggs,real,bottom_index,top_index,err,wavemin,wavemax,n_use,hdulist,Prior_values,dwt_type,dwt_number,resolution,data,cut_regions=None):
	    #tic = time.time()
	    teff,mh,logg,vrot,beta = theta
	    ##################
	    #interpolate model
	    ##################
	    teff_2_closest = nsmallest(2, list(set(temperature)), key=lambda x: abs(x-teff)); teff_2_closest.sort() #G
	    mh_2_closest = nsmallest(2, list(set(metalicity)), key=lambda x: abs(x-mh));mh_2_closest.sort() #G
	    logg_2_closest = nsmallest(2, list(set(loggs)), key=lambda x: abs(x-logg));logg_2_closest.sort() #G
	    # finding spectra
	    models = [[i,j,k] for i in teff_2_closest for j in mh_2_closest for k in logg_2_closest]
	    models_index = [get_index(i,temperature,metalicity,loggs) for i in models]
	    T,M,L = [i[0] for i in models],[i[1] for i in models],[i[2] for i in models] # temperature coordinates
	    # set values to interpolate with and the distance to the box
	    t = teff-min(T) # x
	    tstep = max(T)-min(T)
	    m = mh-min(M) #y
	    mstep = max(M)-min(M)
	    l=logg-min(L) # z
	    lstep=max(L)-min(L)
	    models1 = models
	    models =[hdulist[0].data[i] for i in models_index]
	    # main interp (Vijk in example)
	    interp_spectra = np.array(models[0])*(tstep-t)*(mstep-m)*(lstep-l) \
	    + np.array(models[4])*t*(mstep-m)*(lstep-l)\
	    + np.array(models[2])*(tstep-t)*m*(lstep-l)\
	    + np.array(models[1])*(tstep-t)*(mstep-m)*l\
	    + np.array(models[5])*t*(mstep-m)*l \
	    +np.array(models[3])*(tstep-t)*m*l \
	    + np.array(models[6])*t*m*(lstep-l) \
	    + np.array(models[7])*t*l*m
	    model =  interp_spectra/(tstep*lstep*mstep)
	    if rtc=='spectrum.fits':
		    model_wave = np.linspace(350,800,2**18)
	    else:
		    raise ValueError('No grid found')

	    ##################
	    # Re in-terpolate
	    ##################
	    f = interp1d(model_wave, model)
	    model_wave = np.linspace(wavemin,wavemax,2**n_use)
	    model = f(model_wave)

	    ###################################
	    # Get vmac values from calibration
	    ###################################
	    '''
	    if teff > 6400:
		    macroturbulence = 6
	    #elif teff < 5200:
	    #	macroturbulence = 2
	    else:
		    macroturbulence = 3.21 + 2.33*(teff-5777)*(10**-3) + 2.00*(((teff-5777)**2)*10**(-6))-2.00*(logg-4.44) # Km/s Amanda vmic calibration

	    if macroturbulence>6:
		    macroturbulence = 6
	    '''
	    microturbulence_vel, macroturbulence = estimate_turbulence_from_doyle(teff,logg)
	    #macroturbulence = ispec.estimate_vmac(teff,logg,mh) # km/s
	    #print(macroturbulence)
	    ######################
	    # broaden model
	    ######################
	    #model_wave, model = vmac_broadening(model_wave,model, macroturbulence)
	    #model_wave, model = broaden_vsini(model_wave, model, vrot)
	    #model_wave, model = instrumentally_broaden(model_wave, model , resolution=resolution)
	    model_wave, model = broaden_all_fastn(model_wave,model,macroturbulence,vrot,resolution)
	    #plt.plot(data['waveobs'],data['flux'])
	    #plt.plot(model_wave, model,'r')
	    #plt.show()

	    ###################################
	    # Now cut line regions if needed
	    ###################################
	    '''
	    if cut_regions != None:
		    # First cut the spectrum
		    data = ispec.create_spectrum_structure(model_wave)
		    data['flux'] = model
		    data = cut_spectrum_from_segments(data, cut_regions)

		    # Now we have to re-interpolate in pixel space so it has 2**n
		    model = np.interp(range(2**n_use),range(len(data)),data['flux'])
		    #plt.plot(model,label='model')
		    #plt.plot(flux,label='flux')
		    #plt.show()
	    '''	

	    #################################
	    # Wavelet transform in log space
	    #################################
	    model = w.dwt(np.log10(model),dwt_type,dwt_number)

	    #######################################
	    # To add a log g prior from photometry
	    #######################################
	    # here is an list of wasp planets and logg
	    # wasp	P_prior	P_prior_error
	    # 4	4.45	0.029
	    # 5 	4.403	0.039
	    # 6	4.50	0.06

	    # Methodd1 - error matches that of DWT coeffs
	    wt = 1./(np.power(err[bottom_index:top_index],2)*beta) # must be real
	    #print('SIGMA: {}'.format(sig))
	    # method 2
	    #wt = 1./(np.power(err[bottom_index:top_index],2) + np.power(sig,2))
	
	    #plt.semilogy(real[bottom_index:top_index],'b')
	    #plt.semilogy(model[bottom_index:top_index],'r')
	    #plt.show()
	    log_like_data = 0.5*(np.sum(np.log(wt)))  -0.5*np.sum(wt*np.power(real[bottom_index:top_index]-model[bottom_index:top_index],2))
	    #print('reduceschi_from_data: {}'.format(chi_from_data/len(err[bottom_index:top_index])))
	    #########################################
	    # Add prior Values (remember - negative) teff,mh,logg,vrot,beta
	    #########################################
	    chi_T = 0
	    if Prior_values[0]!=False and Prior_values[1]!=False:
		    chi_T = np.power(teff-Prior_values[0],2)/np.power(Prior_values[1],2)

	    chi_M = 0
	    if Prior_values[2]!=False and Prior_values[3]!=False:
		    chi_M = np.power(mh-Prior_values[2],2)/np.power(Prior_values[3],2)

	    chi_L = 0
	    if Prior_values[4]!=False and Prior_values[5]!=False:
		    chi_L = np.power(logg-Prior_values[4],2)/np.power(Prior_values[5],2)

	    chi_V = 0
	    if Prior_values[6]!=False and Prior_values[7]!=False:
		    chi_V = np.power(vrot-Prior_values[6],2)/np.power(Prior_values[7],2)
	
	    log_like_priors = -0.5*(chi_T +chi_M + chi_L + chi_V)

	    #print('log_lik priors: {}'.format(log_like_priors))
	    #print('log like data : {}'.format(log_like_data))
	    #print('weight : {}'.format(err))
	    # imperfect models atomic data, scales with solar metaliicity. 
	    # beta holds some sins.
	    ############################
	    # Return the log likly hood
	    ############################
	    #toc = time.time()
	    #print('Time: '+str(toc-tic))
	    #print(str(chi_from_data + chi_from_priors))
	    log_like = log_like_data + log_like_priors
	    #print('{:.0f} {:.4f} {:.4f} {:.4f} {:.2f} {:.2f}'.format(teff,mh,logg,vrot,beta,log_like))
	    #print('{} {} {} {} {} {}'.format(teff,mh,logg,vrot,sig,log_like))
	    return log_like

    ##############################
    # define the lnprior function
    ##############################	

    def lnprior(theta,rtc):
        # teff,mh,logg,vrot,sig = theta
        T,M,L,V,sig  = theta
        if (T<4000) or (T>8000):
	        #print('tef')
	        return -np.inf
        if (M<-1) or (M>1):
	        #print('MH')
	        return -np.inf
        #print(rtc)
        if rtc=='spectrum.fits':
	        if (L<3.5) or (L>5):
		        #print('L')
		        return -np.inf
        else:
	        print('FAILS')
	        return -np.inf
        if (V<0) or (V>100):
	        #print('Ve')
	        return -np.inf
        if (sig<0):
	        #print('sigVe')
	        return -np.inf

        #print(T,M,L,V,sig)
        return 0.0

    ##############################
    # define the lnprob function
    ##############################	
    def lnprob(theta,home,temperature,metalicity,loggs,real,bottom_index,top_index,err,wavemin,wavemax,n,hdulist,Prior_values,dwt_type,dwt_number,resolution,rtc,data,cut_regions):
        lp = lnprior(theta,rtc)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlike(theta,home,temperature,metalicity,loggs,real,bottom_index,top_index,err,wavemin,wavemax,n,hdulist,Prior_values,dwt_type,dwt_number,resolution,data,cut_regions)




    ##############################
    # now set up emcee ensembler #
    ##############################
    ndim = 5 # 4 free parameters
    chains = 2*ndim 
    # ...and a positive definite, non-trivial covariance matrix.
    cov  = 0.5-np.random.rand(ndim**2).reshape((ndim, ndim))
    cov  = np.triu(cov)
    cov += cov.T - np.diag(cov.diagonal())
    cov  = np.dot(cov,cov)
    # Invert the covariance matrix first.
    icov = np.linalg.inv(cov)
    # We'll sample with 250 walkers.





    ###################################
    # Now cut line regions if needed
    ###################################

    # Now we have to re-interpolate in pixel space so it has 2**n
    flux = np.interp(np.linspace(wavemin,wavemax,2**n_use),data['waveobs'],data['flux'])
    fluxe = np.interp(np.linspace(wavemin,wavemax,2**n_use),data['waveobs'],data['err'])

    ###########################################
    # Now get coeffs
    ############################################
    real = w.dwt(np.log10(flux),dwt_type,dwt_number)



    ###########
    # Weighting
    ###########
    if not silent:
	    print('Performing weighting... ',end='')
    weights =np.array([flux for i in range(Nweight)]) # pre-allocate
    if not silent:
	    bar = Bar('Performing weighting',max = Nweight)
    for i in range(Nweight):
	    try:
		    weights[i] = w.dwt( np.log10(100 + weights[i] + np.random.normal(0,fluxe)),dwt_type,dwt_number)
		    if not silent:
			    bar.next()
	    except:
		    weights[i] = w.dwt( np.log10(100 + weights[i] + np.random.normal(0,0.1)),dwt_type,dwt_number) # incase no flux err
		    if not silent:
			    bar.next()
    if not silent:
	    bar.finish()

    weight = weights.std(axis=0)
    del weights
    where_are_NaNs = np.isnan(weight)
    where_are_inf = np.isinf(weight)
    weight[where_are_NaNs]=0.003 # botch
    weight[weight<0.0005]=0.0005 # botch dwt_type='d',dwt_number=4
    bottom_index,top_index = (2**n_low) ,(2**n_high) - 1

    ########################
    # Preliminaries for MCMC
    ########################
    Trange,Mrange,Lrange = [],[],[]
    for i in range(0,len(hdulist[1].data)):
	    Trange.append(hdulist[1].data[i][0])
	    Mrange.append(hdulist[1].data[i][1])
	    Lrange.append(hdulist[1].data[i][2])
    Trange,Mrange,Lrange = sorted(list(set(Trange))),sorted(list(set(Mrange))),sorted(list(set(Lrange)))
    step_sizes = [ (Trange[-1]-Trange[0])/(len(Trange)-1),(Mrange[-1]-Mrange[0])/(len(Mrange)-1),(Lrange[-1]-Lrange[0])/(len(Lrange)-1),2] # Step size is originaly grid
    Trange = [min(Trange),max(Trange)]; Mrange = [min(Mrange),max(Mrange)]; Lrange = [min(Lrange),max(Lrange)]; #set inital range as whole grid
    temperature,metalicity,loggs = [i[0] for i in hdulist[1].data],[i[1] for i in hdulist[1].data],[i[2] for i in hdulist[1].data]; #put this outside - shave 0.6 seconds



    def get_index(i,T,M,L):
	    Tindicies = np.where(T==i[0])[0]
	    Mindicies = np.where(M==i[1])[0]
	    Lindicies = np.where(L==i[2])[0]
	    p =  list(set(Tindicies) & set(Mindicies) & set(Lindicies))
	    return p[0] 



    if method == 'emcee':
        print('EMCEE method')
        # Initialse Walkers roun solar value unless otherwis stated	model=0
        model=0
        sampler,pos=0,0

        chains = 100
        sampler = emcee.EnsembleSampler(chains, ndim, lnprob, args = [ispec_dir,temperature,metalicity,loggs,real,bottom_index,top_index,weight,wavemin,wavemax,n_use,hdulist,Prior_values,dwt_type,dwt_number,resolution,rtc,cut_regions,data], threads=threads)


        pos = []
        if not silent:
	        bar = Bar('Bulding starting positions',max = chains)
        while len(pos)!= chains:
	        if T == None:
		        T_tmp = np.random.uniform(5500,6000)
		        #T_tmp = np.log10(T_tmp)
	        else:
		        T_tmp =   np.random.normal(T,20)
		        #T_tmp = np.log10(T_tmp)
	        if M == None:
		        M_tmp = np.random.uniform(-0.1,0.1)
	        else:
		        M_tmp =   np.random.normal(M,0.01)
	        if L == None:
		        L_tmp = np.random.uniform(4.4,4.5)
	        else:
		        L_tmp =  np.random.normal(L,0.01)
	        if V == None:
		        V_tmp = np.random.uniform(5,6)
	        else:
		        V_tmp =  np.random.normal(V,0.01)

	        pos_tmp = [T_tmp, M_tmp, L_tmp, V_tmp, np.random.uniform(200,300)]
	        #print(pos_tmp)
	        lnlike_trial = lnprob(pos_tmp,ispec_dir,temperature,metalicity,loggs,real,
		        bottom_index,top_index,weight,wavemin,wavemax,n_use,
		        hdulist,Prior_values,dwt_type,dwt_number,resolution,rtc,data,cut_regions)
	        #print('Trial function lnlike: {}'.format(lnlike_trial))
	        if lnlike_trial != -np.inf:
		        pos.append(pos_tmp)
		        if not silent:
			        bar.next()
	
        if not silent:
            	bar.finish()
        #pos =np.array([[T+np.random.uniform(0,50),M+np.random.uniform(0,0.001),L+np.random.uniform(0,0.01),V+np.random.uniform(0,1),np.random.uniform(0,15)] for i in range(chains)])


        #################################
        # Run the sampler for a burn in
        #################################
        #pos, prob, state = sampler.run_mcmc(pos,draws)
        #print('done.')

        width=30
        if not silent:
	        start_time = time.time()
	        bar = Progressbar('MCMC', max=draws)
	        print('\n\n')
	        print('    Progrss bar                         Teff          [Fe/H]         Logg           Vsini         e.t.a           Step')
	        print('-----------------------------------------------------------------------------------------------------------------------------')
        for i, result in enumerate(sampler.sample(pos, iterations=draws)):
            #return result
            bar_vals['step'] = result[0]
            n = int((width+1) * float(i) / draws)
            delta_t = time.time()-start_time#  time to do float(i) / n_steps  % of caluculations
            time_incr = delta_t/(float(i+1) / draws) # seconds per increment
            time_left = time_incr*(1- float(i) / draws)
            m, s = divmod(time_left, 60)
            h, m = divmod(m, 60)
            bar_vals['timeleft'] = "{0}h:{1}m:{2}s".format(str(int(h)).zfill(2), str(int(m)).zfill(2), '{:.2f}'.format(s).zfill(2))
            #sys.stdout.write("\r[{0}{1}] {2:.2f}% - {3}h:{4}m:{5:.2f}s\tstep: {6}".format('#' * n, ' ' * (width - n), 100*float(i) / draws ,h, m, s,i))
            if not silent:
            	bar.next()

        if not silent:
            bar.finish()
        #sys.stdout.write("\n")

        ############################################
        # Append the sampler to data for future use
        ############################################
        #data['burn_in_sampler'] = sampler

        ########################
        # Now save to fits file
        ########################

        try:
            os.remove(chain_file)
        except:
            pass

        t = Table(sampler.flatchain, names=['T_eff','MH','Logg','Vsini','beta'])
        t.add_column(Column(sampler.flatlnprobability,name='loglike'))
        indices = np.mgrid[0:chains,0:draws]
        step = indices[1].flatten()

        walker = indices[0].flatten()


        t.add_column(Column(step,name='step'))
        t.add_column(Column(walker,name='walker'))
        if chain_file=='table':
	        return t
        else:
	        t.write(chain_file)

    elif method=='single_walker':
        print('Single walker method')

        if not silent:
            bar = Progressbar('MCMC', max=draws)

        lnlike_trial = -np.inf
        tests = [T,M,L,V]
        while lnlike_trial == -np.inf:
                if tests[0] == None:
	                T_tmp = np.random.uniform(5500,6000)
	                #T_tmp = np.log10(T_tmp)
                else:
	                T_tmp = T#  np.random.normal(T,20)
	                #T_tmp = np.random.uniform(5500,6000)
                if tests[1] == None:
	                M_tmp = np.random.uniform(-0.1,0.1)
                else:
	                M_tmp =  M # np.random.normal(M,0.01)
                if tests[2] == None:
	                L_tmp = np.random.uniform(4.4,4.5)
                else:
	                L_tmp = L # np.random.normal(L,0.01)

                if tests[3] == None:
	                V_tmp = np.random.uniform(5,6)
                else:
	                V_tmp = V# np.random.normal(V,0.01)

                pos_tmp = [T_tmp, M_tmp, L_tmp, V_tmp, np.random.uniform(5000,6000)]   
                lnlike_trial = lnprob(pos_tmp,ispec_dir,temperature,metalicity,loggs,real,bottom_index,top_index,weight,wavemin,wavemax,n_use,  hdulist, Prior_values, dwt_type, dwt_number, resolution, rtc,data, cut_regions)

        acceptance = 0.
        acceptance_count = 0
        acc_frac = 0

        T,M,L,V,B,Log = np.zeros(draws),np.zeros(draws),np.zeros(draws),np.zeros(draws),np.zeros(draws),np.zeros(draws)
        T[0], M[0], L[0], V[0], B[0] = pos_tmp; Log[0] = lnlike_trial
        Tstep, Mstep, Lstep, Vstep, Bstep = 200,0.1,0.1,0.5,50 
        start_time = time.time()
       
        for draw in range(1,draws):
            #print(tests)

            if tests[0] == None:
                t = np.random.normal(T[draw-1], Tstep)
            else:
                t = T[draw-1]

            if tests[1] == None:
                m = np.random.normal(M[draw-1], Mstep)
            else:
                m = M[draw-1]


            if tests[2] == None:
                l = np.random.normal(L[draw-1], Lstep)
            else:
                l = L[draw-1]


            if tests[3] == None:
                v = np.random.normal(V[draw-1], Vstep)
            else:
                v = V[draw-1]
            #print(v)
            b = np.random.normal(B[draw-1], Bstep)

            lnlike_trial = lnprob([t,m,l,v,b],ispec_dir,temperature,metalicity,loggs,real,
	        bottom_index,top_index,weight,wavemin,wavemax,n_use,
	        hdulist,Prior_values,dwt_type,dwt_number,resolution,rtc,data,cut_regions)


            if (lnlike_trial > Log[draw-1]):
                #print('accept')
                T[draw] = t
                M[draw] = m
                L[draw] = l
                V[draw] = v
                B[draw] = b
                Log[draw] = lnlike_trial
                acceptance += 1
            else:
                u = np.random.uniform(0.,1.)
                if (u < np.exp(lnlike_trial - Log[draw-1])):
                    #print('lucky')
                    T[draw] = t
                    M[draw] = m
                    L[draw] = l
                    V[draw] = v
                    B[draw] = b
                    Log[draw] = lnlike_trial
                    acceptance += 1
                else:
                    #print('reject')
                    T[draw] = T[draw-1]
                    M[draw] = M[draw-1]
                    L[draw] = L[draw-1]
                    V[draw] = V[draw-1]
                    B[draw] = B[draw-1]
                    Log[draw] = Log[draw-1]


            # Monitor acceptance
            acceptance_count +=1
            if (acceptance_count == 100):
                acc_frac = acceptance / 100
                if (acc_frac > 0.7):
                    if tests[0] == None:
                        Tstep = Tstep*np.random.uniform(1,2)
                    if tests[1] == None:
                        Mstep = Mstep*np.random.uniform(1,2)
                    if tests[2] == None:
                        Lstep = Lstep*np.random.uniform(1,2)
                    if tests[3] == None:
                        Vstep = Vstep*np.random.uniform(1,2)

                    acceptance_count,acceptance = 0 ,0  
                elif (acc_frac < 0.2):
                    if tests[0] == None:
                        Tstep = Tstep/np.random.uniform(1,2)
                    if tests[1] == None:
                        Mstep = Mstep/np.random.uniform(1,2)
                    if tests[2] == None:
                        Lstep = Lstep/np.random.uniform(1,2)
                    if tests[3] == None:
                        Vstep = Vstep/np.random.uniform(1,2)
                    acceptance_count,acceptance = 0 ,0  
            #print(B[draw])
            bar_vals['step'] = [[T[draw], M[draw], L[draw], V[draw], B[draw]], [T[draw], M[draw], L[draw], V[draw], B[draw]]]
            delta_t = time.time()-start_time#  time to do float(i) / n_steps  % of caluculations
            time_incr = delta_t/(float(draw+1) / draws) # seconds per increment
            time_left = time_incr*(1- float(draw) / draws)
            minutes, seconds = divmod(time_left, 60)
            hours, minutes = divmod(minutes, 60)
            bar_vals['timeleft'] = "{0}h:{1}m:{2}s".format(str(int(hours)).zfill(2), str(int(minutes)).zfill(2), '{:.2f}'.format(seconds).zfill(2))
            #sys.stdout.write("\r[{0}{1}] {2:.2f}% - {3}h:{4}m:{5:.2f}s\tstep: {6}".format('#' * n, ' ' * (width - n), 100*float(i) / draws ,h, m, s,i))
            if not silent:
            	bar.next()                    
                        
                 
             
            #print(T[i], M[i], L[i], V[i], B[i])
        print('Acceptance fraction = {}'.format(acceptance / draws))


        t = Table(np.array([T.tolist(), M.tolist(), L.tolist(), V.tolist(), B.tolist()]).T, names=['T_eff','MH','Logg','Vsini','beta'])
        t.add_column(Column(Log,name='loglike'))
        t.add_column(Column(np.arange(0,draws,1),name='step'))
        t.add_column(Column(np.ones(draws),name='walker'))
        if chain_file=='table':
	        return t
        else:
	        t.write(chain_file)

    elif method=='single_walker+emcee':
        print('\n\nSingle walker method then emcee')

        print('\nFirst, we run a burn in of 5000 steps with a single walker...')
        if not silent:
            print('\n\n')
            print('    Progrss bar                         Teff          [Fe/H]         Logg           Vsini         e.t.a           Step')
            print('-----------------------------------------------------------------------------------------------------------------------------')
            bar = Progressbar('MCMC', max=5000)

        lnlike_trial = -np.inf
        tests = [T,M,L,V]
        while lnlike_trial == -np.inf:
                if tests[0] == None:
	                T_tmp = np.random.uniform(5500,6000)
	                #T_tmp = np.log10(T_tmp)
                else:
	                T_tmp = T#  np.random.normal(T,20)
	                #T_tmp = np.random.uniform(5500,6000)
                if tests[1] == None:
	                M_tmp = np.random.uniform(-0.1,0.1)
                else:
	                M_tmp =  M # np.random.normal(M,0.01)
                if tests[2] == None:
	                L_tmp = np.random.uniform(4.4,4.5)
                else:
	                L_tmp = L # np.random.normal(L,0.01)

                if tests[3] == None:
	                V_tmp = np.random.uniform(5,6)
                else:
	                V_tmp = V# np.random.normal(V,0.01)

                pos_tmp = [T_tmp, M_tmp, L_tmp, V_tmp, np.random.uniform(5000,6000)]   
                lnlike_trial = lnprob(pos_tmp,ispec_dir,temperature,metalicity,loggs,real,bottom_index,top_index,weight,wavemin,wavemax,n_use,  hdulist, Prior_values, dwt_type, dwt_number, resolution, rtc,data, cut_regions)

        acceptance = 0.
        acceptance_count = 0
        acc_frac = 0

        T,M,L,V,B,Log = np.zeros(5000),np.zeros(5000),np.zeros(5000),np.zeros(5000),np.zeros(5000),np.zeros(5000)
        T[0], M[0], L[0], V[0], B[0] = pos_tmp; Log[0] = lnlike_trial
        Tstep, Mstep, Lstep, Vstep, Bstep = 200,0.1,0.1,0.5,50 
        start_time = time.time()
       
        for draw in range(1,5000):
            #print(tests)

            if tests[0] == None:
                t = np.random.normal(T[draw-1], Tstep)
            else:
                t = T[draw-1]

            if tests[1] == None:
                m = np.random.normal(M[draw-1], Mstep)
            else:
                m = M[draw-1]


            if tests[2] == None:
                l = np.random.normal(L[draw-1], Lstep)
            else:
                l = L[draw-1]


            if tests[3] == None:
                v = np.random.normal(V[draw-1], Vstep)
            else:
                v = V[draw-1]
            #print(v)
            b = np.random.normal(B[draw-1], Bstep)

            lnlike_trial = lnprob([t,m,l,v,b],ispec_dir,temperature,metalicity,loggs,real,
	        bottom_index,top_index,weight,wavemin,wavemax,n_use,
	        hdulist,Prior_values,dwt_type,dwt_number,resolution,rtc,data,cut_regions)


            if (lnlike_trial > Log[draw-1]):
                #print('accept')
                T[draw] = t
                M[draw] = m
                L[draw] = l
                V[draw] = v
                B[draw] = b
                Log[draw] = lnlike_trial
                acceptance += 1
            else:
                u = np.random.uniform(0.,1.)
                if (u < np.exp(lnlike_trial - Log[draw-1])):
                    #print('lucky')
                    T[draw] = t
                    M[draw] = m
                    L[draw] = l
                    V[draw] = v
                    B[draw] = b
                    Log[draw] = lnlike_trial
                    acceptance += 1
                else:
                    #print('reject')
                    T[draw] = T[draw-1]
                    M[draw] = M[draw-1]
                    L[draw] = L[draw-1]
                    V[draw] = V[draw-1]
                    B[draw] = B[draw-1]
                    Log[draw] = Log[draw-1]


            # Monitor acceptance
            acceptance_count +=1
            if (acceptance_count == 100):
                acc_frac = acceptance / 100
                if (acc_frac > 0.7):
                    if tests[0] == None:
                        Tstep = Tstep*np.random.uniform(1,2)
                    if tests[1] == None:
                        Mstep = Mstep*np.random.uniform(1,2)
                    if tests[2] == None:
                        Lstep = Lstep*np.random.uniform(1,2)
                    if tests[3] == None:
                        Vstep = Vstep*np.random.uniform(1,2)

                    acceptance_count,acceptance = 0 ,0  
                elif (acc_frac < 0.2):
                    if tests[0] == None:
                        Tstep = Tstep/np.random.uniform(1,2)
                    if tests[1] == None:
                        Mstep = Mstep/np.random.uniform(1,2)
                    if tests[2] == None:
                        Lstep = Lstep/np.random.uniform(1,2)
                    if tests[3] == None:
                        Vstep = Vstep/np.random.uniform(1,2)
                    acceptance_count,acceptance = 0 ,0  
            #print(B[draw])
            bar_vals['step'] = [[T[draw], M[draw], L[draw], V[draw], B[draw]], [T[draw], M[draw], L[draw], V[draw], B[draw]]]
            delta_t = time.time()-start_time#  time to do float(i) / n_steps  % of caluculations
            time_incr = delta_t/(float(draw+1) / 5000) # seconds per increment
            time_left = time_incr*(1- float(draw) / 5000)
            minutes, seconds = divmod(time_left, 60)
            hours, minutes = divmod(minutes, 60)
            bar_vals['timeleft'] = "{0}h:{1}m:{2}s".format(str(int(hours)).zfill(2), str(int(minutes)).zfill(2), '{:.2f}'.format(seconds).zfill(2))
            #sys.stdout.write("\r[{0}{1}] {2:.2f}% - {3}h:{4}m:{5:.2f}s\tstep: {6}".format('#' * n, ' ' * (width - n), 100*float(i) / draws ,h, m, s,i))
            if not silent:
            	bar.next()                    
                        
                 
             
            #print(T[i], M[i], L[i], V[i], B[i])
        print('Acceptance fraction = {}'.format(acceptance / draws))


        t = Table(np.array([T.tolist(), M.tolist(), L.tolist(), V.tolist(), B.tolist()]).T, names=['T_eff','MH','Logg','Vsini','beta'])
        #return T,M,L,V,Log
        t.add_column(Column(Log,name='loglike'))
        t.add_column(Column(np.arange(0,5000,1),name='step'))
        t.add_column(Column(np.ones(5000),name='walker'))
        try:
            os.remove(chain_file[:-5]+'_single_walker.fits')
        except:
	        t.write(chain_file[:-5]+'_single_walker.fits', format='fits')

        # Now get best position for emcee
        best_pos = np.argmax(t['loglike'])
        T_best,M_best,L_best,V_best,B_best, l_best,s_best,w_best = [i for i in t[best_pos] ]

        print('\n\nThe best parameters from burn-in are:')
        print('Teff:   {} K'.format(int(T_best)))
        print('[Fe/H]: {:.2f} dex'.format(M_best))
        print('Log g: {:.2f} dex'.format(L_best))
        print('Vsini: {:.2f} km/s\n\n'.format(V_best))

        print('Now passing best positions over to emcee...\n')


        model=0
        sampler,pos=0,0

        chains = 100
        sampler = emcee.EnsembleSampler(chains, ndim, lnprob, args = [ispec_dir,temperature,metalicity,loggs,real,bottom_index,top_index,weight,wavemin,wavemax,n_use,hdulist,Prior_values,dwt_type,dwt_number,resolution,rtc,cut_regions,data], threads=threads)


        pos = []
        print('\n')
        if not silent:
	        bar = Bar('Bulding starting positions',max = chains)
        while len(pos)!= chains:
            if tests[0] ==None:
                T_tmp =   np.random.normal(T_best,20)
            else:
                T_tmp = T_best
            
            if tests[1] == None:
                M_tmp =   np.random.normal(M_best,0.01)
            else:
                M_tmp  = M_best

            if tests[2] == None:
                L_tmp =  np.random.normal(L_best,0.01)
            else:
                L_tmp = L_best
            
            if tests[3] == None:
                V_tmp =  np.random.normal(V_best,0.01)
            else:
                V_tmp = V_best

            B_tmp = np.random.normal(B_best,0.01)

            pos_tmp = [T_tmp, M_tmp, L_tmp, V_tmp, B_tmp]
            #print(pos_tmp)
            lnlike_trial = lnprob(pos_tmp,ispec_dir,temperature,metalicity,loggs,real,
                bottom_index,top_index,weight,wavemin,wavemax,n_use,
                hdulist,Prior_values,dwt_type,dwt_number,resolution,rtc,data,cut_regions)
            #print('Trial function lnlike: {}'.format(lnlike_trial))
            if lnlike_trial != -np.inf:
                pos.append(pos_tmp)
                if not silent:
                    bar.next()
	
        if not silent:
            	bar.finish()
        #pos =np.array([[T+np.random.uniform(0,50),M+np.random.uniform(0,0.001),L+np.random.uniform(0,0.01),V+np.random.uniform(0,1),np.random.uniform(0,15)] for i in range(chains)])


        #################################
        # Run the sampler for a burn in
        #################################
        #pos, prob, state = sampler.run_mcmc(pos,draws)
        #print('done.')

        width=30
        if not silent:
	        start_time = time.time()
	        bar = Progressbar('MCMC', max=draws)
	        print('\n\n')
	        print('    Progrss bar                         Teff          [Fe/H]         Logg           Vsini         e.t.a           Step')
	        print('-----------------------------------------------------------------------------------------------------------------------------')
        for i, result in enumerate(sampler.sample(pos, iterations=draws)):
            #return result
            bar_vals['step'] = result[0]
            n = int((width+1) * float(i) / draws)
            delta_t = time.time()-start_time#  time to do float(i) / n_steps  % of caluculations
            time_incr = delta_t/(float(i+1) / draws) # seconds per increment
            time_left = time_incr*(1- float(i) / draws)
            m, s = divmod(time_left, 60)
            h, m = divmod(m, 60)
            bar_vals['timeleft'] = "{0}h:{1}m:{2}s".format(str(int(h)).zfill(2), str(int(m)).zfill(2), '{:.2f}'.format(s).zfill(2))
            #sys.stdout.write("\r[{0}{1}] {2:.2f}% - {3}h:{4}m:{5:.2f}s\tstep: {6}".format('#' * n, ' ' * (width - n), 100*float(i) / draws ,h, m, s,i))
            if not silent:
            	bar.next()

        if not silent:
            bar.finish()
        #sys.stdout.write("\n")

        ############################################
        # Append the sampler to data for future use
        ############################################
        #data['burn_in_sampler'] = sampler

        ########################
        # Now save to fits file
        ########################

        try:
            os.remove(chain_file)
        except:
            pass

        t = Table(sampler.flatchain, names=['T_eff','MH','Logg','Vsini','beta'])
        t.add_column(Column(sampler.flatlnprobability,name='loglike'))
        indices = np.mgrid[0:chains,0:draws]
        step = indices[1].flatten()

        walker = indices[0].flatten()


        t.add_column(Column(step,name='step'))
        t.add_column(Column(walker,name='walker'))
        if chain_file=='table':
	        return t
        else:
	        t.write(chain_file)



    if email_notification==True:
        try:
            import samemail
            samemail.email(SUBJECT='Wavelet analysis is complete', TEXT = 'The wavelet analysis in folder:\n{}\n is complete.'.format(os.getcwd()))
            print('Email notification sent')
        except:
            print('Email notification could not be sent')
    return data


