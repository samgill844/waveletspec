
################################################################################################################
# SECTION 0 - IMPORTS 
################################################################################################################

#from waveletspec.input import *

import waveletspec,os,sys
this_dir, this_filename = os.path.split(waveletspec.__file__)


import waveletspec,os,sys
this_dir, this_filename = os.path.split(waveletspec.__file__)
print('ISPEC_1: {}/iSpec_v20160930_py2'.format(this_dir))
# First search for the iSpec tar (eithe Python2 or 3 version)
try:
    if (sys.version_info < (3, 0)):
        os.stat(this_dir+'/iSpec_v20160930_py2')
        print('Found python 2 version of iSpec')
        sys.path.insert(0, os.path.abspath((this_dir+'/iSpec_v20160930_py2')))
    if (sys.version_info > (3, 0)):
        os.stat(this_dir+'/iSpec_v20160930_py3')
        print('Found python 3 version of iSpec')
        sys.path.insert(0, os.path.abspath((this_dir+'/iSpec_v20160930_py3')))
except OSError as e:
    # export the tar
    print('I cant find the iSpec directory, searching for the tar files...')
    if (sys.version_info > (3, 0)):
        os.system('tar -xf {}/iSpec_v20160930_py3.tar.gz -C {}'.format(this_dir,this_dir))
        print('Exctracted to {}/iSpec_v20160930_py3'.format(this_dir))
        sys.path.insert(0, os.path.abspath((this_dir+'/iSpec_v20160930_py3')))
        
    elif (sys.version_info < (3, 0)):
        os.system('tar -xf {}/iSpec_v20160930_py2.tar.gz -C {}'.format(this_dir,this_dir))
        print('Exctracted to {}/iSpec_v20160930_py2'.format(this_dir))
        sys.path.insert(0, os.path.abspath((this_dir+'/iSpec_v20160930_py2')))
    else:
        print('I CANT FIND ISPEC')
import ispec
from numba import autojit, jit, vectorize
import time, os, sys
from heapq import nsmallest
from scipy.signal import fftconvolve, convolve
import numpy as np
from scipy.interpolate import interp1d
from astropy.io import fits



################################################################################################################
# SECTION 1 - Jitted functions 
################################################################################################################

@autojit
def vsini_jitted(vsini, velocity_step):
    epsilon=0.6
    #wave_ = np.log(wavelength)
    #velo_ = np.linspace(wave_[0],wave_[-1],len(wave_))
    #flux_ = np.interp(lnwave_,wave_,flux)
    #dvelo = velo_[1]-velo_[0]
    dvelo = velocity_step/(299792458.*1e-3)
    vrot = vsini/(299792458.*1e-3) # in Km/s
    #-- compute the convolution kernel and normalise it
    # PFLM/SG - changed calculation of kernel to use an odd number of pixels so that 
    # there is at least one pixel in the kernel - more reliable for low vrot
    n = np.int8(1+2*vrot/dvelo) # always have at least one pixel in the kernel !
    #n = int(3+2*vrot/dvelo) # tried at least 3 elements - it fails
    velo_k = np.arange(0,n)*dvelo
    velo_k -= velo_k[-1]/2.
    y = 1 - (velo_k/vrot)**2 # transformation of velocity
    vsini_kernel = (2*(1-epsilon)*np.sqrt(y)+np.pi*epsilon/2.*y)/(np.pi*vrot*(1-epsilon/3.0))  # the kernel
    vsini_kernel /= vsini_kernel.sum()
    return vsini_kernel

@autojit
def resolution_jit(resolution,velocity_step ):
    v_fwhm =(299792458.*1e-3)/resolution # n km/s
    n = int(1+4*v_fwhm/velocity_step) # always have at least one pixel in the kernel !
    velo_k = np.linspace(-n*velocity_step,n*velocity_step,2*n + 1)
    v_sigma = v_fwhm/(2*np.sqrt(2*np.log(2)))
    resolution_kernel =  np.exp(-0.5*velo_k**2/(v_sigma**2)) # the kernel
    resolution_kernel /= resolution_kernel.sum()
    return resolution_kernel



################################################################################################################
# SECTION 2 - Main broadening function 
################################################################################################################

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
    vmac_kernel,vsini_kernel,resolution_kernel = [],[],[]
    #print('\t- Logspace: {}'.format(time.time()-	tic_b_))
    tic_b = time.time()


    kernals_to_use=[]

    if (vmac <=0) and (vsini <= 0) and (resolution <= 0):
        #print('No kernals used')
        return wavelength,flux

    #########################
    # Construct vmac kernel
    #########################
    tic = time.time()
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

    #print('Vmac: {}'.format(time.time()-tic))

    #########################
    # Construct vsini kernel
    #########################
    tic = time.time()
    if vsini > 0:
        tic = time.time()
        kernals_to_use.append(vsini_jitted(vsini, velocity_step))
        #print('Vsini jit: {}'.format(time.time()-tic))

    else:
        pass

    ##############################
    # Construct resolution kernel
    ##############################
    if resolution > 0:
        tic = time.time()
        kernals_to_use.append(resolution_jit(resolution,velocity_step ))
        #print('Reolution_jit: {}'.format(time.time()-tic))

    else:
        pass


    ##########################
    # Create the master kernel
    ##########################
    #print('Number of kernels: {}'.format(len(kernals_to_use)))
    #print(kernals_to_use)

    kernals_to_use = sorted(kernals_to_use, key=len)
    kernals_to_use = kernals_to_use[::-1] # list of kernels - longest first

    master_kernel = kernals_to_use[0]
    tic_b = time.time()
    for i in range(1,len(kernals_to_use)):
        master_kernel = convolve(master_kernel,kernals_to_use[i],mode='same')


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
	




    
################################################################################################################
# SECTION 3 - Main function 
################################################################################################################

def interpolate_grid(teff,mh,logg,vrot,resolution=55000,n_values=17,grid_name='spectrum',vmac=2):
    rtc = grid_name+'.fits'
    path_to_grid = os.path.join(this_dir, "Grids", rtc)

    try:
	    hdulist = fits.open(path_to_grid)
	    #print('Using grid from:')
	    #print(path_to_grid)
    except:
	    print('I cant find the grid.')
	    print(path_to_grid)
	    return



    ########################
    # Preliminaries for interpolation
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

    ##################
    #interpolate model
    ##################
    tic_ = time.time()

    teff_2_closest = nsmallest(2, list(set(temperature)), key=lambda x: abs(x-teff)); teff_2_closest.sort() #G
    mh_2_closest = nsmallest(2, list(set(metalicity)), key=lambda x: abs(x-mh));mh_2_closest.sort() #G
    logg_2_closest = nsmallest(2, list(set(loggs)), key=lambda x: abs(x-logg));logg_2_closest.sort() #G
    # finding spectra
    models = [[i,j,k] for i in teff_2_closest for j in mh_2_closest for k in logg_2_closest]

    ticc = time.time()    
    models_index = [get_index(i,temperature,metalicity,loggs) for i in models]
    #print('get_index old: {}'.format(time.time()-ticc))

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
    #print(tstep,lstep,mstep)
    model =  interp_spectra/(tstep*lstep*mstep)
    if rtc=='spectrum.fits':
        model_wave = np.linspace(350,800,2**18)
    else:
        raise ValueError('No grid found')
    ######################
    # broaden model
    ######################
    model_wave, model = broaden_all_fast(model_wave,model,vmac,vrot,resolution)
    ##################
    # Re in-terpolate
    ##################
    f = interp1d (model_wave,model)

    if rtc=='spectrum.fits':wave = np.linspace(350,800,2**18)
    else:                          raise ValueError('No grid found')
    #print('Re-interplating: {}'.format(time.time()-tic))
    model = np.interp(wave,model_wave,model)

    tic = time.time()

    data = ispec.create_spectrum_structure(wave)
    data['flux'] = model
    data['err'] = np.array([0.05]*len(wave))
    return data
       


