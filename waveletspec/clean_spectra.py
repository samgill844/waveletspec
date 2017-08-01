import numpy as np
import waveletspec,os,sys
import os,sys
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

def clean_spectrum(data):
    flux_base = 0.0
    flux_top = 1.0
    err_base = 0.0
    err_top = 1.0
    ffilter = (data['flux'] > flux_base) & (data['flux'] <= flux_top)
    efilter = (data['err'] > err_base) & (data['err'] <= err_top)
    wfilter = np.logical_and(ffilter, efilter)
    return data[wfilter]


def clean_telluric_regions(data):
    telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"
    telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.0)

    models, ccf = ispec.cross_correlate_with_mask(data, telluric_linelist, \
                            lower_velocity_limit=-100, upper_velocity_limit=100, \
                            velocity_step=0.5, mask_depth=0.01, \
                            fourier = False,
                            only_one_peak = True)

    bv = np.round(models[0].mu(), 2) # km/s
    bv_err = np.round(models[0].emu(), 2) # km/s

    telluric_linelist_file = ispec_dir + "/input/linelists/CCF/Synth.Tellurics.500_1100nm/mask.lst"
    telluric_linelist = ispec.read_telluric_linelist(telluric_linelist_file, minimum_depth=0.0)

    # - Filter regions that may be affected by telluric lines
    #bv = 0.0
    min_vel = -30.0
    max_vel = +30.0
    # Only the 25% of the deepest ones:
    dfilter = telluric_linelist['depth'] > np.percentile(telluric_linelist['depth'], 75)
    tfilter = ispec.create_filter_for_regions_affected_by_tellurics(data['waveobs'], \
                                telluric_linelist[dfilter], min_velocity=-bv+min_vel, \
                                max_velocity=-bv+max_vel)
    return data[~tfilter]



def filter_cosmic_rays(data):
    cosmics = ispec.create_filter_cosmic_rays(data, star_continuum_model, \
                                            resampling_wave_step=0.001, window_size=15, \
                                            variation_limit=0.01)
    return data[~cosmics]


