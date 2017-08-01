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
import numpy as np


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


def load_spectra(name_of_file, rv_correct=False,velocity_step=1.0):
	data = 0
	try:
		data = ispec.read_spectrum(name_of_file) # using ispec
	except:
		print('File: '+name_of_file+' cannot be opend by iSpec')
		print('This may be corrept or not in the right form.')
		return
	if min(data['waveobs'])>1000:
		data['waveobs'] = data['waveobs']/10 # convert to nm

	if rv_correct == True:
		rv, rv_err = determine_radial_velocity_with_template(data,velocity_step)
		print('{}: {} +- {}'.format(name_of_file,rv,rv_err))
		data = ispec.correct_velocity(data, rv)

	data['err'] = [0.01]*len(data)

	return data
