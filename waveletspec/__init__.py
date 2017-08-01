import waveletspec,os,sys

if (sys.version_info > (3, 0)):
    from urllib.request import urlretrieve
elif (sys.version_info < (3, 0)):
    from urllib import urlretrieve

def reporthook(blocknum, blocksize, totalsize):
    readsofar = blocknum * blocksize
    if totalsize > 0:
        percent = readsofar * 1e2 / totalsize
        s = "\r%5.1f%% %*d / %d" % (
            percent, len(str(totalsize)), readsofar, totalsize)
        sys.stderr.write(s)
        if readsofar >= totalsize: # near the end
            sys.stderr.write("\n")
    else: # total size is unknown
        sys.stderr.write("read %d\n" % (readsofar/1e6,))




current_dir = os.getcwd()
this_dir, this_filename = os.path.split(waveletspec.__file__)


##############################
# First search for the iSpec
##############################
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
    # Downlaod the right one, untar it and link it
    print('I cant find the iSpec directory, searching for the tar files...')
    if (sys.version_info > (3, 0)):
        try:
            os.stat(this_dir+'/iSpec_v20160930_py3.tar.gz')
            print('Found iSpec_v20160930_py3.tar.gz, unpacking ...')
            os.system('tar -xf {}/iSpec_v20160930_py3.tar.gz -C {}'.format(this_dir,this_dir))
            print('Exctracted to {}/iSpec_v20160930_py3'.format(this_dir))
        except:
            print('No cached version of iSpec (py3) found, downloading from sams server . . .')
            # First download iSpec in waveletspec directory
            os.chdir(this_dir)
            urlretrieve("http://www.astro.keele.ac.uk/~sgill/iSpec_v20160930_py3.tar.gz", 'iSpec_v20160930_py3.tar.gz', reporthook)
            #download_file("http://www.astro.keele.ac.uk/~sgill/iSpec_v20160930_py3.tar.gz","iSpec_v20160930_py3.tar")
            os.chdir(current_dir)


            os.stat(this_dir+'/iSpec_v20160930_py3.tar.gz')
            print('Found iSpec_v20160930_py3.tar.gz, unpacking ...')
            os.system('tar -xf {}/iSpec_v20160930_py3.tar.gz -C {}'.format(this_dir,this_dir))
            print('Exctracted to {}/iSpec_v20160930_py3'.format(this_dir))
        sys.path.insert(0, os.path.abspath((this_dir+'/iSpec_v20160930_py3')))
        
    elif (sys.version_info < (3, 0)):
        try:
            os.stat(this_dir+'/iSpec_v20160930_py2.tar.gz')
            print('Found iSpec_v20160930_py2.tar.gz, unpacking ...')
            os.system('tar -xf {}/iSpec_v20160930_py2.tar.gz -C {}'.format(this_dir,this_dir))
            print('Exctracted to {}/iSpec_v20160930_py2'.format(this_dir))
        except:
            print('No cached version of iSpec (py2) found, downloading from sams server . . .')
            # First download iSpec in waveletspec directory
            os.chdir(this_dir)
            urlretrieve("http://www.astro.keele.ac.uk/~sgill/iSpec_v20160930_py2.tar.gz", 'iSpec_v20160930_py2.tar.gz', reporthook)
            os.chdir(current_dir)


            os.stat(this_dir+'/iSpec_v20160930_py2.tar.gz')
            print('Found iSpec_v20160930_py2.tar.gz, unpacking ...')
            os.system('tar -xf {}/iSpec_v20160930_py2.tar.gz -C {}'.format(this_dir,this_dir))
            print('Exctracted to {}/iSpec_v20160930_py2'.format(this_dir))
        sys.path.insert(0, os.path.abspath((this_dir+'/iSpec_v20160930_py2')))
        
    else:
        print('I CANT FIND ISPEC')




#################################
# Now need to check for the grid
#################################
try:
    os.stat(this_dir+'/Grids/spectrum.fits')
    print('Found the default grid, "spectrum.fits"')
except:
    print('No default grid found,  downloading from sams server . . .')
    # First download iSpec in waveletspec directory
    os.chdir(this_dir+'/Grids')
    urlretrieve("http://www.astro.keele.ac.uk/~sgill/waveletspec_grids/spectrum.fits","spectrum.fits", reporthook)
    print('spectrum.fits downloaded to {}/Grids'.format(this_dir))
    os.chdir(current_dir)




import ispec
from waveletspec.check_fit import check_fit
from waveletspec.coadd_spectra import coadd_spectra
from waveletspec.interpolate_grid import interpolate_grid
from waveletspec.load_spectra import load_spectra
from waveletspec.normalise_spectra import normalise_spectra
from waveletspec.plot_spectra import plot_spectra
from waveletspec.resample import resample
from waveletspec.wavelet_analysis import wavelet_analysis,estimate_turbulence_from_doyle
from waveletspec.write_spectra import write_spectra
from waveletspec.check_fit import check_fit
from waveletspec.wavelet_analysis_plot import wavelet_analysis_plot
from waveletspec.clean_spectra import clean_spectrum, clean_telluric_regions, filter_cosmic_rays


