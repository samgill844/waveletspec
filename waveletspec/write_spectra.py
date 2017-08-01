import numpy as np

def write_spectra(data,name='spectrum'):
	np.savetxt(name, np.array([data['waveobs'],data['flux'],data['err']]).transpose(),fmt='%.4f')
