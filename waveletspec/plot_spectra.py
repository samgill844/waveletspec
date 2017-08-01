import matplotlib.pyplot as plt

def plot_spectra(data):
	plt.ion()
	plt.plot(data['waveobs'],data['flux'])
	plt.xlabel('Wavelength (nm)')
	plt.ylabel('Flux')
	plt.show()

def plot_spectra_with_errorbars(data):
	plt.ion()
	plt.errorbar(data['waveobs'],data['flux'],yerr=data['fluxe'],fmt='k')
	plt.xlabel('Wavelength (nm')
	plt.ylabel('Flux')
	plt.show()
