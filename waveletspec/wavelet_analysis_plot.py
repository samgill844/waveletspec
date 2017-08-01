import corner
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, Column


def wavelet_analysis_plot(filename='chain.fits',burn_in=None,print_vals=True,plot_vals=True):


	# Load in data
	t = Table.read(filename)
	max_step = np.max(np.array(t['step']))
	print('Read {} steps from {}'.format(max_step,filename))

	# Mask burn_in
	if burn_in != None:
		if burn_in > max_step-1:
			print('Burn-in value cannot exceed total draws ({})'.format(max_step))

		mask = np.zeros(len(t),dtype='bool')
		for i in range(len(t)):
			if t['step'][i] > burn_in:
				mask[i] = True

		t = t[mask]

	# Find the best step
	best_indexd = np.argmax(np.array(t['loglike']))
	
	# extract best
	T,M,L,V,B,l,s,w = t[best_indexd]
	dT,dM,dL,dV,dB,dl,ds,dw = np.array(t.groups.aggregate(np.std))[0]

	# Print if wanted
	if  print_vals==True:
		print('\nAtmospheric parameters:')
		print('\tT_eff: {} +- {} K'.format(int(T),int(dT)))
		print('\t[Fe/H]: {:.2f} +- {:.2f} (dex)'.format(M,dM))
		print('\tLog g: {:.2f} +- {:.2f} (dex)'.format(L,dL))
		print('\tVsini: {:.2f} +- {:.2f} km/s'.format(V,dV))


	# Corner plot if needed
	if plot_vals == True:
		t.remove_columns(['loglike', 'step', 'walker'])
		tt = t.as_array().view(np.float64).reshape(t.as_array().shape + (-1,)) # convert from structuered array to normal array
		fig = corner.corner(tt, labels =  ['T$_{eff}$ (K)','[Fe/H]','Log g','Vsin$i$ (km s$^{-1}$)',r'$\beta$'])
		return fig,T,M,L,V,B,l,s,w
	else:
		return T,M,L,V,B,l,s,w


def torres(T_eff,Logg,MH):
    # Now use torres relation to get R1 for first guess
    a1 = 1.5689
    a2 = 1.3787
    a3 = 0.4243
    a4 = 1.139
    a5 = -0.1425
    a6 = 0.01969
    a7 = 0.1010
    b1 = 2.4427
    b2 = 0.6679
    b3 = 0.1771
    b4 = 0.705
    b5 = -0.21415
    b6 = 0.02306
    b7 = 0.04173
    X = np.log10(T_eff)-4.1
    logM = a1 + a2*X + a3*(X**2) + a4*(X**3) + a5*(Logg)**2 + a6*(Logg)**3 + a7*MH
    logR = b1 + b2*X + b3*(X**2) + b4*(X**3) + b5*(Logg)**2 + b6*(Logg)**3 + b7*MH
    # Torrest estimates a 6.4% eror with mass and 3.2% error with radii
    return 10**logM, 10**logR
