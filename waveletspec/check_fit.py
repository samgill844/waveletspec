import matplotlib.pyplot as plt
def check_fit(list_of_data, list_of_labels):
	
	for i in range(len(list_of_data)):
		plt.plot(list_of_data[i]['waveobs'],list_of_data[i]['flux'], label = list_of_labels[i])


	plt.xlabel('Wavelength')
	plt.ylabel('Flux')
	plt.grid()

	plt.legend()
	plt.show()
