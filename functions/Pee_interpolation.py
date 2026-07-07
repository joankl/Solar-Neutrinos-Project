'''
Script specially designed to return the solar electron
neutrino survival probability given an array of neutrino
energies.

Edit Dates:
- 07/07/2026: Creation of script
'''

import numpy as np
from scipy.interpolate import interp1d
import glob
import re
import os

def orden_natural(archivo):
	"""Función para ordenar archivos naturalmente por número de run/subrun."""
	return [int(texto) if texto.isdigit() else texto.lower() for texto in re.split('(\d+)', archivo)]


def Pee_interp(energy_nu_dir, pselmaa_data_dir, save_dir, f_idx = 0):

	'''
	Function designe to compute and the survival probability of 
	solar electron neutrinos given the neutrino energies.

	Parameters:
	- energy_nu_dir: directory where the neutrino energy np.array is.
	- pselmaa_data_dir: directory where the PSELMAA data is
	- save_dir: directory to save the survival probability array.
	- f_idx: file counter
	'''

	os.makedirs(save_dir, exist_ok=True)

	# ======= Load the data =======
	PSelmaa_data = np.loadtxt(pselmaa_data_dir, skiprows=1)
	Pee = PSelmaa_data[:,1]
	Pee_energy = PSelmaa_data[:,0]

	E_nu = np.load(energy_nu_dir) # Simulated neutrino energy

	# ======= Interpolate Pee with the MC data energy =======
	Pee_f = interp1d(Pee_energy, Pee, kind='linear', bounds_error=False, fill_value='extrapolate')  #function that describes the interpolation
	Pee_int = Pee_f(E_nu) # Interpolation Values

	print(f'Interpolated Pee with values: \n {Pee_int}')

	print(f'Saving results on {save_dir}')

	np.save(os.path.join(save_dir, f'Pee_{f_idx}'), Pee_int)

if __name__ == "__main__":

	energy_nu_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/2p2_ppo/solar_8BNue/ratds_output/np_files/'
	pselmaa_data_dir = '/lstore/sno/joankl/solar_analysis/flux_prediction/PSelmaa/output_files/pselmaa_test_sun_pee_B16_GS98_b8.txt'
	save_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/2p2_ppo/solar_8BNue/ratds_output/np_files/'

	energy_nu_list = glob.glob(energy_nu_dir + 'parent_energy_*.npy')
	energy_nu_list.sort(key=orden_natural)

	print(f'list of files: {energy_nu_list}')

	for idx, file_i in enumerate(energy_nu_list):
		print(f'reading file {file_i}')

		# Busca un guion bajo '_', seguido de dígitos '(\d+)', seguido de '.npy' al final '$'
		match = re.search(r'_(\d+)\.npy$', file_i)
		if match:
			file_index = int(match.group(1)) # El grupo 1 es el número capturado
			print(f'file_index: {file_index}')

		Pee_interp(file_i, pselmaa_data_dir, save_dir, f_idx = file_index)

	print('Done :)')