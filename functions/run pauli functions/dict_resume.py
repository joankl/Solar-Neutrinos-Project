'''
Python script to read files in pickel format and save the resume of many
chunk data from partition jobs (exm: HS and atm analysis).
This script is designed to receive the main directory and names of the
dictionaries to be resumed.
'''

import glob
import os
import pickle
import numpy as np

main_dir = '/lstore/sno/joankl/solar_analysis/real_data/2p2ppo/ntuple_output/veto_HS_atm_candidates/output_*/'
save_dir = '/lstore/sno/joankl/solar_analysis/real_data/2p2ppo/ntuple_output/resume_files/'

# Formato de salida: 'npz' o 'pkl'
save_format = 'npz'

#dict_list = ['prompt_coinc_dict', 'delay_coinc_dict'] # Name of the dictionaries to be readen
#dict_keys = ['runID', 'eventID', 'energy', 'posx', 'posy', 'posz', 'time']

dict_list = ['atm_dict', 'hs_dict'] # Name of the dictionaries to be readen
dict_keys = ['eventID_prompt', 'eventID_followers', 'runID_prompt', 'runID_followers']


def _to_1d_arrays(chunks):
	"""Convierte una lista de chunks (listas, tuplas o ndarrays) a una lista de arrays 1D."""
	arrays = []
	for chunk in chunks:
		arr = np.asarray(chunk)
		if arr.size > 0:
			arrays.append(arr.ravel())
	return arrays


if __name__ == "__main__":
	if save_format not in ('pkl', 'npz'):
		raise ValueError(f"Format not supported: '{save_format}'. Use 'npz' o 'pkl'.")

	os.makedirs(save_dir, exist_ok=True)

	#Loop over each dictionary to save its observables
	for dict_i in dict_list:

		print(f'reading dictionary {dict_i}.pkl')

		#Empty dictionary to be saved. It will contain the resumed dict information for each dict. in dict list.
		dict_to_save = {var_i: [] for var_i in dict_keys}

		# List of directories of the dict_i
		full_dict_dir_list = sorted(glob.glob(os.path.join(main_dir, dict_i + '.pkl')))
		#print(f'list of files to be readen: {full_dict_dir_list}')

		# Loop over the directory list
		for dir_i in full_dict_dir_list:

			# Load Dictionary file
			with open(dir_i, 'rb') as f:
				dict_file = pickle.load(f)

			# Loop over the dictionary keys and append to the dict_to_save
			try:
				for var_i in dict_keys:
					dict_to_save[var_i].append(dict_file[var_i])

			except KeyError:
				print(f'No key {var_i} on the loaded dictionary: {dir_i}. Skipping dictionary.')
				continue

		if save_format == 'pkl':
			# Flat the dictionary lists
			for var_i in dict_keys:
				val = dict_to_save[var_i]

				# Only flat list of lists. Otherwise is unnecesary.
				if (isinstance(val, list) and len(val) > 0 and isinstance(val[0], (list, tuple))):
					dict_to_save[var_i] = [x for sublist in dict_to_save[var_i] for x in sublist]

				print(f'final variable {var_i} with dict entry {dict_to_save[var_i]}')

			# Save the resumed dict
			out_path = os.path.join(save_dir, dict_i + '.pkl')
			with open(out_path, 'wb') as f:
				pickle.dump(dict_to_save, f)

			print(f'dictionary {dict_i}.pkl saved succesfully!')

		elif save_format == 'npz':
			arrays_to_save = {}
			for var_i in dict_keys:
				chunk_arrays = _to_1d_arrays(dict_to_save[var_i])
				if len(chunk_arrays) > 0:
					merged_arr = np.concatenate(chunk_arrays)
				else:
					merged_arr = np.array([])
				arrays_to_save[var_i] = merged_arr
				print(f'final variable {var_i} | shape={merged_arr.shape} | dtype={merged_arr.dtype}')

			out_path = os.path.join(save_dir, dict_i + '.npz')
			np.savez(out_path, **arrays_to_save)

			print(f'dictionary {dict_i}.npz saved succesfully!')



