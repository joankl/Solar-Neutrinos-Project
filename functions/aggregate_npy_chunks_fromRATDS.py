import os
import re
import glob
import numpy as np

def aggregate_npy_chunks(input_dir, output_npz, selected_observables=None):
	"""
	Aggregate np files saven in chunks based on their name and chunk index. It will
	concatenate the arrays in an unique array and save the results on a npz file

	Parameters*:

	- input_dir (str): directory where the chunk numpy data is
	- output_npz (str): directory to save the npz files
	- selected_observables (list or None): List with the name of observables of interest to be summarized.
										   if None, it will summarize all the files available in input_dir
	"""
	# Regular patter of the chunk files to separate observable name from its index
	# example: 'cos_alpha_105.npy' -> obs='cos_alpha', idx='105'
	pattern = re.compile(r'^(?P<obs>.+?)_(?P<idx>\d+)\.npy$')
	
	obs_files = {}
	
	# look for the .npy file on the directory
	for file_path in glob.glob(os.path.join(input_dir, '*.npy')):
		filename = os.path.basename(file_path)
		match = pattern.match(filename)
		if match:
			obs_name = match.group('obs')
			index = int(match.group('idx'))
			
			# Filtrer the observables indicated in the function argument 'selected_observables'
			if selected_observables is not None and obs_name not in selected_observables:
				continue
				
			if obs_name not in obs_files:
				obs_files[obs_name] = []
			obs_files[obs_name].append((index, file_path))
			
	data_to_save = {}
	
	# Process each chunk group found
	for obs_name, files_with_indices in obs_files.items():
		# sort by numerical index order
		files_with_indices.sort(key=lambda x: x[0])
		
		print(f"Processing observable: '{obs_name}' ({len(files_with_indices)} chunks)...")
		
		chunk_arrays = []
		for idx, fpath in files_with_indices:
			arr = np.load(fpath)
			chunk_arrays.append(arr)
		
		# Concatenate all the chunks of the observable as a flat array
		if chunk_arrays:
			data_to_save[obs_name] = np.concatenate(chunk_arrays, axis=0)
			print(f"  -> Final shape of the array '{obs_name}': {data_to_save[obs_name].shape}")
			
	# Save the results on a unique npz file
	if data_to_save:
		np.savez_compressed(output_npz, **data_to_save)
		print(f"\nFile .npz successfully generated at: {output_npz}!")
	else:
		print("No valid filed found to perform the aggregation task!.")

if __name__ == "__main__":

	input_dir = "/lstore/sno/joankl/solar_analysis/mc_data/2p2_ppo/solar_8BNue/ratds_output/np_files/chunk_files/"                        
	output_dir = "/lstore/sno/joankl/solar_analysis/mc_data/2p2_ppo/solar_8BNue/ratds_output/np_files/resume_files/"
	output_fname = 'resume_ratDS_analysis_solar8BNue.npz'

	selected_observables= ['energy_corr', 'hit_qhs', 'posr', 'cos_alpha', 'hit_residual'] 
	
	aggregate_npy_chunks(input_dir, output_dir + output_fname, selected_observables)