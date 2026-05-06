'''
Python module useful for reading the results saved in the folders output_*
and return the concatenated arrays which will be saved in a specific resume
directory
'''

import numpy as np
import glob
import os


def generate_resume(fdir_pattern, var_name, save_dir):

	'''
	Function that takes the directory patter of the var_name to be resumed in
	a concatenated form.

	Parameters:
	-fdir_pattern: Main directory of the folders that contain the files to be resumed.
	               This parameter is asumed to has a pattern form of output_*.
	-var_name: Name of the variable to be resumed.
	-save_dir: Directory where the resume files will be saved.
	'''

	# Construct the list of var_name np.array() directories to be readen

	os.mkdir(save_dir)

	f_list = sorted(glob.glob(fdir_pattern + var_name + '.npy'))

	print(f'list of files to be readen {f_list}')

	resumed_array = [] # Empty list to be filled with all the array values within the file dir. loop and saved at the end

	# Loop over the file list
	for f_i in f_list:
		array_var = np.load(f_i)
		resumed_array.append(array_var)

	resumed_array = np.concatenate(resumed_array)	

	# Save the array
	print('Saving resumed array')
	np.save(save_dir + str(var_name) + '.npy', resumed_array)

	return print('resume files obtined!')


if __name__ == "__main__":

    fdir_pattern = '/lstore/sno/joankl/solar_analysis/real_data/2p2ppo/ntuple_output/output_files/output_*/'
    var_name_list = ['energy', 'posr_av', 'posx', 'posy', 'posz_av', 'itr', 'runID', 'subrunID', 'eventID']
    save_dir = '/lstore/sno/joankl/solar_analysis/real_data/2p2ppo/ntuple_output/resume_files/'

    for var_i in var_name_list:
    	generate_resume(fdir_pattern, var_i, save_dir)


