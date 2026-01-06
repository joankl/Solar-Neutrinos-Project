'''
Python script to read and save the resume of the
hotspot, atmospheric and IBD coincidence analysis.
This script is designed to receive the main directory and names of the
dictionaries to be resumed. 
'''

import glob
import pickle


main_dir = '/lstore/sno/joankl/solar_analysis/real_data/bisMSB/Analysis15/results/bkg_candidates/output_*/'
save_dir = '/lstore/sno/joankl/solar_analysis/real_data/bisMSB/Analysis15/results/resume_files/'

#dict_list = ['prompt_IBD_dict', 'delay_IBD_dict'] # Name of the dictionaries to be readen
#dict_keys = ['runID', 'eventID', 'energy', 'posx', 'posy', 'posz', 'time']

dict_list = ['atm_dict', 'hs_dict'] # Name of the dictionaries to be readen
dict_keys = ['counter', 'eventID', 'runID']


#Loop over each dictionary to save its observables
for dict_i in dict_list:

	print(f'reading dictionary {dict_i}.pkl')

	#Empty dictionary to be saved. It will contain the resumed dict information for each dict. in dict list.
	dict_to_save = {var_i: [] for var_i in dict_keys}

	# List of directories of the dict_i
	full_dict_dir_list = glob.glob(main_dir + dict_i + '.pkl')
	#print(f'list of files to be readen: {full_dict_dir_list}')  

	# Loop over the directory list
	for dir_i in full_dict_dir_list:

		# Load Dictionary file
		with open(dir_i, 'rb') as f:
			dict_file = pickle.load(f)

		# Loop over the dictionary keys and append to the dict_to_save
		for var_i in dict_keys:
			dict_to_save[var_i].append(dict_file[var_i])

	# Flat the dictionary lists
	for var_i in dict_keys:
		val = dict_to_save[var_i]

		# Only flat list of lists. Otherwise is unnecesary.
		if (isinstance(val, list) and len(val) > 0 and isinstance(val[0], (list, tuple))):
			dict_to_save[var_i] = [x for sublist in dict_to_save[var_i] for x in sublist]

		print(f'final variable {var_i} with dict entry {dict_to_save[var_i]}')

	# Save the resumed dict
	with open(save_dir + dict_i + '.pkl', 'wb') as f:
		pickle.dump(dict_to_save, f)

	print(f'dictionary {dict_i}.pkl saved succesfully!')



