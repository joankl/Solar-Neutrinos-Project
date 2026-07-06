'''
Python script designed to read and concatenate a series of
pickle files. This script will be used to concatenate the 
information from HS analysis, and save it in a resumed dictionary.

Created on: 21/05/2026
'''

import numpy as np 
import pickle
import glob
import os

# ====== Directories Definition ======

fdir = '/lstore/sno/joankl/solar_analysis/real_data/2p2ppo/ntuple_output/HS_results/pkl_files/'

HS_fname = 'hs_prompt_dict_*.pkl'
followers_fname = 'hs_delay_dict_*.pkl'

save_dir = '/lstore/sno/joankl/solar_analysis/real_data/2p2ppo/ntuple_output/HS_results/pkl_resume/'

os.makedirs(save_dir, exist_ok=True)

HS_flist = glob.glob(fdir + HS_fname)
followers_flist = glob.glob(fdir + followers_fname)

# ======= Set data to be Readen and Saved =======

dict_keys = ['eventID', 'runID']

hs_prompt_resume = {var_i: [] for var_i in dict_keys}
hs_delay_resume = {var_i: [] for var_i in dict_keys}

# ======= Load and Save the Data =======

for hs_prompt_fi, hs_delay_fi in zip(HS_flist, followers_flist):

	with open(hs_prompt_fi, "rb") as f:
		hs_prompt_dict = pickle.load(f)

		#keys_list = list(hs_prompt_dict.keys())
		for key_i in dict_keys:
			hs_prompt_resume[key_i].append(hs_prompt_dict[key_i])

	with open(hs_delay_fi, "rb") as f:
		hs_delay_dict = pickle.load(f)

		#keys_list = list(hs_delay_dict.keys())
		for key_i in dict_keys:
			hs_delay_resume[key_i].append(hs_delay_dict[key_i])

# Transform dict entries to np.array
hs_prompt_resume = {k: np.concatenate(v) for k, v in hs_prompt_resume.items() if len(v) > 0}
hs_delay_resume = {k: np.concatenate(v) for k, v in hs_delay_resume.items() if len(v) > 0}

N_HS = len(hs_prompt_resume['eventID'])
print(f'total number of HS: {N_HS}')

print('saving the dictionary in pickle format')
with open(save_dir + f'hs_prompt_resume.pkl', 'wb') as f:
	pickle.dump(hs_prompt_resume, f)

with open(save_dir + f'hs_delay_resume.pkl', 'wb') as f:
	pickle.dump(hs_delay_resume, f)

print('Saved Sucessfully')





