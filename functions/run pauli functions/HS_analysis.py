'''
Python script designed to receive ntuple files and perform
the Hotspot (HS) event analysis. Typically used during the 2.2PPO 
phase. This script should save the GTID, runID and subrunID of 
the event compatible with HS and followers. 

Creation date: 14/05/2026

Edit:	- (20/05/2026) Add the section to submit jobs to the Grid
'''

import uproot
import numpy as np
import glob
import re
import os
import pickle
import time


# =================== Auxiliar Functions ===================

def read_files_txt(file_txt_dir):
    '''Read and file.txt and returns a list with the names of the files to be readen'''
    with open(file_txt_dir, "r", encoding="utf-8") as f:
        file_names = [line.strip() for line in f if line.strip()]
    return file_names

def natural_order(file):
    """Function to arange the file names by order of run and subrun"""
    return [int(text) if text.isdigit() else text.lower() for text in re.split('(\d+)', file)]

# ============================================================

def HS_analysis(read_dir, file_txt_dir, save_dir, fcounter = 0):

	'''
	Parameters:
	- read_dir: Directory where the data.ntuple.root is.
	- file_txt_dir: Directory where the files.txt with the names of 
					files to be read are.
	- save_dir: Directiry where the analysis result will be saved.
	-fcounter:  Keep count of the readen files
	'''

	# ----- Set the variables to be saved for the HS analysis -----

	var_list = ['eventID', 'energy', 'clockCount50', 'nhits', 'posx', 'posy', 'posz_av', 'runID']
	data_list = {obs: [] for obs in var_list}

	# ----- Reading -----

	fname_list = read_files_txt(file_txt_dir)
	fname_list.sort(key=natural_order)

	print(f'file name list: {fname_list}')

	for i_dx, fname_i in enumerate(fname_list):
		fdir_i = os.path.join(read_dir, fname_i)

		file_i = uproot.open(fdir_i)
		print(f'Keys of root file= {file_i.keys()} on iteration {i_dx}')
		data = file_i['output;1']

		# ==== Cut Conditions ====

		scintFit = np.array(data['scintFit'])
		fitValid = np.array(data['fitValid'])
		fit_condition = (scintFit & fitValid)

		posr_cut = 5500
		posr_av = np.array(data['posr_av'])
		posr_condition = (posr_av <= posr_cut)

		general_condition_cut = fit_condition & posr_condition

		# ==== Load, cut, and save the data ====

		for obs in var_list:
			data_obs = np.array(data[obs])
			data_list[obs].append(data_obs[general_condition_cut])

	print('concatenating data ...')
	data_dict = {obs: np.concatenate(data_list[obs]) for obs in var_list}
	del data_list


	# ====== HS Analysis ======

	hs_prompt_dict = {'eventID': [], 'runID':[]}
	hs_delay_dict = {'eventID': [], 'runID':[]}

	nhits = data_dict['nhits']
	time = data_dict['clockCount50']*20/1000 #microseconds
	posx = data_dict['posx']
	posy = data_dict['posy']
	posz = data_dict['posz_av']
	runID = data_dict['runID']
	eventID = data_dict['eventID']

	# ----- HS Cut Conditions -----

	dt_hs_remove = 1e6  # Remove 1s of Data after HS

	nhits_hs_number = 1500
	nhits_hs_condition = (nhits > nhits_hs_number)

	pos_condition_1 = (posx > 500) & (posx < 900) & (posy > -900) & (posy < -500) & (posz > 5300)
	pos_condition_2 = (posx > 0) & (posx < 500) & (posy > -1200) & (posy < -800) & (posz > 5300)
	pos_condition = (pos_condition_1) | (pos_condition_2)

	hs_condition = nhits_hs_condition & pos_condition

	# Save the HS prompt information
	hs_ev_index = np.where(hs_condition)[0]

	if len(hs_ev_index) > 0:
		print('HS Found. Saving prompt info. and followers')

		hs_prompt_dict['eventID'].extend(eventID[hs_ev_index])
		hs_prompt_dict['runID'].extend(runID[hs_ev_index])

		for hs_idx in hs_ev_index:
			follower_idx = hs_idx + 1

			hs_time = time[hs_idx]

			try:
				runID_follower = runID[follower_idx]
				eventID_follower = eventID[follower_idx]
				time_follower = time[follower_idx]

				dt = time_follower - hs_time

			except IndexError:
				continue

			while (dt > 0) and (dt <= dt_hs_remove):

				hs_delay_dict['runID'].append(runID_follower)
				hs_delay_dict['eventID'].append(eventID_follower)

				follower_idx += 1

				try:
					runID_follower = runID[follower_idx]
					eventID_follower = eventID[follower_idx]
					time_follower = time[follower_idx]
					dt = time_follower - hs_time

				except IndexError:
					break

		N_hs = len(hs_prompt_dict['eventID'])
		print(f'Nº of total Hotspot events found = {N_hs}')

		hs_prompt_dict = {k: np.array(v) for k, v in hs_prompt_dict.items()}
		hs_delay_dict = {k: np.array(v) for k, v in hs_delay_dict.items()}


		print('saving the dictionary in pickle format')
		with open(save_dir + f'hs_prompt_dict_{fcounter}.pkl', 'wb') as f:
			pickle.dump(hs_prompt_dict, f)

		with open(save_dir + f'hs_delay_dict_{fcounter}.pkl', 'wb') as f:
			pickle.dump(hs_delay_dict, f)

	else:
		print('No Hotspot were found')

	print('Done :)')

'''
# ====== Section to read only one file ======
if __name__ == '__main__':

	read_dir = '/share/neutrino/snoplus/Data/FullFill_2p2/rat_801/PPO/'
	file_txt_dir = '/lstore/sno/joankl/solar_analysis/real_data/2p2ppo/file_name_list/ntuple/sublist_0.txt'
	save_dir = '/lstore/sno/joankl/solar_analysis/real_data/2p2ppo/proof/'

	HS_analysis(read_dir, file_txt_dir, save_dir)
'''

# ====== Section to read various files and launch jobs ======

if __name__ == '__main__':

	data_type = "2.2PPO_HS_analysis"

	source_path = '/lstore/sno/joankl/.venv/bin/activate'
	read_dir = '/share/neutrino/snoplus/Data/FullFill_2p2/rat_801/PPO/'
	file_txt_dir = '/lstore/sno/joankl/solar_analysis/real_data/2p2ppo/file_name_list/ntuple/sublist_*.txt'
	save_dir = '/lstore/sno/joankl/solar_analysis/real_data/2p2ppo/ntuple_output/HS_results/'

	flist = glob.glob(file_txt_dir)
	flist = sorted(flist, key = natural_order)

	print('Ordered List:')
	print(' ')
	print(flist)

	os.makedirs(f'logs_{data_type}', exist_ok=True)

	for fi_dx, file_i in enumerate(flist):

		print(f'reding file with file {fi_dx}: {file_i}')

		script_name = f"run_analysis_{fi_dx}.sh"

		py_cmd = (
			f"from HS_analysis import HS_analysis; "
			f"HS_analysis ('{read_dir}', '{file_i}', '{save_dir}', fcounter = {fi_dx})"
			)
		
		script_cont = f"""#!/bin/bash
#SBATCH --job-name={data_type}_{fi_dx}
#SBATCH --output=logs_{data_type}/job_{fi_dx}.out
#SBATCH --error=logs_{data_type}/job_{fi_dx}.err
#SBATCH --partition=lipq
#SBATCH --mem=6G

echo "Running on host: $(hostname)"

source {source_path}

echo "Input file: {file_i}"

python -c "{py_cmd}"

echo "Job finished"
"""

		with open(script_name, "w") as f:
		    f.write(script_cont)

		# Enviar a la cola
		print(f"Enviando job {fi_dx}...")
		os.system(f"sbatch {script_name}")

		# Borrar el script generado para no llenar la carpeta de basura (opcional)
		# os.remove(script_name) 

		time.sleep(0.5) # Pequeña pausa para no saturar al scheduler


