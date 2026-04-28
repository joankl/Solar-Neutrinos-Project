'''
Function designed to read ROOT files from RATROOT files.
This will output the results of the analysis in a numpy format.
The script basicly pass from root format to numoy format

Edits:
27/04/2026: Implement a reader to select the most recent TTree braches of Data. 
28/04/2026: Correction on cos_alpha calculation due to vertex correction on PMT position
'''

import uproot
import numpy as np 
import glob
import os 
import time

# Useful Function

def magnitude(vector): 
	#Function to Compute the radial position of events
	x = vector[:,0]
	y = vector[:,1]
	z = vector[:,2]

	r = np.sqrt(x**2 + y**2 + z**2)
	r = r.astype(np.float32)
	return r

def read_root(fin_dir, fout_dir, fcounter = 0):

	'''
	Function to read, perfomr cuts, and save observables of
	interest.

	Parameters:
	- fin_dir: directory + input file name (.root format)
	- fout_dir : output directory
	- fcounter: counter of the readen files. Useful when reading and
	            saving various files in order to distiguish the ouput
	            file names
	'''

	# ======== Data Cuts Settings ========
	energy_inf_cut = 2.5
	energy_sup_cut = 15

	posr_cut = 5500

	time_res_inf_cut = -1
	time_res_sup_cut = 50

	qhs_inf_cut = 2000
	qhs_sup_cut = 3400

	# ===================================

	load_data = uproot.open(fin_dir)

	#select the Key names for the PMT info and event data
	all_keys = load_data.keys()

	TTree_pmt_info_name = next((k for k in all_keys if k.startswith('pmt')), None) 

	t_keys = [k for k in all_keys if k.startswith('T;')]
	if t_keys:
		TTree_data_name = sorted(t_keys, key=lambda x: int(x.split(';')[1]))[-1]
	else:
		# If the key has no cicle, then pick the key 'T' which select by default the most recent key verion
		TTree_data_name = 'T'

	event_data = load_data[TTree_data_name]
	pmt_data = load_data[TTree_pmt_info_name]

	#event info to be used:
	var_event_list = ['evtid', 'energy_corr', 'position', 'momentum_mc',
	                  'hit_pmtid', 'hit_residual']  #list the name of the varibles to be extracted and used for the solarnu analysis.

	#pmt info to be used
	var_pmt_list = ['pmt_id', 'pmt_pos_xyz', 'pmt_type']

	#Observables to save
	var_name_save_list = ['evtid', 'energy_corr', 'posr', 'hit_residual']

	# Create empty lists to save observables thar are not part of the three, or are coordinated variables 
	multi_cos_alpha = []   # cos_alpha for the multiple PMTs record
	multi_hit_pmt_xyz = [] # hit pmt coordinates
	multi_mc_momentum = [] # simulated event direction
	n_init_evs = []        # Initial Nº of Solar Nue MC events

	# Extract the variables with the name of the var_event_list in numpy.array from the .root file
	observables = {}

	# Load variables form the Event Info TTree Branch
	for var_name_i in var_event_list:
		observables[var_name_i] = np.array(event_data[var_name_i])

	# ============= Count the initial Nº of events =============
	evtid = observables['evtid']
	evIDi_unique = []  #Empty list to be filled with the unique (non-redundant) values of the initial evIDs. len(evIDi_unique) = Nº of events

	N_data = len(evtid)
	indices_to_delete = [] #Indiced to remove from the observable arrays due to nhits cut

	# ----- Calculation of Data break due to repeated non-perhit observales -----
	data_break_mask = (np.diff(evtid) != 0)
	data_break_index = np.where(data_break_mask)[0]+1
	data_break_index = np.concatenate(([0], data_break_index, [len(evtid) - 1])) # Add the initial and last index

	N_terms = len(data_break_index)

	for i_dx in range(N_terms - 2):
		init_i = data_break_index[i_dx]
		final_i = data_break_index[i_dx + 1]
		evIDi_unique.append(evtid[init_i: final_i][0])

	N_init_evs = len(evIDi_unique)
	print(f'Nº of initial events: {N_init_evs}')
	n_init_evs.append(N_init_evs)

	# Compute posr
	observables['posr'] = magnitude(observables['position'])

	# Load the PMT info TTree Branch
	for var_name_i in var_pmt_list:
		observables[var_name_i] = np.array(pmt_data[var_name_i])

	# Filtering of valid PMT id through PMT type
	pmt_type_condition = (observables['pmt_type'] == 1)
	pmt_id_valid = observables['pmt_id'][pmt_type_condition]

	# general cut conditions
	hit_pmt_id_condition = np.in1d(observables['hit_pmtid'], pmt_id_valid)
	energy_condition = (observables['energy_corr'] >= energy_inf_cut) & (observables['energy_corr'] <= energy_sup_cut)
	time_res_condition = (observables['hit_residual'] >= time_res_inf_cut) & (observables['hit_residual'] <= time_res_sup_cut)
	#qhs_condition = (observables['hit_pmtQHS'] >= qhs_inf_cut) & (observables['hit_pmtQHS'] <= qhs_sup_cut)
	posr_condition = (observables['posr'] <= posr_cut)

	general_condition = hit_pmt_id_condition & energy_condition & time_res_condition & posr_condition

	# Apply the general cut conditions to observables
	for var_name_i in var_event_list:
		observables[var_name_i] = observables[var_name_i][np.array(general_condition)]
	observables['posr'] = observables['posr'][general_condition]

	print(f'saving observables {var_name_save_list}')
	for var_name_i in  var_name_save_list:
		np.save(fout_dir + var_name_i + f'_{fcounter}.npy', observables[var_name_i])

	# ========== cos_alpha computation ==========

	N_samples = len(observables['hit_residual'])

	print('Computing cos_alpha')

	sun_dir = observables['momentum_mc']
	pmt_hit_id = observables['hit_pmtid']
	pmt_hit_pos_xyz = observables['pmt_pos_xyz'][pmt_hit_id]
	event_pos_xyz = observables['position']

	# Correction due to reconstructed vertex position
	photon_dir = pmt_hit_pos_xyz - event_pos_xyz

	norm_sun = np.linalg.norm(sun_dir, axis=1, keepdims=True)
	norm_photon = np.linalg.norm(photon_dir, axis=1, keepdims=True)

	sun_dir_unit = (sun_dir / norm_sun).astype(np.float32)
	photon_dir_unit = (photon_dir / norm_pmt).astype(np.float32)

	multi_cos_alpha = np.sum(sun_dir_unit * photon_dir_unit, axis=1)

#	for sample_idx in range(N_samples):
#
#		sun_dir = observables['momentum_mc'][sample_idx]
#		pmt_hit_id = observables['hit_pmtid'][sample_idx]
#		pmt_hit_xyz = observables['pmt_pos_xyz'][pmt_hit_id]
#
#		norm1 = np.linalg.norm(sun_dir)
#		norm2 = np.linalg.norm(pmt_hit_xyz)
#
#		sun_dir = sun_dir / norm1
#		sun_dir = sun_dir.astype(np.float32)
#
#		pmt_hit_xyz = pmt_hit_xyz / norm2
#		pmt_hit_xyz = pmt_hit_xyz.astype(np.float32)
#
#		dot_prod = np.dot(sun_dir, pmt_hit_xyz)
#		cos_alpha = dot_prod
#
#		multi_hit_pmt_xyz.append(pmt_hit_xyz)
#		multi_mc_momentum.append(sun_dir)
#		multi_cos_alpha.append(cos_alpha)
#
#	multi_cos_alpha = np.array(multi_cos_alpha)
	#multi_hit_pmt_xyz = np.array(multi_hit_pmt_xyz)
	#multi_mc_momentum = np.array(multi_mc_momentum)

	#print('saving hit PMT xyz coordinates')
	#np.save(fout_dir + '_hitpmt_xyz' + f'_{fcounter}.npy', multi_hit_pmt_xyz)

	#print('saving event direction')
	#np.save(fout_dir + '_mc_momentum' + f'_{fcounter}.npy', multi_mc_momentum)

	print('saving cos_alpha')
	np.save(fout_dir + 'cos_alpha' + f'_{fcounter}.npy', multi_cos_alpha)

	print('Analysis Done!')


if __name__ == '__main__':

	data_type = "8B_Nue_MC_bisMSB"

	source_path = '/lstore/sno/joankl/.venv/bin/activate'
	fin_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/bisMSB/B8_solar_Nue/ratDS_output/root_files/*.root'
	fout_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/bisMSB/B8_solar_Nue/ratDS_output/np_files/'
	flist = glob.glob(fin_dir)

	os.makedirs(f'logs_{data_type}', exist_ok=True)

	for fi_dx, file_i in enumerate(flist):

		print(f'reding file with file {fi_dx}: {file_i}')

		script_name = f"run_analysis_{fi_dx}.sh"

		py_cmd = (
			f"from analysis_fROOT import read_root; "
			f"read_root ('{file_i}', '{fout_dir}', fcounter = {fi_dx})"
			)
		
		script_cont = f"""#!/bin/bash
#SBATCH --job-name={data_type}_{fi_dx}
#SBATCH --output=logs_{data_type}/job_{fi_dx}.out
#SBATCH --error=logs_{data_type}/job_{fi_dx}.err
#SBATCH --partition=lipq
#SBATCH --mem=12G

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

		
