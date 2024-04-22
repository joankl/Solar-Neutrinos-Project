	#Necessary library---------------------------------------------------------------------------
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import cm
import seaborn as sn
from tqdm import tqdm

import numpy as np
from numpy import array, where, shape, reshape

import pandas as pd

import itertools
#Useful functions for calculations and plots--------------------------------------------------

def magnitude(vector): 
	
	return np.sqrt(sum(pow(element, 2) for element in vector))

def hist(x, bins, xtitle, ytitle = None, title = None):
  
	plt.figure(figsize=(8,6))
	sn.histplot(x, bins = bins)
	plt.xlabel(xtitle)

	if ytitle == None:
		plt.ylabel('events')
	else:
		plt.ylabel(ytitle)
	  
	if title == None:
		plt.title(xtitle)
	else:
		plt.title(title)
	  
	plt.yscale('log')
	plt.savefig('figs/'+str(xtitle) + str(title)+'.pdf', format = 'pdf')

	#----------------------------------------------------------------------------------------------
	#Definitive function

def Analysis_simulation(file, evID = False):

	"""
	Parametros:
	:file -> type(file) = uproot.reading.ReadOnlyDirectory
	:evID -> type(evID) = list
	
	Funtion which recive the simu. file, and eventID and analize each ID in a given subset of simulation

	"""

	# Simulation Info
	data1 = file['T;7']

	#Extract Simulation Data
	mcID = np.array(data1['mcID'])
	energy = np.array(data1['energy'])
	mc_position = np.array(data1['mc_position'])
	position = np.array(data1['position'])
	hit_pmtid = np.array(data1['hit_pmtid'])
	hit_pmttime = np.array(data1['hit_pmttime'])
	time_residual = np.array(data1['hit_residual'])
	hit_type = np.array(data1['hit_type'])

	#PMT Information
	pmt_info = file['pmt;1']
	#Extract PMT info.
	pmt_id = array(pmt_info['pmt_id'])
	pmt_pos_xyz = array(pmt_info['pmt_pos_xyz'])
	pmt_pos_sph = array(pmt_info['pmt_pos_sph'])
	pmt_type = array(pmt_info['pmt_type'])


	#Separation of Jobs Info (subsets of run simulations) -----------------------------------------------------------

	break_i = [0]                                      # indices where the subsets start and end

	for i in range(len(mcID)-1):
		if mcID[i+1] < mcID[i]:  
			break_i.append(i)

	subsets = len(break_i) 

	#create variables that 
	for i in range(0, subsets):
		locals()['evtid_'+str(i)] = []
		locals()['mcID_'+str(i)] = []
		locals()['energy_'+str(i)] = []
		locals()['mc_position_'+str(i)] = []
		locals()['mc_momentum_'+str(i)] = []
		locals()['position_'+str(i)] = []
		locals()['momentum_'+str(i)] = []
		locals()['hit_pmtid_'+str(i)] = []
		locals()['hit_pmttime_'+str(i)] = []
		locals()['time_residual_'+str(i)] = []
		locals()['hit_type_'+str(i)] = []

		#proof if dividing the break_i
		locals()['break_i_'+str(i)] = []

	for i_dx, elem in enumerate(break_i):
		#i_dx -> iterable over break_i
		#print(i_dx)
		# start counting from 0 when extracting the subsets of info.
		if i_dx == 0:
			init = 0
			final = break_i[i_dx+1]

	        #print(init, final) (fine)

	        #check structure
			locals()['break_i_'+str(i_dx)].append(break_i[i_dx:i_dx+2])

			#add variables
			#locals()['evtid_'+str(i_dx)].append(mcID[init:final])
			locals()['mcID_'+str(i_dx)].append(mcID[init:final])
			locals()['energy_'+str(i_dx)].append(energy[init:final])
			locals()['mc_position_'+str(i_dx)].append(mc_position[init:final])	
			locals()['position_'+str(i_dx)].append(position[init:final])
			locals()['hit_pmtid_'+str(i_dx)].append(hit_pmtid[init:final])
			locals()['hit_pmttime_'+str(i_dx)].append(hit_pmttime[init:final])
			locals()['time_residual_'+str(i_dx)].append(time_residual[init:final])
			locals()['hit_type_'+str(i_dx)].append(hit_type[init:final])

			continue
	        
	    # extract middel subsets as being initiated by the next element in the before break_i and it ends in the next element of break_i
		if i_dx > 0 and i_dx < (len(break_i) - 2):
			#print(i_dx)
			init = break_i[i_dx] + 1
			final = break_i[i_dx+1]

	        #print(init, final)

	        #check structure
			locals()['break_i_'+str(i_dx)].append(break_i[i_dx:i_dx+2])

			#add variables
			#locals()['evtid_'+str(i_dx)].append(evtid[init:final])
			locals()['mcID_'+str(i_dx)].append(mcID[init:final])
			locals()['energy_'+str(i_dx)].append(energy[init:final])
			locals()['mc_position_'+str(i_dx)].append(mc_position[init:final])
			locals()['position_'+str(i_dx)].append(position[init:final])
			locals()['hit_pmtid_'+str(i_dx)].append(hit_pmtid[init:final])
			locals()['hit_pmttime_'+str(i_dx)].append(hit_pmttime[init:final])
			locals()['time_residual_'+str(i_dx)].append(time_residual[init:final])
			locals()['hit_type_'+str(i_dx)].append(hit_type[init:final])


		 # extract the last subset
		if i_dx == (len(break_i) - 2):
			#print(i_dx)
			init = break_i[i_dx] + 1
			final = break_i[i_dx+1]

	        #print(init, final) 

	        #check structure
			locals()['break_i_'+str(i_dx)].append(break_i[i_dx:i_dx+2])

	        #add variables
			#locals()['evtid_'+str(i_dx)].append(evtid[init:final])
			locals()['mcID_'+str(i_dx)].append(mcID[init:final])
			locals()['energy_'+str(i_dx)].append(energy[init:final])
			locals()['mc_position_'+str(i_dx)].append(mc_position[init:final])
			locals()['position_'+str(i_dx)].append(position[init:final])
			locals()['hit_pmtid_'+str(i_dx)].append(hit_pmtid[init:final])
			locals()['hit_pmttime_'+str(i_dx)].append(hit_pmttime[init:final])
			locals()['time_residual_'+str(i_dx)].append(time_residual[init:final])
			locals()['hit_type_'+str(i_dx)].append(hit_type[init:final])


	#Start analysis -------------------------------------------------------------------------------------------------
	#evID = np.array(evID) 
	job = np.arange(0, subsets - 1)

	job_number = []
	evID_ID = [] 											# Monte-Carlo ID
	energy_ID = []											# recons. Energy
	mc_rad_ID = []											# mc event radio UA!
	recons_rad_ID = []    									# reconstructed radio  UA!
	reconstructed_error_ID = []							# reconstruction error (magnitude of diference os position vectors) UA!
	hit_pmtid_ID = []									# ID of PMTs giving hits
	hit_type_ID = []									# hit_type of the eventID
	mc_postion_vertex = []								# coordinates of vertex generation
	time_res_ID = []									# (cleaned) time residual of the event ID
	pmt_xyz_ID = []										# pmt xyz coordintas for every ID
	pmt_sph_ID = []										# pmt sph coordintas for every ID

	#hit_type separation
	pmt_xyz_hit_1_ID = []
	pmt_xyz_hit_2_ID = []

	pmt_sph_hit_1_ID = []
	pmt_sph_hit_2_ID = []

	time_residual_hit1_ID = []							
	time_residual_hit2_ID = []							

	hitpmt_hit_1_ID = []							
	hitpmt_hit_2_ID = []


	pmt_ID_control = []

	for job_j in job:

		print('job:', job_j)

		#data to be used per job
		evID_job = locals()['mcID_'+str(job_j)][0]
		energy_job = locals()['energy_'+str(job_j)][0]
		mc_position_job = locals()['mc_position_'+str(job_j)][0]
		position_job = locals()['position_'+str(job_j)][0]
		hit_pmtid_job = locals()['hit_pmtid_'+str(job_j)][0]
		time_residual_job = locals()['time_residual_'+str(job_j)][0]
		hit_type_job = locals()['hit_type_'+str(job_j)][0]


		#DataFrame information to be saved - Is gonna be refereshed!:

		job_job = [] 											# Nº of job for each simulated data
		mc_rad_ID_job = np.array([])							# mc event radio
		recons_rad_ID_job = np.array([])    					# reconstructed radio
		reconstructed_error_ID_job = np.array([])				# reconstruction error (magnitude of diference os position vectors)
		hit_pmtid_ID_job = []									# ID of PMTs giving hits
		hit_type_ID_job = []									# hit_type of the eventID
		mc_postion_vertex_job = []								# coordinates of vertex generation 
		time_res_ID_job = []									# (cleaned) time residual of the event ID in job
		pmt_xyz_ID_job = []										# pmt xyz coordintas for every ID in job
		pmt_sph_ID_job = []										# pmt sph coordintas for every ID in job

		#hit_type separation
		pmt_xyz_hit_1_ID_job = []
		pmt_xyz_hit_2_ID_job = []

		pmt_sph_hit_1_ID_job = []
		pmt_sph_hit_2_ID_job = []

		time_residual_hit1_ID_job = []							
		time_residual_hit2_ID_job = []							

		hitpmt_hit_1_ID_job = []							
		hitpmt_hit_2_ID_job = []


		pmt_ID_control_job = []

		for i_id in np.unique(evID_job):

			#print('ID:', i_id)

			condition_ev = (evID_job == i_id)
			pmtid_ev = np.extract(condition_ev, hit_pmtid_job) 			#Indices for which elements are True in eventID then extract the hit_pmtID
			time_res_ev = np.extract(condition_ev, time_residual_job)
			#pmttime_ev = np.extract(condition_ev, hit_pmttime)
			hit_type_ev = np.extract(condition_ev, hit_type_job)

			#Positions and radi
			mc_position_ev = []
			recons_position_ev = []
			mc_radius_ev = np.array([])
			recons_radius_ev = np.array([])

			for i in where(condition_ev)[0]:
				mc_position_ev.append(mc_position_job[i])
				recons_position_ev.append(position_job[i])

			mc_position_ev = np.array(mc_position_ev)
			recons_position_ev = np.array(recons_position_ev)

			for i in range(shape(mc_position_ev)[0]):
				MC_radius = magnitude(mc_position_ev[i])
				RECONS_radius = magnitude(recons_position_ev[i])
				mc_radius_ev = np.append(mc_radius_ev, MC_radius)
				recons_radius_ev = np.append(recons_radius_ev, RECONS_radius)
			

			#APPEND ----------------------------------------------------------------
			recons_rad_ID_job = np.append(recons_rad_ID_job, recons_radius_ev[0])
			mc_rad_ID_job = np.append(mc_rad_ID_job, mc_radius_ev[0])
			# -----------------------------------------------------------------------

			#Calculate error of reconstruction
			mc_xyz = mc_position_ev[0]
			recons_xyz = recons_position_ev[0]
			dif_v = mc_xyz - recons_xyz											#vector diference
			error = magnitude(dif_v)											#recontruction error - a simple calculation
			#APPEND
			reconstructed_error_ID_job = np.append(reconstructed_error_ID_job, error)

			#PMT type selection
			valid_type = 1
			valid_pmt_id = []													#valid PMTid

			valid_info_pmt_i = where(pmt_type == valid_type)[0]

			for i in valid_info_pmt_i:
				valid_pmt_id.append(pmt_id[i])

			    
			# Buscar información valida dentro del eventos:
			# indices de info. valida
			valid_id_info_ev_i = where(np.in1d(pmtid_ev, valid_pmt_id))[0]

			#información que deseo observar
			valid_time_res_ev = []
			valid_pmt_id_ev =[]
			valid_hit_type_ev = []
			valid_xyz_pmt_ev = []
			valid_sph_pmt_ev = []

			for j in valid_id_info_ev_i:
				valid_time_res_ev.append(time_res_ev[j])
				valid_pmt_id_ev.append(pmtid_ev[j])
				valid_hit_type_ev.append(hit_type_ev[j])

			for k in valid_pmt_id_ev:
				valid_xyz_pmt_ev.append(pmt_pos_xyz[k])
				valid_sph_pmt_ev.append(pmt_pos_sph[k])


			#Separate hit types:

			index_hit_1 = where(np.in1d(valid_hit_type_ev, [1]))[0]
			index_hit_2 = where(np.in1d(valid_hit_type_ev, [2]))[0]

			pmt_xyz_hit_1_ev = []
			pmt_sph_hit_1_ev = []

			pmt_xyz_hit_2_ev = []
			pmt_sph_hit_2_ev = []

			time_residual_hit1_ev = []
			time_residual_hit2_ev = []

			hitpmt_ID_hit_1_ev = []
			hitpmt_ID_hit_2_ev = []


			for (i,j) in zip(index_hit_1, index_hit_2):

			    pmt_xyz_hit_1_ev.append(valid_xyz_pmt_ev[i])
			    pmt_sph_hit_1_ev.append(valid_sph_pmt_ev[i])
			    time_residual_hit1_ev.append(valid_time_res_ev[i])
			    hitpmt_ID_hit_1_ev.append(valid_pmt_id_ev[i])

			    pmt_xyz_hit_2_ev.append(valid_xyz_pmt_ev[j])
			    pmt_sph_hit_2_ev.append(valid_sph_pmt_ev[j])
			    time_residual_hit2_ev.append(valid_time_res_ev[j])
			    hitpmt_ID_hit_2_ev.append(valid_pmt_id_ev[j])

			#SUB_APPEND--- Is gonna be refersh for each job_j!----------------------------------------------------------------

			job_job.append(job_j)
			hit_pmtid_ID_job.append(valid_pmt_id_ev)
			hit_type_ID_job.append(valid_hit_type_ev)
			mc_postion_vertex_job.append(mc_position_ev[0])
			time_res_ID_job.append(valid_time_res_ev)
			pmt_xyz_ID_job.append(valid_xyz_pmt_ev)
			pmt_sph_ID_job.append(valid_sph_pmt_ev)

			pmt_xyz_hit_1_ID_job.append(pmt_xyz_hit_1_ev)
			pmt_sph_hit_1_ID_job.append(pmt_sph_hit_1_ev)
			time_residual_hit1_ID_job.append(time_residual_hit1_ev)
			hitpmt_hit_1_ID_job.append(hitpmt_ID_hit_1_ev)

			pmt_xyz_hit_2_ID_job.append(pmt_xyz_hit_2_ev)
			pmt_sph_hit_2_ID_job.append(pmt_sph_hit_2_ev)
			time_residual_hit2_ID_job.append(time_residual_hit2_ev)
			hitpmt_hit_2_ID_job.append(hitpmt_ID_hit_2_ev)


		evID_ID.append(evID_job)
		job_number.append(job_job)
		#job_number
		hit_pmtid_ID.append(hit_pmtid_ID_job)
		hit_type_ID.append(hit_type_ID_job)
		mc_postion_vertex.append(mc_postion_vertex_job)
		time_res_ID.append(time_res_ID_job)
		pmt_xyz_ID.append(pmt_xyz_ID_job)
		pmt_sph_ID.append(pmt_sph_ID_job)

		mc_rad_ID.append(mc_rad_ID_job)
		recons_rad_ID.append(recons_rad_ID_job)
		reconstructed_error_ID.append(reconstructed_error_ID_job)

		pmt_xyz_hit_1_ID.append(pmt_xyz_hit_1_ID_job)
		pmt_sph_hit_1_ID.append(pmt_sph_hit_1_ID_job)
		time_residual_hit1_ID.append(time_residual_hit1_ID_job)
		hitpmt_hit_1_ID.append(hitpmt_hit_1_ID_job)

		pmt_xyz_hit_2_ID.append(pmt_xyz_hit_2_ID_job)
		pmt_sph_hit_2_ID.append(pmt_sph_hit_2_ID_job)
		time_residual_hit2_ID.append(time_residual_hit2_ID_job)
		hitpmt_hit_2_ID.append(hitpmt_hit_2_ID_job)

	data = {'job': job,
		'eventID': evID_ID,
		'hitpmt ID': hit_pmtid_ID,
		'hit type': hit_type_ID,
		'time residual': time_res_ID,
		'mc coordinates': mc_postion_vertex,
		'mc radius': mc_rad_ID, 
		'reconst radius': recons_rad_ID, 
		'reconst error': reconstructed_error_ID, 
		'PMT xyz': pmt_xyz_ID,
		'PMT spherical': pmt_sph_ID,
		'hitpmt ID hit 1': hitpmt_hit_1_ID,
		'hitpmt ID hit 2': hitpmt_hit_2_ID,
		'time residual hit 1': time_residual_hit1_ID,
		'time residual hit 2': time_residual_hit2_ID,
		'xyz hit 1': pmt_xyz_hit_1_ID,
		'spherical hit 1': pmt_sph_hit_1_ID,
		'xyz hit 2': pmt_xyz_hit_2_ID,
		'spherical hit 2': pmt_sph_hit_2_ID,
		}

	DataFrame_analysis = pd.DataFrame(data)

		
	return DataFrame_analysis
