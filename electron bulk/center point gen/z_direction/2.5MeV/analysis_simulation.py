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

def Analysis_simulation(file, evID):

	"""
	Parametros:
	:file -> type(file) = uproot.reading.ReadOnlyDirectory
	:evID -> type(evID) = list
	
	Funtion which recive the simu. file, and eventID and analize each ID in a given subset of simulation

	"""

	# Simulation Info
	data1 = file['T;1']

	#Extract Simulation Data
	mcID = np.array(data1['mcID'])
	mc_position = np.array(data1['mc_position'])
	position = np.array(data1['position'])
	hit_pmtid = np.array(data1['hit_pmtid'])
	hit_pmttime = np.array(data1['hit_pmttime'])
	hit_residual = np.array(data1['hit_residual'])
	hit_type = np.array(data1['hit_type'])

	#PMT Information
	pmt_info = file['pmt;1']
	#Extract PMT info.
	pmt_id = array(pmt_info['pmt_id'])
	pmt_pos_xyz = array(pmt_info['pmt_pos_xyz'])
	pmt_pos_sph = array(pmt_info['pmt_pos_sph'])
	pmt_type = array(pmt_info['pmt_type'])


	#Start analysis -------------------------------------------------------------------------------------------------
	evID = np.array(evID)

	#DataFrame information to be saved:
	#ID
	mc_rad_ID = np.array([])							# mc event radio
	recons_rad_ID = np.array([])    					# reconstructed radio
	reconstructed_error_ID = np.array([])				# reconstruction error (magnitude of diference os position vectors)
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

	for i_id in evID:

		condition_ev = (mcID == i_id)
		pmtid_ev = np.extract(condition_ev, hit_pmtid) 			#Indices for which elements are True in eventID then extract the hit_pmtID
		time_res_ev = np.extract(condition_ev, hit_residual)
		pmttime_ev = np.extract(condition_ev, hit_pmttime)
		hit_type_ev = np.extract(condition_ev, hit_type)

		#Positions and radi
		mc_position_ev = []
		recons_position_ev = []
		mc_radius_ev = np.array([])
		recons_radius_ev = np.array([])

		for i in where(condition_ev)[0]:
			mc_position_ev.append(mc_position[i])
			recons_position_ev.append(position[i])

		mc_position_ev = np.array(mc_position_ev)
		recons_position_ev = np.array(recons_position_ev)

		for i in range(shape(mc_position_ev)[0]):
			MC_radius = magnitude(mc_position_ev[i])
			RECONS_radius = magnitude(recons_position_ev[i])
			mc_radius_ev = np.append(mc_radius_ev, MC_radius)
			recons_radius_ev = np.append(recons_radius_ev, RECONS_radius)
		

		#APPEND
		recons_rad_ID = np.append(recons_rad_ID, recons_radius_ev[0])
		mc_rad_ID = np.append(mc_rad_ID, mc_radius_ev[0])

		#Calculate error of reconstruction
		mc_xyz = mc_position_ev[0]
		recons_xyz = recons_position_ev[0]
		dif_v = mc_xyz - recons_xyz											#vector diference
		error = magnitude(dif_v)											#recontruction error - a simple calculation
		#APPEND
		reconstructed_error_ID = np.append(reconstructed_error_ID, error)

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

		#APPEND-------------------------------------------------------------------------------------------
		hit_pmtid_ID.append(valid_pmt_id_ev)
		hit_type_ID.append(valid_hit_type_ev)
		mc_postion_vertex.append(mc_position_ev[0])
		time_res_ID.append(valid_time_res_ev)
		pmt_xyz_ID.append(valid_xyz_pmt_ev)
		pmt_sph_ID.append(valid_sph_pmt_ev)

		pmt_xyz_hit_1_ID.append(pmt_xyz_hit_1_ev)
		pmt_sph_hit_1_ID.append(pmt_sph_hit_1_ev)
		time_residual_hit1_ID.append(time_residual_hit1_ev)
		hitpmt_hit_1_ID.append(hitpmt_ID_hit_1_ev)

		pmt_xyz_hit_2_ID.append(pmt_xyz_hit_2_ev)
		pmt_sph_hit_2_ID.append(pmt_sph_hit_2_ev)
		time_residual_hit2_ID.append(time_residual_hit2_ev)
		hitpmt_hit_2_ID.append(hitpmt_ID_hit_2_ev)


	data = {'eventID': evID,
		'hitpmt ID': hit_pmtid_ID,
		'hit type': hit_type_ID,
		'time residual': time_res_ID,
		'mc coordinates': mc_postion_vertex,
		'mc radius': mc_rad_ID, #P
		'reconst radius': recons_rad_ID, #P
		'reconst error': reconstructed_error_ID, #P
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
