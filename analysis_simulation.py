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
	:evID -> types(file) = list
	
	Funtion which recive the simu. file, and eventID and analize each ID in a given subset of simulation

	"""

	# Simulation Info
	data1 = file['T;7']

	#Extract Simulation Data
	evtid = np.array(data1['evtid'])
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

	for i_id in evID:

		#When simulation data includes many simulations the IDs are repeated, then execute the following to separate subsets:
		condition_1 = (evtid == i_id)
		index = where(condition_1)[0]               
		break_data = np.array([index[0]], dtype = np.int64)    #indices where data Breaks - Indices where subsets start
		# extract the index where data breaks
		for i in range(len(index)-1):
		    if index[i+1] != index[i]+1:
		    	break_data = np.append(break_data,index[i+1])
		        
		break_data = np.append(break_data,index[-1]+1)

		subset = [0,1]  #first subset

		init =  break_data[subset[0]]
		final_ = break_data[subset[1]]

		#Extract new data of subset: 
		sub_evtid = evtid[init:final_]
		sub_mc_position = mc_position[init:final_]
		sub_position = position[init:final_]
		sub_hit_pmtid = hit_pmtid[init:final_]
		sub_hit_pmttime = hit_pmttime[init:final_]
		sub_hit_residual = hit_residual[init:final_]
		sub_hit_type = hit_type[init:final_]

		#Retaking desired data of eventID  in subset
		condition_2 = (sub_evtid == i_id)
		pmtid_ev = np.extract(condition_2, sub_hit_pmtid) 			#Indices for which elements are True in eventID then extract the hit_pmtID
		time_res_ev = np.extract(condition_2, sub_hit_residual)
		pmttime_ev = np.extract(condition_2, sub_hit_pmttime)
		hit_type_ev = np.extract(condition_2, sub_hit_type)

		#Positions and radi
		mc_position_ev = []
		recons_position_ev = []
		mc_radius_ev = np.array([])
		recons_radius_ev = np.array([])

		for i in where(condition_2)[0]:
			mc_position_ev.append(sub_mc_position[i])
			recons_position_ev.append(sub_position[i])
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
		valid_pmt_id = np.array([],dtype = np.int64)	#valid PMTid

		valid_info_pmt_i = where(pmt_type == valid_type)[0]

		for i in valid_info_pmt_i:
			valid_pmt_id = np.append(valid_pmt_id, pmt_id[i])

		    
		# Buscar información valida dentro del eventos:
		#indices de info. valida
		valid_info_ev_i = where(np.in1d(pmtid_ev, valid_pmt_id))

		#información que deseo observar
		time_res_ev_valid = np.array([])
		pmt_id_ev_valid = np.array([], dtype = np.int64)
		hit_type_valid = np.array([], dtype = np.int64)

		for j in valid_info_ev_i:
			time_res_ev_valid = np.append(time_res_ev_valid, time_res_ev[j])
			pmt_id_ev_valid = np.append(pmt_id_ev_valid, pmtid_ev[j])
			hit_type_valid = np.append(hit_type_valid, hit_type_ev[j])

		#APPEND
		hit_pmtid_ID.append(pmt_id_ev_valid)
		hit_type_ID.append(hit_type_valid)
		mc_postion_vertex.append(mc_position_ev[0])

		#Save useful plots
		#matplotlib.use('Agg') #-----------> Dont show in screen Plots
		#hist(time_res_ev_valid, 100, 'time_res_ev', title = 'ID=' +str(i_id))

	#Construct DataFrame

	#AGREGAR COORDENADAS DE GENERACIÓN DEL EVENTO
	data = {'eventID': evID,
			'hitpmt ID': hit_pmtid_ID,
			'hit type': hit_type_ID,
			'mc coordinates': mc_postion_vertex,
			'mc radius': mc_rad_ID,
			'reconst radius': recons_rad_ID,
			'reconst error': reconstructed_error_ID,
			}

	DataFrame_analysis = pd.DataFrame(data)

		
	return DataFrame_analysis

# cuantos más motivos mejor :)
