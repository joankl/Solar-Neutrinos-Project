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


	pmt_ID_control = []

	for i_id in evID:

		condition_2 = (mcID == i_id)
		pmtid_ev = np.extract(condition_2, hit_pmtid) 			#Indices for which elements are True in eventID then extract the hit_pmtID
		time_res_ev = np.extract(condition_2, hit_residual)
		pmttime_ev = np.extract(condition_2, hit_pmttime)
		hit_type_ev = np.extract(condition_2, hit_type)

		#Positions and radi
		mc_position_ev = []
		recons_position_ev = []
		mc_radius_ev = np.array([])
		recons_radius_ev = np.array([])

		for i in where(condition_2)[0]:
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

		#APPEND-------------------------------------------------------------------------------------------
		hit_pmtid_ID.append(pmt_id_ev_valid)
		hit_type_ID.append(hit_type_valid)
		mc_postion_vertex.append(mc_position_ev[0])
		time_res_ID.append(time_res_ev_valid)

		#ORGANIZE DATA IN A DATA FRAME TO BE USED BELLOW

	pre_data = {'eventID': evID,
		'hitpmt ID': hit_pmtid_ID,
		'hit type': hit_type_ID,
		'mc coordinates': mc_postion_vertex,
		'mc radius': mc_rad_ID,
		'reconst radius': recons_rad_ID,
		'reconst error': reconstructed_error_ID,
		'time residual': time_res_ID,
		}

	df = pd.DataFrame(pre_data)

    #Now, separate the hit_type data: ------------------------------------------------------------------

    #Information to be saved
	pmt_coordinates_xyz_hit1 = []					#Coordinates of pmts giving hit type = 1
	pmt_coordinates_sph_hit1 = []

	pmt_coordinates_xyz_hit2 = []					#Coordinates of pmts giving hit type = 2
	pmt_coordinates_sph_hit2 = []

	time_residual_hit1 = []							#Time residual of hit type = 1
	time_residual_hit2 = []							#Time residual of hit type = 2

	hitpmt_ID_hit_1 = []							#ID of pmt giving hit type = 1
	hitpmt_ID_hit_2 = []

	for evID in np.array(df['eventID']):

		evt_id_n = df.loc[df['eventID'] == evID]

		#Extract variables of interest from pre_data

		hit_type_ev = np.array(evt_id_n['hit type'])[0]
		pmt_ID_ev = np.array(evt_id_n['hitpmt ID'])[0]
		vertex_coord = np.array(evt_id_n['mc coordinates'])[0]
		time_res = np.array(evt_id_n['time residual'])[0]
		hitpmt_ID = np.array(evt_id_n['hitpmt ID'])[0]

		condition_hitype_1 = (hit_type_ev == 1)
		condition_hitype_2 = (hit_type_ev == 2)

		pmtID_hit_1 = np.array([], dtype = np.int64)      #ID of PMTs detecting hit_type = 1
		pmtID_hit_2 = np.array([], dtype = np.int64)      #ID of PMTs detecting hit_type = 2

		time_res_hit_1_ID_ev = np.array([])                     #Time residuals of hit type = 1
		hitpmt_ID_hit_1_ev = np.array([], dtype = np.int64)  #ID of PMTS that recorded hit_type = 1

		time_res_hit_2_ID_ev = np.array([])                     #Time residuals of hit type = 2
		hitpmt_ID_hit_2_ev = np.array([], dtype = np.int64)  #ID of PMTS that recorded hit_type = 2


		for (i,j) in zip(where(condition_hitype_1), where(condition_hitype_2)):

		    pmtID_hit_1 = np.append(pmtID_hit_1, pmt_ID_ev[i])
		    pmtID_hit_2 = np.append(pmtID_hit_2, pmt_ID_ev[j])
		    
		    time_res_hit_1_ID_ev = np.append(time_res_hit_1_ID_ev, time_res[i])
		    hitpmt_ID_hit_1_ev = np.append(pmtID_hit_1, hitpmt_ID[i])
		    
		    time_res_hit_2_ID_ev = np.append(time_res_hit_2_ID_ev, time_res[j])
		    hitpmt_ID_hit_2_ev = np.append(pmtID_hit_2, hitpmt_ID[j])

		pmt_ID_control.append(hitpmt_ID_hit_2_ev)

		hitpmt_ID_hit_1.append(hitpmt_ID_hit_1_ev)
		hitpmt_ID_hit_2.append(hitpmt_ID_hit_2_ev)


		time_residual_hit1.append(time_res_hit_1_ID_ev)
		time_residual_hit2.append(time_res_hit_2_ID_ev)


		#Now, find the PMT coordinates given the PMTid:

		condition_hit1 = np.in1d(pmt_id, hitpmt_ID_hit_1_ev)
		condition_hit2 = np.in1d(pmt_id, hitpmt_ID_hit_2_ev)


		pmt_coordinates_xyz_hit1_ID = []
		pmt_coordinates_sph_hit1_ID = []

		pmt_coordinates_xyz_hit2_ID = []
		pmt_coordinates_sph_hit2_ID = []


		for (i,j) in zip(where(condition_hit1), where(condition_hit2)):

			#xyz_coord_ev.append(pmt_pos_xyz[i])
			#sphe_coord_ev.append(pmt_pos_sph[i])
			pmt_coordinates_xyz_hit1_ID.append(pmt_pos_xyz[i])
			pmt_coordinates_sph_hit1_ID.append(pmt_pos_sph[i])

			pmt_coordinates_xyz_hit2_ID.append(pmt_pos_xyz[j])
			pmt_coordinates_sph_hit2_ID.append(pmt_pos_sph[j])


		#APPEND-------------------------------------------------------------------------------------------
		pmt_coordinates_xyz_hit1.append(pmt_coordinates_xyz_hit1_ID)
		pmt_coordinates_sph_hit1.append(pmt_coordinates_sph_hit1_ID)

		pmt_coordinates_xyz_hit2.append(pmt_coordinates_xyz_hit2_ID)
		pmt_coordinates_sph_hit2.append(pmt_coordinates_sph_hit2_ID)

		#Save useful plots: Choose a variable and then plot it - by default is time_res_ev_valid

		#matplotlib.use('Agg') #-----------> Dont show Plots in screen 
		#hist(time_res_ev_valid, 100, 'time_res_ev', title = 'ID=' +str(i_id))

		#cut in time_residual
		#hit_pmtid_ID_cut = 

	#Construct DataFrame

	#AGREGAR COORDENADAS DE GENERACIÓN DEL EVENTO
	data = {'eventID': df['eventID'].tolist(),
			'hitpmt ID': hit_pmtid_ID,
			'hit type': hit_type_ID,
			'time residual': time_res_ID,
			'mc coordinates': mc_postion_vertex,
			'mc radius': mc_rad_ID,
			'reconst radius': recons_rad_ID,
			'reconst error': reconstructed_error_ID,
			'hitpmt ID hit 1': hitpmt_ID_hit_1,
			'hitpmt ID hit 2': hitpmt_ID_hit_2,
			'time residual hit 1': time_residual_hit1,
			'time residual hit 2': time_residual_hit2,
			'xyz hit 1': pmt_coordinates_xyz_hit1,
			'spherical hit 1': pmt_coordinates_sph_hit1,
			'xyz hit 2': pmt_coordinates_xyz_hit2,
			'spherical hit 2': pmt_coordinates_sph_hit2,
			}

	DataFrame_analysis = pd.DataFrame(data)

		
	return DataFrame_analysis

# cuantos más motivos mejor :)
