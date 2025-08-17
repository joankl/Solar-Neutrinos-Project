"""
-logica:

usemos los resultados de Analysis.py para separar los hit types


"""

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

import uproot

def Chev_Scint_analysis(df): 

	"""
	:df -> type(Dataframe) = pandas.core.frame.DataFrame

	"""

	#Load and extract PMT info:
	file = uproot.open("/sno/electron bulk/simu_Analysis_elec_z_10MeV.root")
	pmt_info = file['pmt;1']

	#Extract PMT info.
	pmt_id = array(pmt_info['pmt_id'])
	pmt_pos_xyz = array(pmt_info['pmt_pos_xyz'])
	pmt_pos_sph = array(pmt_info['pmt_pos_sph'])

	#What I want to be the output:

	pmt_coordinates_xyz_hit1 = []
	pmt_coordinates_sph_hit1 = []

	pmt_coordinates_xyz_hit2 = []
	pmt_coordinates_sph_hit2 = []

	#Iterate over evtIDs:

	for evID in np.array(df['eventID']):

		evt_id_n = df.loc[df['eventID'] == evID]

		#Extract variables from Analysis_simulation:

		hit_type_ev = np.array(evt_id_n['hit type'])[0]
		pmt_ID_ev = np.array(evt_id_n['hitpmt ID'])[0]
		vertex_coord = np.array(evt_id_n['mc coordinates'])[0]
		time_res = np.array(evt_id_n['time residual'])[0]
		hitpmt_ID = np.array(evt_id_n['hitpmt ID'])[0]

		#Separate PMTid relative to the hit_type value:

		condition_1 = (hit_type_ev == 1)
		condition_2 = (hit_type_ev == 2)

		pmtID_hit_1 = np.array([], dtype = np.int64)      #ID of PMTs detecting hit_type = 1
		pmtID_hit_2 = np.array([], dtype = np.int64)      #ID of PMTs detecting hit_type = 2

		time_res_hit_1 = np.array([])                     #Time residuals of hit type = 1
		hitpmt_ID_hit_1 = np.array([], dtype = np.int64)  #ID of PMTS that recorded hit_type = 1

		time_res_hit_2 = np.array([])                     #Time residuals of hit type = 2
		hitpmt_ID_hit_2 = np.array([], dtype = np.int64)  #ID of PMTS that recorded hit_type = 2


		for (i,j) in zip(where(condition_1), where(condition_2)):

		    pmtID_hit_1 = np.append(pmtID_hit_1, pmt_ID_ev[i])
		    pmtID_hit_2 = np.append(pmtID_hit_2, pmt_ID_ev[j])
		    
		    time_res_hit_1 = np.append(time_res_hit_1, time_res[i])
		    hitpmt_ID_hit_1 = np.append(hitpmt_ID_hit_1, hitpmt_ID[i])
		    
		    time_res_hit_2 = np.append(time_res_hit_2, time_res[j])
		    hitpmt_ID_hit_2 = np.append(hitpmt_ID_hit_2, hitpmt_ID[j])


		#Now, find the PMT coordinates given the PMTid:

		condition_hit1 = np.in1d(pmt_id, hitpmt_ID_hit_1)
		condition_hit2 = np.in1d(pmt_id, hitpmt_ID_hit_2)

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

		pmt_coordinates_xyz_hit1.append(pmt_coordinates_xyz_hit1_ID)
		pmt_coordinates_sph_hit1.append(pmt_coordinates_sph_hit1_ID)

		pmt_coordinates_xyz_hit2.append(pmt_coordinates_xyz_hit2_ID)
		pmt_coordinates_sph_hit2.append(pmt_coordinates_sph_hit2_ID)


		#print(evID)
		#print(pmt_coordinates_xyz_hit1)

	# pmt_coordinates_xyz_hit1 = np.array(pmt_coordinates_xyz_hit1)
	# pmt_coordinates_sph_hit1 = np.array(pmt_coordinates_sph_hit1)

	# pmt_coordinates_xyz_hit2 = np.array(pmt_coordinates_xyz_hit2)
	# pmt_coordinates_sph_hit2 = np.array(pmt_coordinates_sph_hit2)


	# x_hit1 = pmt_coordinates_xyz_hit1[:,0]
	# y_hit1 = pmt_coordinates_xyz_hit1[:,1]
	# z_hit1 = pmt_coordinates_xyz_hit1[:,2]

	# zen_hit1 = pmt_coordinates_sph_hit1[:,0]
	# azim_hit1 = pmt_coordinates_sph_hit1[:,1]
	# rad_hit1 = pmt_coordinates_sph_hit1[:,2]

	# x_hit2 = pmt_coordinates_xyz_hit2[:,0]
	# y_hit2 = pmt_coordinates_xyz_hit2[:,1]
	# z_hit2 = pmt_coordinates_xyz_hit2[:,2]

	# zen_hit2 = pmt_coordinates_sph_hit2[:,0]
	# azim_hit2 = pmt_coordinates_sph_hit2[:,1]
	# rad_hit2 = pmt_coordinates_sph_hit2[:,3]

	#Construct DataFrame:

	data = {'xyz hit 1': pmt_coordinates_xyz_hit1,
			'spherical hit 1': pmt_coordinates_sph_hit1,
			'xyz hit 2': pmt_coordinates_xyz_hit2,
			'spherical hit 2': pmt_coordinates_sph_hit2,
			}

	# data = {'xyz hit 1': np.array([x_hit1, y_hit1, z_hit1]).tolist(),
	# 		'spherical hit 1': np.array([zen_hit1, azim_hit1, rad_hit1]).tolist(),
	# 		'xyz hit 2': np.array([x_hit2, y_hit2, z_hit2]).tolist(),
	# 		'spherical hit 2': np.array([zen_hit2, azim_hit2, rad_hit2]).tolist(),
	# 		}

	DataFrame_coord = pd.DataFrame(data)

	return DataFrame_coord



