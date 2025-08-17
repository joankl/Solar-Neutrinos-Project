# Create here a set of opereations which simplifies to navegate by Data!!

#look for ID´s
#Clean PMT_type

import numpy as np
from numpy import array, where, shape, reshape, pi, cos, sin, sqrt, linspace

import pandas as pd

def GetPMTCoord(file, pmtID_ev):

	"""
	Generic function that given pmt IDs of an event extracts PMT Coordinates 

	"""

	# Extracts Information

	data1 = file['T;1']

	#PMT Information
	pmt_info = file['pmt;1']
	#Extract PMT info.
	pmt_id = array(pmt_info['pmt_id'])
	pmt_pos_xyz = array(pmt_info['pmt_pos_xyz'])
	pmt_pos_sph = array(pmt_info['pmt_pos_sph'])

	# Iterate over a range of pmtID and extract the coord. of every pmtID

	condition = np.in1d(pmt_id,pmtID_ev)
	xyz_coord_ev = []
	sphe_coord_ev = []

	for i in where(condition)[0]:
		xyz_coord_ev.append(pmt_pos_xyz[i])
		sphe_coord_ev.append(pmt_pos_sph[i])

	xyz_coord_ev = np.array(xyz_coord_ev)
	sphe_coord_ev = np.array(sphe_coord_ev)

	x = xyz_coord_ev[:,0]
	y = xyz_coord_ev[:,1]
	z = xyz_coord_ev[:,2]

	zen = sphe_coord_ev[:,0]
	azim = sphe_coord_ev[:,1]
	rad = sphe_coord_ev[:,2]

	data = {'x': x,
			'y': y,
			'z': z,
			'zenit': zen,
			'azimut': azim,
			'rad':rad,
			}

	DataFrame_coord = pd.DataFrame(data)

	return DataFrame_coord