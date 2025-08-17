# Useful Libraries --------------------------------

import uproot
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import seaborn as sn
from tqdm import tqdm

import numpy as np
from numpy import array, where, shape, reshape, pi

import pandas as pd

#-----------------------------------------------------------------------------------------------

def time_res_hist(df, bins = None, split = False, bin1 = None, bin2 = None):

	"""

	Function that construct histograms of given time residuals

	:df -> type(df) = pandas Data Frame (from analysis simualtion)
	:bins -> type(bins)=int (number of bins of histogram)

	"""


	evIDs = df['eventID'].tolist()

	for ID in evIDs:

		evt_id_n = df.loc[df['eventID'] == ID]

		time_residual = np.array(evt_id_n['time residual'])[0]

		if bins == None:
			bins = 2

		if split == False:

			x_title = 'time residual'
			title = 'Time Residual - evtID =' +str(ID)+'- 10MeV'

			plt.figure(figsize=(8,6))
			sn.histplot(time_residual, binwidth = bins, color = 'blue',edgecolor = 'black', alpha = 0.7)

			#plt.xlim(-50, 300)
			plt.yscale('log')


			plt.xlabel(x_title)
			plt.title(title)

			plt.savefig('Random point gen/figs/'+str(x_title) + str(title)+'.pdf', format = 'pdf')


		if split == True:

			time_residual_hit1 = (evt_id_n['time residual hit 1']).to_numpy()[0].tolist()	
			time_residual_hit2 = (evt_id_n['time residual hit 2']).to_numpy()[0].tolist()

			if bin1 == None:
				bin1 = 2

			if bin2 == None:
				bin2 = 2

			x_title = 'time residual'
			title = 'Compare Time Residual - evtID='+str(ID)+'-10MeV'

			plt.figure(figsize=(8,6))
			sn.histplot(time_residual_hit1, binwidth = bin1, color = 'blue', label = 'hit type 1')
			sn.histplot(time_residual_hit2, binwidth = bin2, color = 'red', alpha = 0.6, label = 'hit type 2')

			#plt.xlim(-50, 300)  #evtID = 3 REALMENTE TIENE UN RESIDUAL A -10 000! Por eso corto la ventana para jugar con el gráfico
			plt.yscale('log')

			plt.xlabel(x_title)
			plt.title(title)
			plt.legend()

			plt.savefig('figs/'+ str(title)+'.pdf', format = 'pdf')

#-----------------------------------------------------------------------------------------------		


def xy_plane_hits(df, bins = None, cut = False, inf_cut = None, up_cut = None):

	"""
	To use only when the electrons are directionated alogn z axis

	"""

	evIDs = df['eventID'].tolist()

	for ID in evIDs:

		evt_id_n = df.loc[df['eventID'] == ID]
		xyz_hit_1 = evt_id_n['xyz hit 1'].tolist()[0][0]

		x_hit_1 = xyz_hit_1[:, 0]
		y_hit_1 = xyz_hit_1[:, 1]
		z_hit_1 = xyz_hit_1[:, 2]

		if cut == False:

			title = 'xy plane - hit type 1 - evtID ='+str(ID)+'- 10MeV'

			if bins == None:
				bins = 30

			plt.hist2d(x = x_hit_1, y = y_hit_1, bins = [bins, bins], density = True)

			plt.ylabel('y_pmt')
			plt.xlabel('x_pmt')
			plt.title(title)

			#equal axis ration
			ax = plt.gca()
			ax.set_aspect('equal', adjustable='box')

			plt.savefig('figs/'+ str(title)+'.pdf', format = 'pdf')


		if cut == True:

			time_residual_hit1 = (evt_id_n['time residual hit 1']).to_numpy()[0].tolist()

			#Filter by means of time residual cut:

			#Data to be taken
			xyz_hit_1_cut = np.array([])                                    # cartesian coordinates
			#sph_hit_1_cut = np.array([])                                    # spherical coordinates

			#cuts
			for i in np.where((np.array(time_residual_hit1) > inf_cut) & (np.array(time_residual_hit1) < up_cut))[0]:

			    xyz_hit_1_cut = np.append(xyz_hit_1_cut, xyz_hit_1[i]).reshape((-1,3))
			    #sph_hit_1_cut = np.append(sph_hit_1_cut, sph_hit_1[i]).reshape((-1,3))
			
			x_hit_1_cut = xyz_hit_1_cut[:, 0]
			y_hit_1_cut = xyz_hit_1_cut[:, 1]
			z_hit_1_cut = xyz_hit_1_cut[:, 2]

			#Plot:

			if bins == None:
				bins = 30

			title = 'xy plane cut - hit type 1 - evtID='+str(ID)+'- 10MeV - cut'

			#plt.figure(figsize=(8,8))
			plt.hist2d(x = x_hit_1_cut, y = y_hit_1_cut, bins = [30,30], density = True)
			plt.ylabel('y_pmt')
			plt.xlabel('x_pmt')
			plt.title(title)

			#equal acis ration
			ax = plt.gca()
			ax.set_aspect('equal', adjustable='box')
			plt.show()

			plt.savefig('figs/'+ str(title)+'.pdf', format = 'pdf')


#-----------------------------------------------------------------------------------------------

def angular_dist(df, bins = None, coseno = False, cut = False, inf_cut = None, up_cut = None):

	"""
	angular distribution on azimuth(zenith) hits

	"""

	evIDs = df['eventID'].tolist()


	for ID in evIDs:

		evt_id_n = df.loc[df['eventID'] == ID]
		sph_hit_1 = evt_id_n['spherical hit 1'].tolist()[0][0]

		zenit_hit_1 = sph_hit_1[:,0]
		azimut_hit_1 = sph_hit_1[:,1]

		if bins == None:
			bins = 30

		if cut == False:

			if coseno == False:

				title = '$\phi(\Theta )$ angular distribution - hit type 1 - ID='+str(ID)+'- 10MeV'
				x_title = ' zenith '
				y_title = ' azimuth '

				plt.figure(figsize=(8,8))
				plt.hist2d(x = zenit_hit_1, y = azimut_hit_1, bins = [bins,bins], density = True)
				plt.xlabel(x_title)
				plt.ylabel(y_title)
				plt.title(title)

				plt.savefig('figs/'+ x_title + y_title+ 'ID=' +str( ID ) + '.pdf', format = 'pdf')

			if coseno == True:

				title = '$\phi(cos(\Theta))$ angular distribution - hit type 1 - evtID='+str(ID)+'- 10MeV'
				x_title = 'cos(zenith)'
				y_title = 'azimuth'

				plt.figure(figsize=(8,8))
				plt.hist2d(x = np.cos(zenit_hit_1), y = azimut_hit_1, bins = [bins,bins], density = True)
				plt.xlabel(x_title)
				plt.ylabel(y_title)
				plt.title(title)

				plt.savefig('figs/'+ x_title + y_title +'- ID='+str( ID ) + '.pdf', format = 'pdf')


		if cut == True:

			if coseno == False:

				time_residual_hit1 = (evt_id_n['time residual hit 1']).to_numpy()[0].tolist()
				sph_hit_1_cut = np.array([])

				for i in np.where((np.array(time_residual_hit1) > inf_cut) & (np.array(time_residual_hit1) < up_cut))[0]:

				    #xyz_hit_1_cut = np.append(xyz_hit_1_cut, xyz_hit_1[i]).reshape((-1,3))
				    sph_hit_1_cut = np.append(sph_hit_1_cut, sph_hit_1[i]).reshape((-1,3))

				zenit_hit_1_cut = sph_hit_1_cut[:,0]
				azimut_hit_1_cut = sph_hit_1_cut[:,1]
				rad_hit_1_cut = sph_hit_1_cut[:,2]

				if bins == None:
					bins = 30

				title = '$\phi(\Theta)$ angular distribution - hit type 1 - evtID='+str(ID)+'-10MeV - cut'
				x_title = ' zenith '
				y_title = ' azimuth '

				plt.figure(figsize=(8,8))
				plt.hist2d(x = zenit_hit_1_cut, y = azimut_hit_1_cut, bins = [bins,bins], density = True)
				plt.xlabel(x_title)
				plt.ylabel(y_title)
				plt.title(title)

				plt.savefig('figs/'+ x_title + y_title+' cut '+' evtID=' +str(ID) + '.pdf', format = 'pdf')

			if coseno == True:

				time_residual_hit1 = (evt_id_n['time residual hit 1']).to_numpy()[0].tolist()
				sph_hit_1_cut = np.array([])                                    				# spherical coordinates

				#cuts
				for i in np.where((np.array(time_residual_hit1) > inf_cut) & (np.array(time_residual_hit1) < up_cut))[0]:

				    #xyz_hit_1_cut = np.append(xyz_hit_1_cut, xyz_hit_1[i]).reshape((-1,3))
				    sph_hit_1_cut = np.append(sph_hit_1_cut, sph_hit_1[i]).reshape((-1,3))

				zenit_hit_1_cut = sph_hit_1_cut[:,0]
				azimut_hit_1_cut = sph_hit_1_cut[:,1]
				rad_hit_1_cut = sph_hit_1_cut[:,2]

				if bins == None:
					bins = 30

				title = '$\phi(\Theta)$ angular distribution - hit type 1 - evtID='+str(ID)+'-10MeV - cut'
				x_title = ' zenith '
				y_title = ' azimuth '

				plt.figure(figsize=(8,8))
				plt.hist2d(x = zenit_hit_1_cut, y = azimut_hit_1_cut, bins = [bins,bins], density = True)
				plt.xlabel(x_title)
				plt.ylabel(y_title)
				plt.title(title)

				plt.savefig('figs/'+ x_title + y_title+' cut '+' evtID=' +str(ID) + '.pdf', format = 'pdf')
