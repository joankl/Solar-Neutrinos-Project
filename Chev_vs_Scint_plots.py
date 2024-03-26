# Useful Libraries --------------------------------

import uproot
import numpy as np

import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import style
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

			plt.savefig('figs/'+str(x_title) + str(title)+'.pdf', format = 'pdf')


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

			plt.figure(figsize=(8,6))
			sn.set_style(rc = {'axes.facecolor': 'black'})
			sn.histplot(x = x_hit_1, y = y_hit_1, bins = [30,30], stat='count', cbar = 'True', cmap = cm.nipy_spectral)

			plt.ylabel('y_pmt')
			plt.xlabel('x_pmt')
			plt.title(title)

			#equal axis ration
			ax = plt.gca()
			ax.set_aspect('equal', adjustable='box')

			plt.savefig('figs/'+ str(title)+'.pdf', format = 'pdf')
			plt.show()


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

			plt.figure(figsize=(8,6))
			sn.set_style(rc = {'axes.facecolor': 'black'})
			sn.histplot(x = x_hit_1_cut, y = y_hit_1_cut, bins = [30,30], stat='count', cbar = 'True', cmap = cm.nipy_spectral)
			plt.ylabel('y_pmt')
			plt.xlabel('x_pmt')
			plt.title(title)

			#equal acis ration
			ax = plt.gca()
			ax.set_aspect('equal', adjustable='box')

			plt.savefig('figs/'+ str(title)+'.pdf', format = 'pdf')
			plt.show()


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

				plt.figure(figsize=(8,6))
				sn.set_style(rc = {'axes.facecolor': 'black'})
				sn.histplot(x = zenit_hit_1, y = azimut_hit_1, bins = [bins,bins] , stat='count', cbar = 'True', cmap = cm.nipy_spectral)
				plt.xlabel(x_title)
				plt.ylabel(y_title)
				plt.title(title)

				plt.savefig('figs/'+ x_title + y_title+ 'ID=' +str( ID ) + '.pdf', format = 'pdf')

			if coseno == True:

				title = '$\phi(cos(\Theta))$ angular distribution - hit type 1 - evtID='+str(ID)+'- 10MeV'
				x_title = 'cos(zenith)'
				y_title = 'azimuth'

				plt.figure(figsize=(8,6))
				sn.set_style(rc = {'axes.facecolor': 'black'})
				sn.histplot(x = np.cos(zenit_hit_1), y = azimut_hit_1, bins = [bins,bins] , stat='count', cbar = 'True', cmap = cm.nipy_spectral)
				plt.xlabel(x_title)
				plt.ylabel(y_title)
				plt.title(title)

				plt.savefig('figs/'+ x_title + y_title +'- ID='+str( ID ) + '.pdf', format = 'pdf')
				plt.show()


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

				plt.figure(figsize=(8,6))
				sn.set_style(rc = {'axes.facecolor': 'black'})
				sn.histplot(x = zenit_hit_1_cut, y = azimut_hit_1_cut, bins = [bins,bins] , stat='count', cbar = 'True', cmap = cm.nipy_spectral)
				plt.xlabel(x_title)
				plt.ylabel(y_title)
				plt.title(title)

				plt.savefig('figs/'+ x_title + y_title+' cut '+' evtID=' +str(ID) + '.pdf', format = 'pdf')
				plt.show()

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

				title = '$\phi(cos(\Theta))$ angular distribution - hit type 1 - evtID='+str(ID)+'-10MeV - cut'
				x_title = ' cos(zenith) '
				y_title = ' azimuth '

				plt.figure(figsize=(8,6))
				sn.set_style(rc = {'axes.facecolor': 'black'})
				sn.histplot(x = np.cos(zenit_hit_1_cut), y = azimut_hit_1_cut, bins = [bins,bins] , stat='count', cbar = 'True',cmap = cm.nipy_spectral)
				plt.xlabel(x_title)
				plt.ylabel(y_title)
				plt.title(title)

				plt.savefig('figs/'+ x_title + y_title+ ' cut '+' evtID=' +str(ID) + '.pdf', format = 'pdf')
				plt.show()

#----------------------------------------------------------------------------------------------------

def angle_hit_direction(df, bins = None, cut = False, inf_cut = None, up_cut = None):

	"""

	Function that return the 2D histogram of the time_residual distribution  vs. the cos of angle
	between the directon of the event and the PMT hits

	IF DIMENSIONAL ERROR MEANS A PROBLEM WITH THE COSINE AND THE TIME RESIDUAL SHAPES OF AN EVENT!

	"""

	evIDs = df['eventID'].tolist()

	for ID in evIDs:

		evt_id_n = df.loc[df['eventID'] == ID]

		#useful variables bellow:

		time_residual = (evt_id_n['time residual']).to_numpy()[0]
		xyz_hit = (evt_id_n['PMT xyz']).to_numpy()[0]

		#xyz_hit_1 = (evt_id_n['xyz hit 1']).to_numpy()[0][0]


		#To compute: cos(angle)--------------------------
		cos_angle = np.array([])

		vec_ev = np.array([0.0, 0.0, -1.0]) #----> Direction of event!

		N = shape(xyz_hit)[0] 
		for i in range(N):
			cos_val = np.dot(xyz_hit[i],vec_ev)/np.linalg.norm(xyz_hit[i])
			cos_angle = np.append(cos_angle, cos_val)
		# -----------------------------------------------

		if bins == None:
			bins = 30

		if cut == False:

			title = 'cos(α)(time res.) - evtID =' + str(ID) + ' - 10MeV'
			x_title = 'cos(α)'
			y_title = 'time residual'

			plt.figure(figsize=(8,8))
			sn.set_style(rc = {'axes.facecolor': 'black'})

			sn.histplot(x = cos_angle, y = time_residual, bins = [bins,bins], stat='count', cbar = 'True', cmap = cm.nipy_spectral)
			plt.xlabel(x_title)
			plt.ylabel(y_title)
			plt.title(title)

			plt.savefig('figs/All hit Type/'+ title + ' evtID=' +str(ID) + '.jpg', format = 'jpg')
			plt.show()


		if cut == True:

			#Data to be taken
			time_residual_cut = np.array([]) 
			cos_angle_cut = np.array([])                                  

			#cuts
			for i in np.where((np.array(time_residual) > inf_cut) & (np.array(time_residual) < up_cut))[0]:
				time_residual_cut = np.append(time_residual_cut, time_residual[i])
				cos_angle_cut = np.append(cos_angle_cut, cos_angle[i])


			title = 'cos(α)(time res.) cut - evtID =' + str(ID) + ' - 10MeV'
			x_title = 'cos(α)'
			y_title = 'time residual'

			plt.figure(figsize=(8,8))
			sn.set_style(rc = {'axes.facecolor': 'black'})

			sn.histplot(x = cos_angle_cut, y = time_residual_cut, bins = [bins,bins], stat='count', cbar = 'True', cmap = cm.nipy_spectral)
			plt.xlabel(x_title)
			plt.ylabel(y_title)
			plt.title(title)

			plt.savefig('figs/All hit Type/'+ title + ' evtID=' +str(ID) + '.jpg', format = 'jpg')
			plt.show()

# ------------------------------------------------------------------------------------------

def angle_hit_direction3D(df, bins = None, cut = False, inf_cut = None, up_cut = None):

	"""
	Function that return a 3D histogram with time res and cos of angle
	between the directon of the event and the PMT hits in the xy plane
	and counts in the z axis

	"""

	evIDs = df['eventID'].tolist()

	for ID in evIDs:

		evt_id_n = df.loc[df['eventID'] == ID]

		#useful variables

		time_residual = (evt_id_n['time residual']).to_numpy()[0]
		xyz_hit = (evt_id_n['PMT xyz']).to_numpy()[0]


		time_residual_hit1 = (evt_id_n['time residual hit 1']).to_numpy()[0]
		xyz_hit_1 = (evt_id_n['xyz hit 1']).to_numpy()[0][0]


		#To compute: cos(angle)--------------------------
		cos_angle = np.array([])

		vec_ev = np.array([0.0, 0.0, -1.0]) #----> Direction of event!

		N = shape(xyz_hit)[0] 
		for i in range(N):
			cos_val = np.dot(xyz_hit[i],vec_ev)/np.linalg.norm(xyz_hit[i])
			cos_angle = np.append(cos_angle, cos_val)
		# -----------------------------------------------

		if bins == None:
			bins = 30

		#to resets default plot style
		mpl.rcParams.update(mpl.rcParamsDefault)

		if cut == False:
			#extract count of events
			counts, x_bin_edges, y_bin_edges = np.histogram2d(x = cos_angle, y = time_residual_hit1, bins = [bins, bins]) #PROBLEMA! No me deja usar nº diferente de bins


			#Obtain coordinates of bins ---------------------------------
			x_bin = []
			y_bin = []

			Nx = len(x_bin_edges)
			Ny = len(y_bin_edges)

			for (x_i, y_i) in zip(range(Nx-1), range(Ny-1)):
				mid_point_x = (x_bin_edges[x_i] + x_bin_edges[x_i+1])/2
				mid_point_y = (y_bin_edges[x_i] + y_bin_edges[x_i+1])/2
			    
				x_bin.append(mid_point_x)
				y_bin.append(mid_point_y)
			#----------------------------------------------------------------

			xx, yy = np.meshgrid(x_bin, y_bin)

			x, y = xx.ravel(), yy.ravel()
			z = 0

			bin_width_x = x_bin_edges[1] - x_bin_edges[0]
			bin_width_y = y_bin_edges[1] - y_bin_edges[0]

			dx = bin_width_x
			dy = bin_width_y
			dz = counts.T.ravel()

			#plotting 3D bars
			#more style: https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html
			style.use('seaborn-v0_8-colorblind')
			fig = plt.figure()
			ax = fig.add_subplot(111, projection='3d')
			#colors: https://matplotlib.org/stable/gallery/color/named_colors.html
			ax.bar3d(x, y, z, dx, dy, dz, color='cornflowerblue' )
			ax.view_init(elev=24., azim=-33)

			title = 'evID = ' + str(ID) + ' - 10MeV ' 
			x_title = ' cos(α) '
			y_title = ' time residual '
			z_title = 'Counts'

			plt.title(title)
			ax.set_xlabel(x_title)
			ax.set_ylabel(y_title)
			ax.set_zlabel(z_title)

			plt.savefig('figs/All hit Type/' + x_title + y_title + ' evID=' + str(ID) +'.pdf', format='pdf')
			plt.show()

			if cut == True: 

				#Data to be taken
				time_residual_cut = np.array([]) 
				time_residual_hit1_cut = np.array([]) 
				cos_angle_cut = np.array([])                                  

				#cuts
				for i in np.where((np.array(time_residual) > inf_cut) & (np.array(time_residual) < up_cut))[0]:
					time_residual_ut = np.append(time_residual_cut, time_residual[i])
					cos_angle_cut = np.append(cos_angle_cut, cos_angle[i])

				counts, x_bin_edges, y_bin_edges = np.histogram2d(x = cos_angle_cut, y = time_residual_cut, bins = [bins, bins])

				#Obtain coordinates of bins ---------------------------------
				x_bin = []
				y_bin = []

				Nx = len(x_bin_edges)
				Ny = len(y_bin_edges)

				for (x_i, y_i) in zip(range(Nx-1), range(Ny-1)):
					mid_point_x = (x_bin_edges[x_i] + x_bin_edges[x_i+1])/2
					mid_point_y = (y_bin_edges[x_i] + y_bin_edges[x_i+1])/2
				    
					x_bin.append(mid_point_x)
					y_bin.append(mid_point_y)
				#----------------------------------------------------------------

				xx, yy = np.meshgrid(x_bin, y_bin)

				x, y = xx.ravel(), yy.ravel()
				z = 0

				bin_width_x = x_bin_edges[1] - x_bin_edges[0]
				bin_width_y = y_bin_edges[1] - y_bin_edges[0]

				dx = bin_width_x
				dy = bin_width_y
				dz = counts.T.ravel()

				#plotting 3D bars
				#more style: https://matplotlib.org/stable/gallery/style_sheets/style_sheets_reference.html
				style.use('seaborn-v0_8-colorblind')
				fig = plt.figure()
				ax = fig.add_subplot(111, projection='3d')
				#colors: https://matplotlib.org/stable/gallery/color/named_colors.html
				ax.bar3d(x, y, z, dx, dy, dz, color='cornflowerblue' )
				ax.view_init(elev=24., azim=-33)

				title = 'Time Res. Cut - evID = ' + str(ID) + ' - 10MeV ' 
				x_title = ' cos(α) '
				y_title = ' time residual '
				z_title = 'Counts'

				plt.title(title)
				ax.set_xlabel(x_title)
				ax.set_ylabel(y_title)
				ax.set_zlabel(z_title)

				plt.savefig('figs/All hit Type/' + x_title + y_title + str(ID) +'.pdf', format='pdf')
				plt.show()











		





