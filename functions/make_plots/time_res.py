'''
Plot maker of the Time residual distribution.
Be careful chosing the t_res_min and t_res_max to construct the binedges
of the time residual distribution in the sense that we must respect the 
cuts implements in time residual when constructing the numpy arrays.

Creation: 05/05/2026
'''

import numpy as np

import glob
import os

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, ScalarFormatter
from matplotlib import font_manager


print('Reading Data ...')

read_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/bisMSB/B8_solar_Nue/ratDS_output/np_files/'
read_dir_list = glob.glob(read_dir)
save_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/bisMSB/B8_solar_Nue/ratDS_output/plots/'

os.makedirs(save_dir, exist_ok=True)

bins = 150
E_cut_list = [5]
R_cut_list = [5500]

t_res_min_cut = -50
t_res_max_cut = 200

# 1 - Prepare the structure to save the counts of the histograms.
bin_edges = np.linspace(t_res_min_cut, t_res_max_cut, bins + 1)
hist_data = {E: {R: np.zeros(bins) for R in R_cut_list} for E in E_cut_list}

# 2 - Available files and correct Correspondence of files
base_files = glob.glob(read_dir + 'cos_alpha_*.npy')
indices = [os.path.basename(f).replace('cos_alpha_', '').replace('.npy', '') for f in base_files] #--> index of files cos_alpha_*.npy to be used when loading


print(f'Files to be readen: \n {base_files}')
print(f'Total de chunks a procesar: {len(indices)}')

# 3 - Loading info by chuncks
print('Reading and Histogramming Data in Chunks...')

for read_dir_i in read_dir_list:
	for idx in indices:
		# Load Chuncks
		print(f'Loading Chunk {idx}')
		en_chunk = np.load(read_dir_i + f'energy_corr_{idx}.npy').astype(np.float32)
		posr_chunk = np.load(read_dir_i + f'posr_{idx}.npy').astype(np.float32)
		tres_chunk = np.load(read_dir_i + f'hit_residual_{idx}.npy').astype(np.float32)

		# Apply Energy and Radial Cuts in chunck Data
		for Ecut_i in E_cut_list:
			en_mask = (en_chunk >= Ecut_i)

			for Rcut_i in R_cut_list:
				mask = en_mask & (posr_chunk <= Rcut_i)

				# Compute the histogram
				counts, hist_bin_edges = np.histogram(tres_chunk[mask], bins=bin_edges)
				# Save the counts
				hist_data[Ecut_i][Rcut_i] += counts

print('Data fully processed. Generating Plots...')

print(f'time residual limits defined as: [{t_res_min_cut},{t_res_max_cut}] ns')
print(f'time residual limits found in histogram as: [{min(hist_bin_edges):.0f},{max(hist_bin_edges):.0f}] ns')

# =========== Plot Construction ===========

font_style_title = {'family':'serif', 'weight': 'normal','color':'black','size':13}
font_style_axis= {'family':'serif', 'weight': 'normal','color':'black','size':12}
font_prop = font_manager.FontProperties(family=font_style_axis['family'], weight=font_style_axis['weight'], size=font_style_axis['size'])

for Ecut_i in E_cut_list:
    
	fig, axes = plt.subplots(1, len(R_cut_list), figsize=(8, 6))
	axes = np.atleast_1d(axes)

	for i_dx, Rcut_i in enumerate(R_cut_list):
		ax = axes[i_dx]

		counts = hist_data[Ecut_i][Rcut_i]

		ax.hist(bin_edges[:-1], bins=bin_edges, weights=counts, density=True, histtype='step', color='blue', linewidth=1.5)

		ax.set_yscale('log')
		ax.set_xlabel(r'$t_{res}$ (ns)', fontdict = font_style_axis)
		ax.set_ylabel(r'Prob. Density', fontdict = font_style_axis)
		ax.set_xlim(t_res_min_cut, t_res_max_cut)

		# --- Markers ---
		ax.xaxis.set_minor_locator(MultipleLocator(0.5))
		ax.xaxis.set_major_locator(MultipleLocator(1))

		#ax.yaxis.set_minor_locator(MultipleLocator(0.01))
		#ax.yaxis.set_major_locator(MultipleLocator(0.05))

		ax.tick_params(which='minor', top=True, bottom=True, left=True, right=True)
		ax.tick_params(which='major', top=True, bottom=True, left=True, right=True)

		ax.set_title(rf'Time Residual Distribution - BisMSB $^8$B-$\nu_e$ MC - E $\geq$ {Ecut_i} (MeV) & R $\leq$ {Rcut_i*10**-3:.1f} (m)', fontdict = font_style_title)

	# Titulo General
	#plt.suptitle(rf'Time Residual Distribution - BisMSB $^8$B-$\nu_e$ Candidates - $t_{{res}}$: [{t_res_min_cut:.0f}, {t_res_max_cut:.0f}] (ns)', fontsize=13, y=1.02)

	#save_path = save_dir + f'time_res_E_{Ecut_i}_MeV_All_Rcuts.png'
	save_path = save_dir + f'time_res_E_{Ecut_i}_MeV_R_{Rcut_i}_mm.png'
	plt.savefig(save_path, dpi=300, bbox_inches='tight')

	# Limpiamos la figura para que no se superpongan en la RAM
	plt.close(fig)
	print(f'Saved: {save_path}')

print('Done :)')