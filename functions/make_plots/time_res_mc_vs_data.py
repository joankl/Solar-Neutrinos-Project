'''
Script to plot time residual distributions for MC and real data
in the same axis.

created on 09/07/2026

Edit Dates:
- 09/07/2026: Script creation. For now it will compare the time
			  residual distributions of data an MC normalized to
			  the same area.
'''

import numpy as np

import glob
import os

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, ScalarFormatter
from matplotlib import font_manager


# ====== Define the directories to read and save the data ======

read_mc_data_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/bisMSB/B8_solar_Nue/ratDS_output/np_files/'
read_real_data_dir = '/lstore/sno/joankl/solar_analysis/real_data/bisMSB/Analysis*/ratDS_output/np_files/'
real_data_dir_list = glob.glob(read_real_data_dir)
save_dir = '/lstore/sno/joankl/solar_analysis/real_data/bisMSB/plots/figures/'

os.makedirs(save_dir, exist_ok=True)

# ====== Data cuts definitions ======

E_cut_list = [5]
R_cut_list = [5500]

t_res_min_cut = -10
t_res_max_cut = 20

# ====== Histogram Definitions ======
bins = 80

ax_x_min_mark = 1
ax_x_maj_mark = 5

font_style_title = {'family':'serif', 'weight': 'normal','color':'black','size':13}
font_style_axis= {'family':'serif', 'weight': 'normal','color':'black','size':12}
font_prop = font_manager.FontProperties(family=font_style_axis['family'], weight=font_style_axis['weight'], size=font_style_axis['size'])

# Prepare the structure to save the counts of the histograms: Dictionaries with energy and radial cuts as keys
bin_edges = np.linspace(t_res_min_cut, t_res_max_cut, bins + 1)
hist_mc = {E: {R: np.zeros(bins) for R in R_cut_list} for E in E_cut_list}
hist_data = {E: {R: np.zeros(bins) for R in R_cut_list} for E in E_cut_list}

print('Reading Files ...')

# 2 - Available MC files and correct ordering of files index 
mc_base_files = glob.glob(read_mc_data_dir + 'cos_alpha_*.npy')
mc_indices = [os.path.basename(f).replace('cos_alpha_', '').replace('.npy', '') for f in mc_base_files] #--> index of files cos_alpha_*.npy to be used when loading

print(f'MC Files to be readen: \n {mc_base_files}')
print(f'Total MC chunks to process: {len(mc_indices)}')

for idx in mc_indices:
	# Load Chuncks
	print(f'Loading MC Chunk {idx}')
	en_chunk = np.load(read_mc_data_dir + f'energy_corr_{idx}.npy').astype(np.float32)
	posr_chunk = np.load(read_mc_data_dir + f'posr_{idx}.npy').astype(np.float32)
	tres_chunk = np.load(read_mc_data_dir + f'hit_residual_{idx}.npy').astype(np.float32)

	# Apply Energy and Radial Cuts in chunck Data
	for Ecut_i in E_cut_list:
		en_mask = (en_chunk >= Ecut_i)

		for Rcut_i in R_cut_list:
			mask = en_mask & (posr_chunk <= Rcut_i)

			# Compute the histogram
			counts, hist_bin_edges = np.histogram(tres_chunk[mask], bins=bin_edges)
			# Save the counts
			hist_mc[Ecut_i][Rcut_i] += counts

# 3 - Available real data files and correct ordering of files index 
data_base_files = glob.glob(read_real_data_dir + 'cos_alpha_*.npy')
data_indices = [os.path.basename(f).replace('cos_alpha_', '').replace('.npy', '') for f in data_base_files] #--> index of files cos_alpha_*.npy to be used when loading

print(f'Real data Files to be readen: \n {data_base_files}')
print(f'Total Real Data chunks to process: {len(data_indices)}')

for read_real_data_dir_i in real_data_dir_list:
	for idx in data_indices:
		# Load Chuncks
		print(f'Loading real data Chunk {idx}')
		en_chunk = np.load(read_real_data_dir_i + f'energy_corr_{idx}.npy').astype(np.float32)
		posr_chunk = np.load(read_real_data_dir_i + f'posr_{idx}.npy').astype(np.float32)
		tres_chunk = np.load(read_real_data_dir_i + f'hit_residual_{idx}.npy').astype(np.float32)

		# Apply Energy and Radial Cuts in chunck Data
		for Ecut_i in E_cut_list:
			en_mask = (en_chunk >= Ecut_i)

			for Rcut_i in R_cut_list:
				mask = en_mask & (posr_chunk <= Rcut_i)

				# Compute the histogram
				counts, hist_bin_edges = np.histogram(tres_chunk[mask], bins=bin_edges)
				# Save the counts
				hist_data[Ecut_i][Rcut_i] += counts

# =========== Plot Construction ===========

print('Data fully processed. Generating Plots...')

# Generating bins center and widths
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
bin_width = bin_edges[1] - bin_edges[0]

for Ecut_i in E_cut_list:
    
	fig, axes = plt.subplots(1, len(R_cut_list), figsize=(8, 6))
	axes = np.atleast_1d(axes)

	for i_dx, Rcut_i in enumerate(R_cut_list):
		ax = axes[i_dx]

		mc_counts = hist_mc[Ecut_i][Rcut_i]
		data_counts = hist_data[Ecut_i][Rcut_i]

		ax.hist(bin_edges[:-1], bins=bin_edges, weights=mc_counts, density=True, color='blue',
				histtype='step', linewidth=1.5, label=r'$^8$B-$\nu_e$ MC')

		data_total = np.sum(data_counts)
		if data_total > 0:
			norm_factor = data_total * bin_width
			data_pdf = data_counts / norm_factor
			data_errors = np.sqrt(data_counts) / norm_factor

			ax.errorbar(bin_centers, data_pdf, yerr=data_errors, fmt='.', color='black', ecolor='black', 
						elinewidth=1.0, capsize=2, markersize=6, zorder=3, label='Real Data (Stat. Unc.)')

		ax.set_yscale('log')
		ax.set_xlabel(r'$t_{res}$ (ns)', fontdict = font_style_axis)
		ax.set_ylabel(r'Prob. Density', fontdict = font_style_axis)
		ax.set_xlim(t_res_min_cut, t_res_max_cut)

		# --- Markers ---
		ax.xaxis.set_minor_locator(MultipleLocator(ax_x_min_mark))
		ax.xaxis.set_major_locator(MultipleLocator(ax_x_maj_mark))

		#ax.yaxis.set_minor_locator(MultipleLocator(0.01))
		#ax.yaxis.set_major_locator(MultipleLocator(0.05))

		ax.tick_params(which='minor', top=True, bottom=True, left=True, right=True)
		ax.tick_params(which='major', top=True, bottom=True, left=True, right=True)

		ax.set_title(rf'Time Residual Distribution - BisMSB Data vs $^8$B-$\nu_e$ MC' + '\n' + rf'E $\geq$ {Ecut_i} (MeV) & R $\leq$ {Rcut_i*10**-3:.1f} (m)', fontdict = font_style_title)

	# Titulo General
	#plt.suptitle(rf'Time Residual Distribution - BisMSB $^8$B-$\nu_e$ Candidates - $t_{{res}}$: [{t_res_min_cut:.0f}, {t_res_max_cut:.0f}] (ns)', fontsize=13, y=1.02)

	#save_path = save_dir + f'time_res_E_{Ecut_i}_MeV_All_Rcuts.png'
	save_path = save_dir + f'time_res_MC_vs_Data_E_{Ecut_i}_MeV_R_{Rcut_i}_mm.png'
	plt.savefig(save_path, dpi=300, bbox_inches='tight')

	# Limpiamos la figura para que no se superpongan en la RAM
	plt.close(fig)
	print(f'Saved: {save_path}')

print('Done :)')