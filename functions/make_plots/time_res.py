'''
Plot maker of the Time residual distribution.
Be careful chosing the t_res_min and t_res_max to construct the binedges
of the time residual distribution

Creation: 05/05/2026
'''

import numpy as np
import glob
import os
import matplotlib.pyplot as plt


print('Reading Data ...')

read_dir = '/lstore/sno/joankl/solar_analysis/real_data/bisMSB/Analysis*/ratDS_output/np_files/'
read_dir_list = glob.glob(read_dir)
save_dir = '/lstore/sno/joankl/solar_analysis/real_data/bisMSB/plots/figures/'

os.makedirs(save_dir, exist_ok=True)

bins = 80
E_cut_list = [5, 6, 8, 10]
R_cut_list = [5500, 4500, 3500]

t_res_min_cut = -150
t_res_max_cut = 250

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

# =========== Plot Construction ===========
print(f'time residual limits defined as: [{t_res_min_cut},{t_res_max_cut}] ns')
print(f'time residual limits found in histogram as: [{min(hist_bin_edges):.0f},{max(hist_bin_edges):.0f}] ns')

for Ecut_i in E_cut_list:
    
	fig, axes = plt.subplots(1, len(R_cut_list), figsize=(24, 5.8))

	for i_dx, Rcut_i in enumerate(R_cut_list):
		ax = axes[i_dx]

		# Recuperamos los conteos que calculamos en el loop de lectura
		counts = hist_data[Ecut_i][Rcut_i]

		ax.hist(bin_edges[:-1], bins=bin_edges, weights=counts, histtype='step', color='blue', linewidth=1.5)

		ax.set_title(rf'E $\geq$ {Ecut_i} (MeV) & R $\leq$ {Rcut_i*10**-3:.1f} (m)')
		ax.set_yscale('log')
		ax.set_xlabel(r'$t_{res}$')
		ax.set_xlim(t_res_min_cut, t_res_max_cut)

	# Titulo General
	plt.suptitle(rf'Time Residual Distribution - BisMSB $^8$B-$\nu_e$ Candidates - $t_{{res}}$: [{t_res_min_cut:.0f}, {t_res_max_cut:.0f}] (ns)', fontsize=13, y=1.02)

	save_path = save_dir + f'time_res_E_{Ecut_i}_MeV_All_Rcuts.png'
	plt.savefig(save_path, dpi=300, bbox_inches='tight')

	# Limpiamos la figura para que no se superpongan en la RAM
	plt.close(fig)
	print(f'Saved: {save_path}')

print('Done :)')