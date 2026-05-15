'''
Python script dedicated to build histograsm and save
the corresponding plots from numpy array data. Until now
it only calculate the cos(alpha) PDF distribution given various
cuts on time residual.

Creation: 28/05/2026

Edit:	
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

bins = 80
E_inf_cut = 5
#E_cut_list = [5, 6, 8, 10]
#R_cut_list = [5500, 4500, 3500]

t_res_cut_list = [(-5,5), (-5,4), (-5,3), (-5,2), (-5,1)]

# 1 - Prepare the structure to save the counts of the histograms.
#Structure is  the dictionary with entries dict[t_res_cut_list]
bin_edges = np.linspace(-1.0, 1.0, bins + 1)
hist_data = {tres: np.zeros(bins) for tres in t_res_cut_list}

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

		cos_chunk = np.load(read_dir_i + f'cos_alpha_{idx}.npy').astype(np.float32)
		en_chunk = np.load(read_dir_i + f'energy_corr_{idx}.npy').astype(np.float32)
		posr_chunk = np.load(read_dir_i + f'posr_{idx}.npy').astype(np.float32)
		tres_chunk = np.load(read_dir_i + f'hit_residual_{idx}.npy').astype(np.float32)

		# Energy Cut
		E_mask = (en_chunk >= E_inf_cut)

		en_chunk = en_chunk[E_mask]
		#posr_chunk = posr_chunk[t_res_mask]
		cos_chunk = cos_chunk[E_mask]
		tres_chunk = tres_chunk[E_mask]

		# Apply Time residual Cuts in chunck cos_alpha

		for t_res_cut_i in t_res_cut_list:
			tres_mask = (tres_chunk >= t_res_cut_i[0]) & (tres_chunk <= t_res_cut_i[1]) 

			# Compute the histogram
			counts, _ = np.histogram(cos_chunk[tres_mask], bins=bin_edges)

			# Save the counts
			hist_data[t_res_cut_i] += counts

			#print(hist_data)

print('Data fully processed. Generating Plots...')

# =========== Plot Construction ===========

font_style_title = {'family':'serif', 'weight': 'normal','color':'black','size':13}
font_style_axis= {'family':'serif', 'weight': 'normal','color':'black','size':12}
font_prop = font_manager.FontProperties(family=font_style_axis['family'], weight=font_style_axis['weight'], size=font_style_axis['size'])

fig, ax = plt.subplots(figsize = (8, 6))

for t_res_cut_i in t_res_cut_list:
    
	counts = hist_data[t_res_cut_i]

	plt_label = f'time res.: [{t_res_cut_i[0]},{t_res_cut_i[1]}] (ns)'

	plt.hist(bin_edges[:-1], bins=bin_edges, weights=counts, density=True, histtype='step', linewidth=1.5, label = plt_label)

plt.yscale('log')
plt.xlabel(r'cos($\alpha$)', fontdict = font_style_axis)
plt.ylabel(r'Prob. Density', fontdict = font_style_axis)
plt.xlim(-1, 1)

# --- Markers ---
ax.xaxis.set_minor_locator(MultipleLocator(0.2))
ax.xaxis.set_major_locator(MultipleLocator(0.5))

#ax.yaxis.set_minor_locator(MultipleLocator(0.01))
#ax.yaxis.set_major_locator(MultipleLocator(0.05))

ax.tick_params(which='minor', top=True, bottom=True, left=True, right=True)
ax.tick_params(which='major', top=True, bottom=True, left=True, right=True)

plt.legend(loc = 'best', frameon = True, edgecolor='black', prop = font_prop)


plt.title(rf'Directionality Distribution - BisMSB $^8$B-$\nu_e$ MC - E $\geq$ {E_inf_cut} (MeV) & R $\leq$ 5.5 (m)', fontdict = font_style_title)

save_path = save_dir + f'cos_alpha_dir_tres_cuts.png'
plt.savefig(save_path, dpi=300, bbox_inches='tight')

# Limpiamos la figura para que no se superpongan en la RAM
plt.close(fig)
print(f'Saved: {save_path}')

print('Done :)')