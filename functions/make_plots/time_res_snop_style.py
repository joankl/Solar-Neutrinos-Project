'''
Plot maker of the Time residual distribution.
Formatted following SNO+ Collaboration standards with Cherenkov hit region highlighted.

Created on 01/09/2026
'''

import glob
import os
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np

# ==========================================
# --- SNO+ Font Registration & Global Settings ---
# ==========================================
font_path = 'Times_New_Roman_Normal.ttf'

if os.path.exists(font_path):
  fm.fontManager.addfont(font_path)
else:
  print(f'Advertencia: Archivo de fuente no encontrado en {font_path}')

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.unicode_minus'] = False

# Standard SNO+ tick labels size
plt.rcParams['xtick.labelsize'] = 21
plt.rcParams['ytick.labelsize'] = 21

plt.rcParams['xaxis.labellocation'] = 'right'
plt.rcParams['yaxis.labellocation'] = 'top'
plt.rcParams['xtick.major.size'] = 6
plt.rcParams['xtick.major.width'] = 1.6
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['xtick.minor.width'] = 1.6
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['ytick.major.width'] = 1.6
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['ytick.minor.width'] = 1.6
plt.rcParams['axes.linewidth'] = 1.6
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['figure.figsize'] = (8.5, 6.5)
plt.rcParams['xtick.major.pad'] = '6'
plt.rcParams['ytick.major.pad'] = '6'
plt.rcParams['axes.titlepad'] = 12
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.top'] = True
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True
plt.rcParams['ytick.right'] = True

# Standard SNO+ single MC histogram style
mc_histogram_style = {
    'histtype': 'step',
    'color': 'blue',
    'linewidth': 2.0,
    'alpha': 1.0,
}
# ==========================================

print('Reading Data ...')

read_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/2p2_ppo/solar_8BNue/ratds_output/np_files/'
read_dir_list = glob.glob(read_dir)
save_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/2p2_ppo/solar_8BNue/ratds_output/plots/'

os.makedirs(save_dir, exist_ok=True)

bins = 150
E_cut_list = [5]
R_cut_list = [5500]

t_res_min_cut = -5
t_res_max_cut = 100

# Cherenkov region boundaries
t_res_cher_min = -3
t_res_cher_max = 1

bin_edges = np.linspace(t_res_min_cut, t_res_max_cut, bins + 1)
hist_data = {E: {R: np.zeros(bins) for R in R_cut_list} for E in E_cut_list}

base_files = glob.glob(read_dir + 'cos_alpha_*.npy')
indices = [
    os.path.basename(f).replace('cos_alpha_', '').replace('.npy', '')
    for f in base_files
]

print(f'Files to be readen: \n {base_files}')
print(f'Total de chunks a procesar: {len(indices)}')

print('Reading and Histogramming Data in Chunks...')

for read_dir_i in read_dir_list:
  for idx in indices:
    print(f'Loading Chunk {idx}')
    en_chunk = np.load(read_dir_i + f'energy_corr_{idx}.npy').astype(np.float32)
    posr_chunk = np.load(read_dir_i + f'posr_{idx}.npy').astype(np.float32)
    tres_chunk = np.load(read_dir_i + f'hit_residual_{idx}.npy').astype(
        np.float32
    )

    for Ecut_i in E_cut_list:
      en_mask = en_chunk >= Ecut_i

      for Rcut_i in R_cut_list:
        mask = en_mask & (posr_chunk <= Rcut_i)
        counts, hist_bin_edges = np.histogram(tres_chunk[mask], bins=bin_edges)
        hist_data[Ecut_i][Rcut_i] += counts

print('Data fully processed. Generating Plots...')

# =========== Plot Construction ===========

for Ecut_i in E_cut_list:
  fig, axes = plt.subplots(1, len(R_cut_list), figsize=(8.5, 6.5))
  axes = np.atleast_1d(axes)

  for i_dx, Rcut_i in enumerate(R_cut_list):
    ax = axes[i_dx]
    ax.minorticks_on()

    # Region Highlight
    ax.axvspan(
        t_res_cher_min,
        t_res_cher_max,
        color='gray',
        alpha=0.4,
    )
    ax.axvline(t_res_cher_min, color='gray', linestyle='--', linewidth=1.4)
    ax.axvline(t_res_cher_max, color='gray', linestyle='--', linewidth=1.4)

    counts = hist_data[Ecut_i][Rcut_i]

    # MC Histogram
    ax.hist(
        bin_edges[:-1],
        bins=bin_edges,
        weights=counts,
        density=True,
        **mc_histogram_style,
    )

    # Labels & limits
    ax.set_xlabel(r'$t_{\rm res}$ [ns]', size=23)
    ax.set_ylabel('Probability Density', size=23)
    ax.set_xlim(t_res_min_cut, t_res_max_cut)

    # Locators
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.xaxis.set_major_locator(MultipleLocator(25))

    # Title
    #ax.set_title(
    #    r'Solar $^8$B-$\nu_e$ - 2.2 g/L PPO MC'
    #    + '\n'
    #    + rf'E $\geq$ {Ecut_i} MeV & R $\leq$ {Rcut_i*1e-3:.1f} m',
    #    size=25,
    #    y=1.02)

    # Custom legend with line proxy
    cherenkov_label = (
        'Most probable Cherenkov hits\n'
        rf'$t_{{\rm res}} \in [{t_res_cher_min}, {t_res_cher_max}]$ ns'
    )
    custom_handles = [
        Line2D(
            [0],
            [0],
            color=mc_histogram_style['color'],
            linestyle='-',
            linewidth=mc_histogram_style['linewidth'],
            label=r'$^8$B-$\nu_e$ MC',
        ),
        Patch(
            facecolor='gray',
            edgecolor='none',
            alpha=0.4,
            label=cherenkov_label,
        ),
    ]
    ax.legend(
        handles=custom_handles, loc='center right', frameon=False, fontsize=16
    )

    # SNO+ Preliminary watermark and legend of data type
    ax.text(0.25,0.90,'SNO+ Preliminary',transform=ax.transAxes,color='black',size=22)

    ax.text(0.25, 0.83, "2.2 g/L PPO MC", transform=ax.transAxes, color='black', size=15)

  # Save vector PDF[cite: 1]
  save_path = os.path.join(save_dir, f'time_res_{t_res_min_cut}_{t_res_max_cut}_ns_E_{Ecut_i}_MeV_R_{Rcut_i}_mm_snopl_style.pdf')
  plt.savefig(save_path, format='pdf', bbox_inches='tight')

  plt.close(fig)
  print(f'Saved: {save_path}')

print('Done :)')