'''
Python script dedicated to build histograms and save
the corresponding plots from numpy array data.
Calculates the cos(alpha) PDF distribution given various
cuts on time residual.

Formatted following SNO+ Collaboration standards.

created on 01/09/2026
'''

import glob
import os
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MultipleLocator
import numpy as np

# ==========================================
# --- SNO+ Global Plot Settings ---
# ==========================================
font_path = 'Times_New_Roman_Normal.ttf'

if os.path.exists(font_path):
  fm.fontManager.addfont(font_path)
else:
  # Alternativa: colocar aquí la ruta absoluta si está en otro directorio del cluster
  print(f'Advertencia: Archivo de fuente no encontrado en {font_path}')

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
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

colors = ['#0072B2', '#D55E00', '#009E73', '#D50000']
linestyles = ['-', '--', '-.', ':']
# ==========================================

print('Reading Data ...')

read_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/2p2_ppo/solar_8BNue/ratds_output/np_files/'
read_dir_list = glob.glob(read_dir)
save_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/2p2_ppo/solar_8BNue/ratds_output/plots/'

os.makedirs(save_dir, exist_ok=True)

bins = 80
E_inf_cut = 5
t_res_cut_list = [(-3, 4), (-3, 3), (-3, 2), (-3, 1)]

bin_edges = np.linspace(-1.0, 1.0, bins + 1)
hist_data = {tres: np.zeros(bins) for tres in t_res_cut_list}

base_files = glob.glob(read_dir + 'cos_alpha_*.npy')
indices = [
    os.path.basename(f).replace('cos_alpha_', '').replace('.npy', '')
    for f in base_files
]

for read_dir_i in read_dir_list:
  for idx in [indices[0]]:
    print(f'Loading Chunk {idx}')
    cos_chunk = np.load(read_dir_i + f'cos_alpha_{idx}.npy').astype(np.float32)
    en_chunk = np.load(read_dir_i + f'energy_corr_{idx}.npy').astype(np.float32)
    tres_chunk = np.load(read_dir_i + f'hit_residual_{idx}.npy').astype(
        np.float32
    )

    E_mask = en_chunk >= E_inf_cut
    cos_chunk = cos_chunk[E_mask]
    tres_chunk = tres_chunk[E_mask]

    for t_res_cut_i in t_res_cut_list:
      tres_mask = (tres_chunk >= t_res_cut_i[0]) & (tres_chunk <= t_res_cut_i[1])
      counts, _ = np.histogram(cos_chunk[tres_mask], bins=bin_edges)
      hist_data[t_res_cut_i] += counts

print('Data fully processed. Generating Plots...')

# =========== Plot Construction ===========

fig, ax = plt.subplots(figsize=(8.5, 6.5))
ax.minorticks_on()

for i, t_res_cut_i in enumerate(t_res_cut_list):
  counts = hist_data[t_res_cut_i]
  plt_label = rf'$t_{{\rm res}} \in [{t_res_cut_i[0]}, {t_res_cut_i[1]}]$ ns'

  ax.hist(
      bin_edges[:-1],
      bins=bin_edges,
      weights=counts,
      density=True,
      histtype='step',
      linewidth=2.0,
      color=colors[i % len(colors)],
      linestyle=linestyles[i % len(linestyles)],
      label=plt_label,
  )

ax.set_yscale('log')

# --- FORMATEADOR SIN MATHTEXT PARA EVITAR GLIFOS CORRUPTOS ---
superscript_map = {
    '-': '⁻',
    '0': '⁰',
    '1': '¹',
    '2': '²',
    '3': '³',
    '4': '⁴',
    '5': '⁵',
    '6': '⁶',
    '7': '⁷',
    '8': '⁸',
    '9': '⁹',
}


def clean_log_formatter(val, pos):
  if val <= 0:
    return ''
  exp = int(np.floor(np.log10(val)))
  factor = val / (10**exp)
  exp_str = ''.join(superscript_map.get(c, c) for c in str(exp))

  if np.isclose(factor, 1.0, atol=1e-2):
    return f'10{exp_str}'
  elif np.isclose(factor, np.round(factor), atol=1e-2):
    return f'{int(np.round(factor))}×10{exp_str}'
  else:
    return f'{factor:.1f}×10{exp_str}'


ax.yaxis.set_major_formatter(FuncFormatter(clean_log_formatter))
# -------------------------------------------------------------

ax.set_xlabel(r'$\cos(\alpha)$', size=26)
ax.set_ylabel('Probability Density', size=26)
ax.set_xlim(-1, 1)

ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.xaxis.set_major_locator(MultipleLocator(0.5))

ax.legend(loc='upper left', frameon=False, fontsize=16)  #[cite: 1, 7]
ax.set_title(
    r'Solar $^8$B-$\nu_e$ Directionality - 2.2PPO MC'
    + '\n'
    + rf'E $\geq$ {E_inf_cut} MeV & R $\leq$ 5.5 m',
    size=18,
    y=1.02,
)

ax.text(
    0.05,
    0.08,
    'SNO+ Preliminary',
    transform=ax.transAxes,
    color='black',
    size=22,
)  #[cite: 1, 7]

save_path = os.path.join(save_dir, 'cos_alpha_dir_tres_cuts_snop_style.pdf')
plt.savefig(save_path, format='pdf', bbox_inches='tight')  #[cite: 1, 7]

plt.close(fig)
print(f'Saved: {save_path}')