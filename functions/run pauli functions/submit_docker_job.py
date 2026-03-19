"""
Python script designed to submit partition jobs on the grid. 

This scrip is speciallized on reading MC ratds files
by parts through grid jobs submissions. The script should
entre to the RAT container and develop al the read files 
code within it.
"""

import os
import glob
import re
import time
import sys

# --- Directory Config ---
read_dir = '/share/neutrino/snoplus/MonteCarlo/FullFill_bismsb_801/ScintFit_2p2ppo_2p2bismsB8_Solar_NueRun/ratds/'  # Directory where the RATDS ROOT files are
files_txt_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/bisMSB/B8_solar_Nue/file_name_list/ratds/*.txt'  # Directory of the sublist files to be readen and submited for each job.                                                        # Pattern name of the input files
save_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/bisMSB/B8_solar_Nue/ratDS_output/'
data_type = '8B_Nue_MC'

os.makedirs(f'logs_{data_type}', exist_ok=True)

# Run Container Commands
CONTAINER_DIR = "/lstore/sno/joankl/udocker/udocker/"               # Directory of Container executable 
volume = ("--volume=/lstore/sno/joankl/RAT/rat:/rat "               # Data to be imported within container
		  "--volume=/lstore/sno/joankl:/lstore/sno/joankl "
		  "--volume=/share:/share")              

UDOCKER_EXEC = f'./udocker run --entrypoint="" --bindhome {volume} rat_final /bin/bash' 
UDOCKER_ENV = ('export UDOCKER_DIR="/lstore/sno/joankl/udocker_storage" &&'
				'export UDOCKER_TMPDIR="/lstore/sno/joankl/udocker_tmp" &&') #Define directories for udocker to find the container within its network
#MY_LIBS = "/lstore/sno/joankl/my_pylibs/" # Complementar python libraries
SCRIPT_DIR = os.getcwd() # Return the current directory path, where the current scripts are.

ftxt_list = glob.glob(files_txt_dir)
# Loop on the ftxt_list. Submit  a job for each file txt.
print(f'file txt list: {ftxt_list}')
for i_dx, filetxt_i in enumerate([ftxt_list[0]]):

	# ------- Construct the output filenames -------
	script_name = f"run_job_{data_type}_{i_dx}.sh"  # Script output name

	fout_name = f'Analysis_{i_dx}.root'
	fout_dir = save_dir + fout_name

	# ------- Define the code to be executed within the container on each job -------

	py_one_liner = ("import sys; "
		"from read_RATROOT_MC_flist import extract_data; "
		"extract_data(sys.argv[1], sys.argv[2], sys.argv[3])")

	container_command = ("source scripts/setup-env.sh && "
		"source /rat/env.sh && "
		f"cd {SCRIPT_DIR} && "
		f"python3 -c '{py_one_liner}' '{read_dir}' '{filetxt_i}' '{fout_dir}'")

	# Script Content
	script_content = f"""#!/bin/bash
#SBATCH --job-name={data_type}_{i_dx}
#SBATCH --output=logs_{data_type}/job_{i_dx}.out
#SBATCH --error=logs_{data_type}/job_{i_dx}.err
#SBATCH --partition=lipq

{UDOCKER_ENV}

echo "Running on host: $(hostname)"
echo "Starting Container..."

# Run udocker
# Within the UDOCKER_EXEC is defined the directories to be readen for the analysis
cd {CONTAINER_DIR}
{UDOCKER_EXEC} -c "{container_command}"

echo "Job finished"
"""

	with open(script_name, "w") as f:
	    f.write(script_content)

	# Enviar a la cola
	print(f"Enviando job {i_dx}...")
	os.system(f"sbatch {script_name}")

	# Borrar el script generado para no llenar la carpeta de basura (opcional)
	# os.remove(script_name) 

	time.sleep(1) # Pequeña pausa para no saturar al scheduler


