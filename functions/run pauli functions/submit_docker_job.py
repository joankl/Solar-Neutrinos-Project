"""
Python script designed to submit partition jobs 
on the grid. 
This scrip is speciallized on reading individual RATDS files with 
the runeventlis.mac and extract the events of interest, saving the 
corresponding RATDS ROOT file. For each ROOT file, submit a grid job.
RAT is executed within a udocker container.
"""
import os
import glob
import re
import time
import sys

# --- Directory Config ---
read_dir = '/share/neutrino/snoplus/Data/FullFill_2p2/rat_801/bisMSB/Analysis15_bMR/ratds/'  # Directory where the RATDS ROOT files are
fin_patter = 'Analysis15_bMR*.root'                                                          # Pattern name of the input files
f_list = glob.glob(read_dir + fin_patter)  # List of full dir + name of the input files
save_dir = '/share/neutrino/snoplus/Data/FullFill_2p2/rat_801/bisMSB/Analysis15_bMR/ratds/run_evlist/'

# Run Container Commands
CONTAINER_DIR = "/lstore/sno/joankl/udocker/udocker/"               # Directory of Container executable 
volume = ("--volume=/lstore/sno/joankl/RAT/rat:/rat "               # Data to be imported within container
		  "--volume=/lstore/sno/joankl:/lstore/sno/joankl "
		  "--volume=/share:/share")              

UDOCKER_EXEC = f'./udocker run --entrypoint="" --hostenv --bindhome {volume} rat_final /bin/bash' 
UDOCKER_ENV = ('export UDOCKER_DIR="/lstore/sno/joankl/udocker_storage" &&'
				'export UDOCKER_TMPDIR="/lstore/sno/joankl/udocker_tmp" &&') #Define directories for udocker to find the container within its network
#MY_LIBS = "/lstore/sno/joankl/my_pylibs/" # Complementar python libraries
SCRIPT_DIR = os.getcwd() # Return the current directory path, where the current scripts are.

data_type = '15_bMR'
os.makedirs(f'logs_{data_type}', exist_ok=True)

# Loop on the file list nad submit jobs
for i_dx, file_i in enumerate([f_list[0]]):

	# ------- Cosntruct the output filenames -------
	script_name = f"run_job_{data_type}_{i_dx}.sh"  # Script output name

	fname = os.path.basename(file_i)
	fname = os.path.splitext(fname)[0]  # Selected filename without format extention. Exm: ('fname', '.root')
	output_i = save_dir + fname + f'_runevelist_{i_dx}.root' # save_dir + output_file name

	# ------- Define the code to be executed within the container on each job -------

	py_one_liner = ("import sys; "
		"from run_eventlist import execute_mac; "
		"execute_mac(sys.argv[1], sys.argv[2])")

	container_command = ("source /home/script/setup-env.sh && "
		"source /rat/env.sh && "
		f"cd {SCRIPT_DIR} && "
		f"python3 -c \"{py_one_liner}\" \"{file_i}\" \"{output_i}\"")

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


