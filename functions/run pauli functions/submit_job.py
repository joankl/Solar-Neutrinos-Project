"""
Python script designed to submit partition jobs 
on the grid. 
This scrip is speciallized on reading individual RATDS files with 
the runeventlis.mac and extract the events of interest, saving the 
corresponding RATDS ROOT file. For each ROOT file, submit a grid job.
"""
import os
import glob
import re
import time
import sys

# --- Directory Config ---
read_dir = '/share/neutrino/snoplus/Data/FullFill_2p2/rat_801/bisMSB/Analysis20_bMR/ratds/'  # Directory where the RATDS ROOT files are
fin_patter = 'Analysis20_bMR*.root'                                                          # Pattern name of the input files
f_list = glob.glob(read_dir + fin_patter)  # List of full dir + name of the input files
save_dir = '/lstore/sno/joankl/solar_analysis/real_data/bisMSB/Analysis20_bMR/ratDS_selection/ratds_data/'

# Container config. for  Workers of python libraries
CONTAINER_SIF = "/lstore/sno/joankl/RAT/rat_8.1.0.sif" # Container Directory
APPTAINER_EXEC = "/cvmfs/oasis.opensciencegrid.org/mis/apptainer/bin/apptainer" # Look apptainer in the wn node
MY_LIBS = "/lstore/sno/joankl/my_pylibs/" # Complementar python libraries
SCRIPT_DIR = os.getcwd() # Return the current directory path, where the script and rat_read_ntuples.py are.

data_type = '20_bMR'
os.makedirs(f'logs_{data_type}', exist_ok=True)

# Loop on the file list nad submit jobs
for i_dx, file_i in enumerate(f_list):

	# ------- Cosntruct the output filenames -------
	script_name = f"run_job_{data_type}_{i_dx}.sh"  # Script output name


	fname = os.path.basename(file_i)
	fname = os.path.splitext(fname)[0]  # Selected filename without format extention. Exm: ('fname', '.root')
	output_i = save_dir + fname + f'_runevelist_{i_dx}.root' # save_dir + output_file name

	# ------- Define the code to be executed on each job -------

	py_one_liner = ("import sys; "
					"from run_eventlist import execute_mac; "
					"execute_mac(sys.argv[1], sys.argv[2])")

	container_command = ("source /usr/local/bin/geant4.sh && "
						"source /root-bin/bin/thisroot.sh && "
						"source /rat/env.sh && "
						f"export PYTHONPATH=$PYTHONPATH:{MY_LIBS}:{SCRIPT_DIR} && "
						f"python3 -c \"{py_one_liner}\" \"{file_i}\" \"{output_i}\"")

	# Script Content
	script_content = f"""#!/bin/bash
#SBATCH --job-name={data_type}_{i_dx}
#SBATCH --output=logs_{data_type}/job_{i_dx}.out
#SBATCH --error=logs_{data_type}/job_{i_dx}.err
#SBATCH --partition=lipq

echo "Running on host: $(hostname)"
echo "Starting Container..."

# Run apptainer
# We initialize the directories to be readen: /share (date), /lstore (libs y output), y $PWD (scripts)
{APPTAINER_EXEC} exec \\
    -B /share/neutrino/snoplus/ \\
    -B /lstore \\
    -B {SCRIPT_DIR} \\
    {CONTAINER_SIF} \\
    bash -c '{container_command}'
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


