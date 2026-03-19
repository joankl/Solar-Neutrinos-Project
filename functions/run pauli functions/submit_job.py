import os
import re
import time
import sys

# --- Directory Config ---
read_dir = '/share/neutrino/snoplus/Data/FullFill_2p2/rat_801/bisMSB/Analysis15/'
sublist_dir = '/lstore/sno/joankl/solar_analysis/real_data/bisMSB/Analysis15/file_name_list/'
save_dir = '/lstore/sno/joankl/solar_analysis/real_data/bisMSB/Analysis15/results/output_files/'

# Container config. for  Workers
CONTAINER_SIF = "/lstore/sno/joankl/RAT/rat_8.1.0.sif" # Container Directory
APPTAINER_EXEC = "/cvmfs/oasis.opensciencegrid.org/mis/apptainer/bin/apptainer" # Look apptainer in the wn node
MY_LIBS = "/lstore/sno/joankl/my_pylibs/" # Complementar python libraries
SCRIPT_DIR = os.getcwd() # Return the current directory path, where the script and rat_read_ntuples.py are.

data_type = '15'
os.makedirs(f'logs_{data_type}', exist_ok=True)

def orden_natural(archivo):
    return [int(texto) if texto.isdigit() else texto.lower() for texto in re.split('(\d+)', archivo)]

if __name__ == "__main__":

    # load the name of files to be analyzed and sort them by natural order (run_subrun)
    sublist_files = [f for f in os.listdir(sublist_dir) if f.startswith("sublist_") and f.endswith(".txt")]
    sublist_files.sort(key=orden_natural)
    
    print(f"Found {len(sublist_files)} sublists to analyze.")

    # Loop over the list of files
    for i_dx, sublist in enumerate(sublist_files):
        sublist_path = os.path.join(sublist_dir, sublist)
        
        # Extract index
        match = re.search(r"sublist_(\d+)\.txt", sublist)
        sublist_index = match.group(1) if match else "unknown"
        
        sub_save_dir = os.path.join(save_dir, f"output_{sublist_index}/")
        
        # Temporal Scrf"python3 -c \"{py_one_liner}\" \"{read_dir}\" \"{sublist_path}\" \"{sub_save_dir}\""ipt Name
        script_name = f"run_job_{data_type}_{i_dx}.sh"

        # Call the Analysis function and pass arguments through sys
        py_one_liner = ("import sys; "
            "from rat_read_ntuples import extract_data; "
            "extract_data(sys.argv[1], sys.argv[2], sys.argv[3])")
        
        # --- HERE IS THE MAGIC	 ---
	# We'll write the script that apptainer will execute
	# We use -B to define our disk on the container and we config. the environment through the command bash -c (open a consol)
        
        # We need to source geant4, root, rat, and  my python libs to run rat_readntuples.py within the cointaner
        command_inside_container = (
            "source /usr/local/bin/geant4.sh && "
            "source /root-bin/bin/thisroot.sh && "
            "source /rat/env.sh && "
            f"export PYTHONPATH=$PYTHONPATH:{MY_LIBS}:{SCRIPT_DIR} && "
            f"python3 -c \"{py_one_liner}\" \"{read_dir}\" \"{sublist_path}\" \"{sub_save_dir}\""
        )

        # The line f"python3 -c \"{py_one_liner}\" \"{read_dir}\" \"{sublist_path}\" \"{sub_save_dir}\""
        # make uses of sys, to define the following consol outputs as the inputs of extract_data functions.

        # Content of script SBATCH
        script_content = f"""#!/bin/bash
#SBATCH --job-name={data_type}_{sublist_index}
#SBATCH --output=logs_{data_type}/job_{sublist_index}.out
#SBATCH --error=logs_{data_type}/job_{sublist_index}.err
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
    bash -c '{command_inside_container}'

echo "Job finished"
"""

        with open(script_name, "w") as f:
            f.write(script_content)
        
        # Enviar a la cola
        print(f"Enviando job {sublist_index}...")
        os.system(f"sbatch {script_name}")
        
        # Borrar el script generado para no llenar la carpeta de basura (opcional)
        # os.remove(script_name) 
        
        time.sleep(1) # Pequeña pausa para no saturar al scheduler
