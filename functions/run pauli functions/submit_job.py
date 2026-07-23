import os
import sys
import re
import time

# --- Directory Config ---
read_dir = '/share/neutrino/snoplus/Data/FullFill_2p2/rat_801/bisMSB/Analysis20_bMR/'  #dir where the real data is
sublist_dir = '/lstore/sno/joankl/solar_analysis/real_data/bisMSB/Analysis20_bMR/file_name_list/ntuple_for_veto/' #dir where the lists with the file names are
save_dir = '/lstore/sno/joankl/solar_analysis/real_data/bisMSB/Analysis20_bMR/results/bkg_candidates/' #dir where we want the output_i/ be generated and save the results

# Container config. for  Workers
CONTAINER_SIF = "/lstore/sno/joankl/RAT/containers/rat_8.1.0.sif" # Container Directory
APPTAINER_EXEC = "/cvmfs/oasis.opensciencegrid.org/mis/apptainer/bin/apptainer" # Look apptainer in the wn node
MY_LIBS = "/lstore/sno/joankl/my_pylibs/" # Complementar python libraries
SCRIPT_DIR = os.getcwd() # Return the current directory path, where the script and rat_read_ntuples.py are.

data_type = 'veto_BisMSB_15bMR'  # Definition of the type of data being processed. Useful to distinguish the job types 

os.makedirs(f'logs_{data_type}', exist_ok=True)

def orden_natural(archivo):
        """Función para ordenar archivos naturalmente por número de run/subrun."""
        return [int(texto) if texto.isdigit() else texto.lower() for texto in re.split('(\d+)', archivo)]

if __name__ == "__main__":

    # Buscar todos los sublist_i.txt
    sublist_files = [f for f in os.listdir(sublist_dir) if f.startswith("sublist_") and f.endswith(".txt")]
    sublist_files.sort(key=orden_natural)  # Orden natural por número i
    print(sublist_files)

    for i_dx, sublist in enumerate(sublist_files):
        sublist_path = os.path.join(sublist_dir, sublist)

        # Extraer el número i
        match = re.search(r"sublist_(\d+)\.txt", sublist)
        sublist_index = match.group(1) if match else "unknown"

        # Crear subdirectorio para guardar outputs: /save_dir/output_i/
        sub_save_dir = save_dir + f"output_{sublist_index}/"
        os.makedirs(sub_save_dir, exist_ok=True)
        
        #Run sbatch scrip Bellow --------------------------------------

        script_name = f"run_job_{data_type}_{i_dx}.sh"

        # Call the Analysis function and pass arguments through sys
        py_one_liner = ("import sys; "
            "from veto_and_coincidence_analysis import veto_and_coincidence_analysis; "
            "veto_and_coincidence_analysis(sys.argv[1], sys.argv[2], sys.argv[3])")

        command_inside_container = (
            "source /usr/local/bin/geant4.sh && "
            "source /root-bin/bin/thisroot.sh && "
            "source /rat/env.sh && "
            f"export PYTHONPATH=$PYTHONPATH:{MY_LIBS}:{SCRIPT_DIR} && "
            f"python3 -c \"{py_one_liner}\" \"{read_dir}\" \"{sublist_path}\" \"{sub_save_dir}\""
        )


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



