'''
Python script to submit the analysis code to study electron simulations. 
It will submit as many jobs as simulation files were generated.
'''

import glob
import re
import os
import time
import subprocess

# --- Directories settings ---
input_dir = "/lstore/sno/joankl/solar_analysis/mc_data/2p2_ppo/electrons/output_root/"
output_dir = "/lstore/sno/joankl/solar_analysis/mc_data/2p2_ppo/electrons/results_npz/"
scripts_dir = "/lstore/sno/joankl/solar_analysis/mc_data/2p2_ppo/electrons/scripts/"
logs_dir = "/lstore/sno/joankl/solar_analysis/mc_data/2p2_ppo/electrons/logs/"

analysis_script = "/lstore/sno/joankl/solar_analysis/mc_data/2p2_ppo/electrons/Analysis_electron_simulations.py"
CONTAINER_SIF = "/lstore/sno/joankl/RAT/containers/rat_8.3.1.sif"
MY_LIBS = "/lstore/sno/joankl/my_pylibs/" # Complementar python libraries
SCRIPT_DIR = os.getcwd() # Return the current directory path, where the script and rat_read_ntuples.py are.

# Si usas CVMFS o enlaces locales de almacenamiento, bindéalos aquí si Apptainer no lo hace por defecto
apptainer_binds = "-B /lstore:/lstore"

def setup_environment():
    for directory in [output_dir, scripts_dir, logs_dir]:
        os.makedirs(directory, exist_ok=True)

def generate_and_submit():
    setup_environment()
    
    # Buscar todos los archivos .root disponibles
    root_files = sorted(glob.glob(os.path.join(input_dir, "*.root")))
    
    if not root_files:
        print(f"No se encontraron archivos .root en: {input_dir}")
        return

    print(f"Se encontraron {len(root_files)} archivos ROOT para analizar.")

    for filepath in root_files:
        base_name = os.path.basename(filepath)
        file_id = os.path.splitext(base_name)[0]
        
        save_path = os.path.join(output_dir, f"{file_id}.npz")
  
        slurm_script_path = os.path.join(scripts_dir, f"ana_{file_id}.sh")

        # Call the Analysis function and pass arguments through sys
        py_one_liner = ("import sys; "
        "from Analysis_electron_simulations import extract_cherenkov_scint_data; "
        "extract_cherenkov_scint_data(sys.argv[1], sys.argv[2])")

        command_inside_container = (
        "source /usr/local/bin/geant4.sh && "
        "source /root-bin/bin/thisroot.sh && "
        "source /rat/env.sh && "
        f"export PYTHONPATH=$PYTHONPATH:{MY_LIBS}:{SCRIPT_DIR} && "
        f"python3 -c \"{py_one_liner}\" \"{filepath}\" \"{save_path}\""
        )
        
        # Estructura del script SBATCH
        slurm_content = f"""#!/bin/bash
#SBATCH --job-name=ana_{file_id}
#SBATCH --output={logs_dir}/ana_{file_id}.out
#SBATCH --error={logs_dir}/ana_{file_id}.err

echo "=== Starting Analysis of {base_name} in $(hostname) ==="

apptainer exec \\
-B /lstore/sno/joankl \\
-B {SCRIPT_DIR} \\
{CONTAINER_SIF} \\
bash -c '{command_inside_container}'

STATUS=$?
echo "=== Finalizing with Code Status $STATUS ==="
exit $STATUS
"""
        with open(slurm_script_path, "w") as f:
            f.write(slurm_content)

        # Enviar a la cola SLURM
        res = subprocess.run(["sbatch", slurm_script_path], capture_output=True, text=True)
        if res.returncode == 0:
            job_id = res.stdout.strip().split()[-1]
            print(f"Enviado {base_name} -> Job ID: {job_id}")
        else:
            print(f"Error al enviar {base_name}: {res.stderr.strip()}")

if __name__ == "__main__":
    generate_and_submit()



