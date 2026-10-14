'''
Python Script designed to submit the analysis code of PMT_charge_analysis
to SLURM Jobs
'''

import glob
import re
import os
import time
import subprocess

# --- Directories settings ---
input_dir = "/lstore/sno/joankl/solar_analysis/mc_data/2p2_ppo/electrons/results_npz/"
output_dir = "/lstore/sno/joankl/solar_analysis/mc_data/2p2_ppo/electrons/charge_analysis/results/"
scripts_dir = "/lstore/sno/joankl/solar_analysis/mc_data/2p2_ppo/electrons/charge_analysis/scripts/"
logs_dir = "/lstore/sno/joankl/solar_analysis/mc_data/2p2_ppo/electrons/charge_analysis/logs/"

analysis_script = "/lstore/sno/joankl/solar_analysis/mc_data/2p2_ppo/electrons/charge_analysis/PMT_charge_analysis.py"

def setup_environment():
    for directory in [output_dir, scripts_dir, logs_dir]:
        os.makedirs(directory, exist_ok=True)

def generate_and_submit():
    setup_environment()

    # Buscar todos los archivos disponibles
    npz_files = sorted(glob.glob(os.path.join(input_dir, "*.npz")))
    
    if not npz_files:
        print(f"No se encontraron archivos en: {input_dir}")
        return

    print(f"Se encontraron {len(npz_files)} archivos para analizar.")

    for filepath in npz_files:
        base_name = os.path.basename(filepath)
        file_id = os.path.splitext(base_name)[0]
        
        save_path = os.path.join(output_dir, f"{file_id}.npz")
  
        slurm_script_path = os.path.join(scripts_dir, f"ana_{file_id}.sh")

        # Call the Analysis function and pass arguments through sys
        py_one_liner = ("import sys; "
        "from PMT_charge_analysis import calculate_integrated_charge_per_event; "
        "calculate_integrated_charge_per_event(npz_file = sys.argv[1], save_output = sys.argv[2])")

        slurm_command = (f"python3 -c \"{py_one_liner}\" \"{filepath}\" \"{save_path}\"")
        
        # Estructura del script SBATCH
        slurm_content = f"""#!/bin/bash
#SBATCH --job-name=ana_{file_id}
#SBATCH --output={logs_dir}/ana_{file_id}.out
#SBATCH --error={logs_dir}/ana_{file_id}.err

echo "=== Starting Analysis of {base_name} in $(hostname) ==="

bash -c '{slurm_command}'

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
