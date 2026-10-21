"""
Python script to submit SLURM jobs using Apptainer with rat_8.3.1.sif
to process MC sublists directly into compressed .npz format.
"""

import os
import glob
import re
import time

# --- Configuración de Directorios ---
read_dir = '/share/neutrino/snoplus/MonteCarlo/FullFill_2p2_801/ScintFit_2p2B8_Solar_NueRun/ratds/'
files_txt_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/2p2_ppo/solar_8BNue/file_name_list/ratds/*.txt'
save_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/2p2_ppo/solar_8BNue/ratds_output/npz_files/'
data_type = '8B_Nue_MC_PPO'

CONTAINER_SIF = "/lstore/sno/joankl/RAT/containers/rat_8.3.1.sif"
APPTAINER_EXEC = "apptainer"
SCRIPT_DIR = os.getcwd()

os.makedirs(f'logs_{data_type}', exist_ok=True)
os.makedirs(save_dir, exist_ok=True)

def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]

if __name__ == "__main__":
    ftxt_list = sorted(glob.glob(files_txt_dir), key=natural_sort_key)
    print(f"Encontradas {len(ftxt_list)} sublistas para enviar.")

    for i_dx, filetxt_path in enumerate(ftxt_list):
        # Extraer índice para nombrar salidas
        match = re.search(r"sublist_(\d+)\.txt", os.path.basename(filetxt_path))
        idx_label = match.group(1) if match else str(i_dx)

        script_name = f"run_job_{data_type}_{idx_label}.sh"
        fout_dir = os.path.join(save_dir, f"Analysis_{idx_label}.npz")

        # Comprobar si ya existe para evitar re-ejecución
        if os.path.exists(fout_dir):
            print(f"Saltando sublista {idx_label}: salida ya existe.")
            continue

        py_one_liner = (
            "import sys; "
            "from read_RATROOT_MC_flist import extract_data; "
            "extract_data(sys.argv[1], sys.argv[2], sys.argv[3])"
        )

        command_inside_container = (
            "source /usr/local/bin/geant4.sh && "
            "source /root-bin/bin/thisroot.sh && "
            "source /rat/env.sh && "
            f"cd {SCRIPT_DIR} && "
            f"python3 -c \"{py_one_liner}\" \"{read_dir}\" \"{filetxt_path}\" \"{fout_dir}\""
        )

        script_content = f"""#!/bin/bash
#SBATCH --job-name={data_type}_{idx_label}
#SBATCH --output=logs_{data_type}/job_{idx_label}.out
#SBATCH --error=logs_{data_type}/job_{idx_label}.err
#SBATCH --partition=lipq

echo "Host: $(hostname)"
echo "Iniciando contenedor Apptainer (RAT 8.3.1)..."

{APPTAINER_EXEC} exec \\
    -B /share/neutrino/snoplus/ \\
    -B /lstore \\
    -B {SCRIPT_DIR} \\
    {CONTAINER_SIF} \\
    bash -c '{command_inside_container}'

STATUS=$?
echo "Job terminado con estado: $STATUS"
exit $STATUS
"""

        with open(script_name, "w") as f:
            f.write(script_content)

        print(f"Enviando sublista {idx_label} a SLURM...")
        os.system(f"sbatch {script_name}")
        time.sleep(0.5)