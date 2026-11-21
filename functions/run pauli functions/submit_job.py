import os
import sys
import re
import time

# --- Directory Config ---
read_dir = '/share/neutrino/snoplus/Data/FullFill_2p2/rat_801/PPO/'  # Directorio con los datos reales
sublist_dir = '/lstore/sno/joankl/solar_analysis/real_data/2p2ppo/file_name_list/ntuple_for_veto/' # Directorio con las listas de archivos
save_dir = '/lstore/sno/joankl/solar_analysis/real_data/2p2ppo/ntuple_output/bkg_candidates/' # Directorio de salida

# Container config
CONTAINER_SIF = "/lstore/sno/joankl/RAT/containers/rat_8.3.1.sif" # Ruta de la imagen Apptainer
APPTAINER_EXEC = "apptainer" 
VENV_ORIGIN = "/lstore/sno/joankl/.venv_py39" # Ruta de tu entorno virtual original
SCRIPT_DIR = os.getcwd() # Directorio actual donde está el script de análisis

data_type = 'HS_atm_PPO'  
os.makedirs(f'logs_{data_type}', exist_ok=True)

def orden_natural(archivo):
        """Función para ordenar archivos naturalmente por número de run/subrun."""
        return [int(texto) if texto.isdigit() else texto.lower() for texto in re.split('(\d+)', archivo)]

if __name__ == "__main__":

    # Buscar todos los sublist_i.txt ordenados
    sublist_files = [f for f in os.listdir(sublist_dir) if f.startswith("sublist_") and f.endswith(".txt")]
    sublist_files.sort(key=orden_natural)
    print(f"Archivos encontrados: {sublist_files}")

    for i_dx, sublist in enumerate([sublist_files[0]]):
        sublist_path = os.path.join(sublist_dir, sublist)

        # Extraer el número i del sublist
        match = re.search(r"sublist_(\d+)\.txt", sublist)
        sublist_index = match.group(1) if match else "unknown"

        # Crear subdirectorio de salida para cada job
        sub_save_dir = save_dir + f"output_{sublist_index}/"
        os.makedirs(sub_save_dir, exist_ok=True)
        
        script_name = f"run_job_{data_type}_{i_dx}.sh"

        # Comando de Python que se ejecutará al final
        py_cmd = (
            "from veto_and_coincidence_analysis import veto_and_coincidence_analysis; "
            f"veto_and_coincidence_analysis('{read_dir}', '{sublist_path}', '{sub_save_dir}')"
        )

        # Contenido del script SBATCH
        script_content = f"""#!/bin/bash
#SBATCH --job-name={data_type}_{sublist_index}
#SBATCH --output=logs_{data_type}/job_{sublist_index}.out
#SBATCH --error=logs_{data_type}/job_{sublist_index}.err
#SBATCH --partition=lipq
#SBATCH --mem=6G

echo "Running on host: $(hostname)" 

echo "Initializing container ..."

{APPTAINER_EXEC} exec \\
    -B /share/neutrino/snoplus/ \\
    -B /lstore/sno/joankl/ \\
    -B {SCRIPT_DIR} \\
    {CONTAINER_SIF} \\
    bash -s << 'EOF'

echo "Setting up Geant4, ROOT, and RAT environment..."
source /usr/local/bin/geant4.sh
cd /root-bin
source bin/thisroot.sh
cd - > /dev/null
source /rat/env.sh

echo "Copy and Source the virtual environment from /tmp..."
cp -r {VENV_ORIGIN} /tmp/
source /tmp/.venv_py39/bin/activate

echo "Verification - Active Python version and path:"
python --version
which python

echo "Starting Analysis Script ... "
cd {SCRIPT_DIR}

python -c "{py_cmd}"

echo "Job finished successfully inside container"
EOF

echo "Job finished on host"
"""

        with open(script_name, "w") as f:
            f.write(script_content)
        
        # Enviar a la cola de Slurm
        print(f"Enviando job {sublist_index}...")
        os.system(f"sbatch {script_name}")
        
        time.sleep(0.5)

