'''
Python script designed to submit electron simulations to the farm.
The electron are generated uniformly within the inner AV and for random
direnctions.
This script will create the macro files, the output ROOT files, and 
the corresponding file.err and file.out of the sumbited Jobs.
The script allows to define different electron enegies, and number of events
per job.

Created on 30/09/2026
'''

import os
import subprocess

# --- Simulation Setting ---
energies_mev = [2.5, 3.0, 5.0]
events_per_job = 1
jobs_per_energy = 1  

base_dir = "/lstore/sno/joankl/solar_analysis/mc_data/2p2_ppo/electrons/" # working directory

RAT_contained_dir = "/lstore/sno/joankl/RAT/containers/rat_8.3.1_dir"

# command to execute the container
container_cmd = f"apptainer exec {RAT_contained_dir} rat"

# Creation of directories where to save the files
def create_directories():
    for folder in ['macros', 'logs', 'output_root', 'scripts']:
        os.makedirs(os.path.join(base_dir, folder), exist_ok=True)

# Function to generate the macros
def generate_macro(energy, job_idx):
    macro_filename = f"{base_dir}/macros/e_{energy}MeV_{job_idx}.mac"
    output_root = f"{base_dir}/output_root/e_{energy}MeV_{job_idx}.root"
    
    # Macro file content
    macro_content = f"""
/rat/physics_list/OmitMuonicProcesses true
/rat/physics_list/OmitHadronicProcesses true

# /rat/physics_list/OmitCerenkov true
# /rat/physics_list/Optical/OmitBoundaryEffects true

/rat/db/set DETECTOR geo_file "geo/snoplusnative.geo"
/rat/db/set GEO[inner_av] material "labppo_2p2_scintillator"


/run/initialize

/rat/proc prune
/rat/procset prune "mc.track" # Exclude track info.
/rat/proc frontend
/rat/proc trigger
/rat/proc eventbuilder
/rat/proc count
/rat/procset update 10
/rat/proc calibratePMT

/rat/proclast outroot
/rat/procset file "{output_root}"


/generator/add combo gun:fill
/generator/vtx/set e- 0. 0. 0. {energy}
/generator/pos/set inner_av
/generator/rate/set 1

/rat/run/start {events_per_job}
exit
"""
    with open(macro_filename, 'w') as f:
        f.write(macro_content)
    return macro_filename

def generate_slurm_script(energy, job_idx, macro_path):
    slurm_filename = f"{base_dir}/scripts/job_e{energy}_{job_idx}.sh"
    
    slurm_content = f"""#!/bin/bash
#SBATCH --job-name=sno_e{energy}_{job_idx}
#SBATCH --output={base_dir}/logs/job_e{energy}_{job_idx}_%j.out
#SBATCH --error={base_dir}/logs/job_e{energy}_{job_idx}_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

echo "Initializing simulation of e- of {energy} MeV (Job {job_idx})"
{container_cmd} {macro_path}
echo "Simulation completed."
"""
    with open(slurm_filename, 'w') as f:
        f.write(slurm_content)
    return slurm_filename

def main():
    create_directories()
    
    for energy in energies_mev:
        for i in range(jobs_per_energy):
            macro_path = generate_macro(energy, i)
            slurm_path = generate_slurm_script(energy, i, macro_path)
            
            # Enviar el trabajo a SLURM
            subprocess.run(["sbatch", slurm_path])
            print(f"Sending Job {i} with {energy} MeV")

if __name__ == "__main__":
    main()