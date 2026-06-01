"""
Python Script designed to automatically run the runeventlist.mac through
python console control. This script will employ a function to read a input file
to be summarized by choosing events by runID and GTID using the runeventlist.mac
rat script.
This script need RAT to be executed.
"""

import subprocess
import glob
import os

def execute_mac(fin_dir, fout_dir):

	"""
	Function to execute the mac file.
	Parameters:
	- fin_dir: Input directory of the file.root ratds structure file to be resumed
	- fout_dir: Save directory + filename.root
	"""

	command_line = f"rat -i {fin_dir} -o {fout_dir} runeventlist.mac"
	subprocess.run(command_line, shell=True, check=True)

# Run to Loop over files
if __name__ == "__main__":

	print('Function activated! DONT SUBMIT JOBS?')
	read_dir = '/share/neutrino/snoplus/Data/FullFill_2p2/rat_801/PPO/ratds/'	# Directory where the RATDS ROOT files are
	fin_patter = 'Analysis20_PPOR*.root'										# Pattern name of the input files
	f_list = glob.glob(read_dir + fin_patter)  									# List of full dir + name of the input files

	save_dir = '/share/neutrino/snoplus/Data/FullFill_2p2/rat_801/PPO/ratds/run_evlist/'

	print(f'files to be readen: {f_list}')

	for i_dx, file_i in enumerate(f_list):
		print(f'On loop number {i_dx} of {len(f_list)}')


		fname = os.path.splitext(os.path.basename(file_i))[0]
		execute_mac(file_i, save_dir + fname + '_resume.root')

	print('Done :-)')
