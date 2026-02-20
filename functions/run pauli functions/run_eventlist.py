"""
Python Script designed to automatically run the runeventlist.mac through
python console control. This script will employ a function to read a input file
to be summarized by choosing events by runID and GTID using the runeventlist.mac
rat script.
This script need RAT to be executed.
"""

import subprocess

def execute_mac(fin_dir, fout_dir):

	"""
	Function to execute the mac file.
	Parameters:
	- fin_dir: Input directory of the file.root ratds structure file to be resumed
	- fout_dir: Save directory + filename.root
	"""

	command_line = f"rat -i {fin_dir} -o {fout_dir} runeventlist.mac"
	subprocess.run(command_line, shell=True, check=True)



if __name__ == "__main__":

	fin_dir = '/share/neutrino/snoplus/Data/FullFill_2p2/rat_801/bisMSB/Analysis20_bMR/ratds/Analysis20_bMR_r0000354099_s014_p005.root'
	fout_dir = '/lstore/sno/joankl/solar_analysis/real_data/bisMSB/Analysis20_bMR/test/out.root'

	execute_mac(fin_dir, fout_dir)