'''
Python Script designed to receive a list of input files.exactly
and return the corresponding file.dat for data download
'''

import glob
import os
import subprocess
import rat


def from_exac_to_dat(f_in, is_data = False, is_mc = False):

	"""

	Function which will execute the cmd command using RAT 
	to extract the file.dat. You should indicate whether
	is real data (is_data = True and is_mc = False) or MC data
	(data = False and mc = True)

	Parameters:
	- f_in: Input file.exactly directory
	- is_data: Indicate if it is real data

	"""

	f_out_name = os.path.basename(f_in)
	f_out_name = os.path.splitext(f_out_name)[0]

	if is_data == is_mc:
		print(f'arguments data = {is_data} and {is_mc} equals')
		print('Specify the data type! Real data or MC data?"')
		return

	if is_data:
		print(' ------- Reading Real Data ------- ')
		cmd_line = f"./processing_list -e {f_in} -f ratds -o {f_out_name}.dat"
		

	else:
		print(' ------- Reading MC Data ------- ')
		cmd_line = f"./production_list -e {f_in} -f ratds -o {f_out_name}.dat"
		

	try:
		subprocess.run(cmd_line, shell=True, check=True)
		print(f'Successfully extracted file {f_out_name}.dat')

	except subprocess.CalledProcessError as e:
		print(f'error trying to execute the cmd_line: {e.returncode}')


if __name__ == "__main__":

	is_data = False
	is_mc = True

	main_fdir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/bisMSB/download_data/exac_files/solars/*.exactly'
	flist = glob.glob(main_fdir)

	for file_i in flist:

		print(f'reading file {file_i}')
		from_exac_to_dat(file_i, is_data, is_mc)





