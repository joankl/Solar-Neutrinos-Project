"""
Python script designed to read  the RATDS ROOT Files
and extract observables of interest. For now, the script
doesn't perform any analysis/cuts on the data. It will read 
data already analized and save RATDS observables.
"""

import rat
import ROOT
from array import array

fdir_load = ''
fdir_save = ''

def extract_data(read_dir, save_dir):

	"""
	Wrapper function to read and save the PMT info and events info.
	Parameters:

	- read_dir: Directory where the RATDS ROOT files are.
	- save_dir: Directory to save the ouput files
	-
	"""

	print('Extracting RATDS ROOT Data')

	def write_pmt_info(fout,reload=False):

		"""
		Function to save the PMT Information
		"""
		
		du = rat.utility()

		if reload:
			du.LoadDBAndBeginRun()

		pmtinfo = du.GetPMTInfo()
		fout.cd()

		# PMT TTree
		tpmt = ROOT.TTree("pmt","PMT information")

		pmt_id = array('i',[0])
		tpmt.Branch('pmt_id',pmt_id,'pmt_id/I')

		pmt_pos_xyz = array('d',3*[0.0])
		tpmt.Branch('pmt_pos_xyz',pmt_pos_xyz,'pmt_pos_xyz[3]/D')

		pmt_pos_sph = array('d',3*[0.0])
		tpmt.Branch('pmt_pos_sph',pmt_pos_sph,'pmt_pos_sph[3]/D')

		pmt_type = array('i',[0])
		tpmt.Branch('pmt_type',pmt_type,'pmt_type/I')

		for i_pmt in range(0,pmtinfo.GetCount()):
			pos = pmtinfo.GetPosition(i_pmt)
			pmt_pos_xyz[0] = pos.x()
			pmt_pos_xyz[1] = pos.y()
			pmt_pos_xyz[2] = pos.z()
			pmt_pos_sph[0] = pos.Theta()
			pmt_pos_sph[1] = pos.Phi()
			pmt_pos_sph[2] = pos.Mag()
			pmt_id[0] = i_pmt
			pmt_type[0] = pmtinfo.GetType(i_pmt)
			tpmt.Fill()
    
    tpmt.Write()



