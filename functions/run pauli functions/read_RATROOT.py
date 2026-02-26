"""
Python script designed to read a RATDS ROOT file
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

	# Create Output file
	fout = ROOT.TFile.Open(save_dir,"RECREATE")
	fout.cd()

	# ====== Define the Branches of the output event data TTree ======
	tree = ROOT.TTree("T", "output_summary")

	energy = array('d', [0.0])
	tree.Branch('energy', energy, 'energy/D')

	energy_corr = array('d', [0.0])
	tree.Branch('energy_corr', energy_corr, 'energy_corr/D')

	evtid = array('i',[0])
	tree.Branch('evtid',evtid,'evtid/I')

	pos_xyz = array('d',3*[0.0])
	tree.Branch('position',pos_xyz, 'position[3]/D')

	sun_dir = array('d',3*[0.0])
	tree.Branch('sun_dir',sun_dir, 'sun_dir[3]/D')

	evt_time_day = array('d', [0.0]) 
	tree.Branch('evt_time_day',evt_time_day, 'evt_time_day/D')

	evt_time_sec = array('d', [0.0])
	tree.Branch('evt_time_sec',evt_time_sec, 'evt_time_sec/D')

	evt_time_nsec = array('d', [0.0])
	tree.Branch('evt_time_nsec',evt_time_nsec, 'evt_time_nsec/D')

	hit_pmtid = array('i',[0])
	tree.Branch('hit_pmtid',hit_pmtid,'hit_pmtid/I')

	hit_type = array('i',[0])
	tree.Branch('hit_type',hit_type,'hit_type/I')

	# ------ Load Utility Tools ------
	util = rat.utility()
	util.LoadDBAndBeginRun()

	# Sun Direction Function
	SunDir = rat.RAT.SunDirection

	# Energy Callibrator Definitions 
	calibrator = util.GetReconCalibrator()
	MATERIAL_NAME = "labppo_2p2_bismsb_2p2_scintillator"
	CORRECTION_VER = 3   # VER = 2 for bisMSB data/MC
	IS_DATA = True   # False if MC
	P3D = ROOT.RAT.DU.Point3D
	av_id = P3D.GetSystemId("av")

	# ====== Reading file root info ======
	reader = ROOT.RAT.DU.DSReader(read_dir)

	# Loop on entries in root file
	for ievent in range(0, reader.GetEntryCount()):
		rDS = reader.GetEntry(ievent)

		# Loop on each triggered event in the entry
		for iev in range(0, rDS.GetEVCount()):

			# Select the Getters
			rEV = rDS.GetEV(iev) # EV Brach
			fResult = rEV.GetFitResult("scintFitter")           # Fitter Branch
			fPosition = fResult.GetVertex(ivtx).GetPosition()
			calibratedPMTs = rEV.GetCalPMTs()                   # PMTs Brach
			rTime = rEV.GetUniversalTime()                      # Event Time Branch

			# ------- Fill Branches -------
			evtid[0] = rEV.GetGTID()  # GTID

			# Reconstructed Position
			pos_xyz[0] = fPosition.x() 
			pos_xyz[1] = fPosition.y()
			pos_xyz[2] = fPosition.z()

			position = P3D(fPosition.x(), fPosition.y(), fPosition.z())

			# Reconstructed Time
			evt_time_day[0] = rTime.GetDays() 
			evt_time_sec[0] = rTime.GetSeconds()
			evt_time_nsec[0] = rTime.GetNanoSeconds()

			# Compute Sun Direction
			sun_dir = SunDir(rTime.GetDays(), rTime.GetSeconds(), rTime.GetNanoSeconds())

			# Reconstructed Energy
			for ivtx in range(0, fResult.GetVertexCount()): 
				fVertex = fResult.GetVertex(ivtx)
				energy[0] = fVertex.GetEnergy()
				print(f'event energy: {energy[0]}')

				energy_corr[0] = calibrator.CalibrateEnergyRTF(IS_DATA, fVertex.GetEnergy(), 
					position, MATERIAL_NAME, CORRECTION_VER)
				print(f'Corrected energy: {energy_corr[0]}')

			# Loop over the hits
			for iPMT in range(0, calibratedPMTs.GetAllCount()):
				pmtCal = calibratedPMTs.GetAllPMT(iPMT)
				pmtid = pmtCal.GetID()

				hit_pmtid[0] = pmtid

				tree.Fill()

	tree.Write()
	write_pmt_info(fout)
	fout.Close()

if __name__ == "__main__":
	read_dir = '/lstore/sno/joankl/solar_analysis/real_data/bisMSB/Analysis20_bMR/ratDS_selection/ratds_data/Analysis20_bMR_r0000358052_s002_p005_runevelist_42.root'
	save_dir = '/lstore/sno/joankl/solar_analysis/real_data/bisMSB/Analysis20_bMR/test/'
	file_out_name = 'x.root'
	extract_data(read_dir, save_dir + file_out_name)




			















