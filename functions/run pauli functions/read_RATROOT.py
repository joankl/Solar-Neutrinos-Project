"""
Python script designed to read a RATDS ROOT file of real data
and extract observables of interest. For now, the script
doesn't perform any analysis/cuts on the data. It will read 
data already analized and save RATDS observables.
This script assumes already processed RATDS! It uses known
reconstructe events, so doesnt employs fit checks or pile-up events!
"""

import rat
import ROOT
from array import array
import glob


def extract_data(read_dir, file_txt_dir, save_dir):

	"""
	Wrapper function to read and save the PMT info and events info.
	Parameters:

	- read_dir: Directory where the RATDS ROOT files are.
	- file_txt_dir: Directory of the txt files which contains the list
	                of the input_files
	- save_dir: Directory to save the ouput files

	"""

	def read_files_txt(file_txt_dir):
		"""
		Function that reads the file_txt_dir and returns 
		the list of the file name to read the RATDS files
		"""
		with open(file_txt_dir, "r", encoding="utf-8") as f:
			file_names = [line.strip() for line in f if line.strip()]
		return file_names

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

	# Save PMT Info.
	write_pmt_info(fout, reload = True)
	fout.Flush()  # Force to write metadata in disk


	fout.cd()

	# ====== Define the Branches of the output event data TTree ======

	print('creating TTree Branches')
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

	hit_residual = array('d',[0.])
	tree.Branch('hit_residual',hit_residual,'hit_residual/D')

	# ------ Load Utility Tools ------
	util = rat.utility()
	util.LoadDBAndBeginRun()

	# Sun Direction Function
	SunDir = rat.RAT.SunDirection

	# Time residual calculator
	#timeResCalc = util.GetTimeResidualCalculator()

	# light Path Calculator
	light_path_cal = util.GetLightPathCalculator()


	# Group Velocity
	group_velocity = util.GetGroupVelocity()

	# pmt Information
	pmtinfo = util.GetPMTInfo()

	# Point3D with PSUP reference
	P3D = ROOT.RAT.DU.Point3D
	psup_id = P3D.GetSystemId("innerPMT") 

	# Energy Callibrator Definitions 
	#print ('Takinkg Energy Calibrator Tool')
	calibrator = util.GetReconCalibrator()
	#MATERIAL_NAME = "labppo_2p2_scintillator"
	MATERIAL_NAME = "labppo_2p2_bismsb_2p2_scintillator"
	CORRECTION_VER = 3   
	IS_DATA = True   # False if MC
	av_id = P3D.GetSystemId("av")

	# ====== Reading file root info ======
	flist = read_files_txt(file_txt_dir)
	assert len(flist) != 0, f"Unexpected number of input files. Got {0}"
	print(f'reading file list {flist}')

	for input_file_i in flist:
		full_dir_flist = os.path.join(read_dir, input_file_i)  #Construct the full directory of the input files
		
		print('Getting in reader ...')
		reader = ROOT.RAT.DU.DSReader(full_dir_flist)
		print(f'reader here: {reader}', flush=True)

		# Loop on entries in root file
		for ievent in range(0, reader.GetEntryCount()):
			print(f'reading root file at entry {ievent}')
			rDS = reader.GetEntry(ievent)

			# Loop on entries in root file
			for iev in range(0, rDS.GetEVCount()):
				print(f'reading triggered event at entry {iev}')

				# Select the Getters
				rEV = rDS.GetEV(iev) # EV Brach
				fResult = rEV.GetFitResult("scintFitter")           # Fitter Branch
				fPosition = fResult.GetVertex(0).GetPosition()      # Use the first fitted position
				calibratedPMTs = rEV.GetCalPMTs()                   # PMTs Brach
				rTime = rEV.GetUniversalTime()                      # Event Time Branch

				# ------- Fill Branches -------
				evtid[0] = rEV.GetGTID()  # GTID

				# Reconstructed Position
				pos_xyz[0] = fPosition.x() 
				pos_xyz[1] = fPosition.y()
				pos_xyz[2] = fPosition.z()

				#position = P3D(fPosition.x(), fPosition.y(), fPosition.z())

				# Reconstructed Time
				evt_time_day[0] = rTime.GetDays() 
				evt_time_sec[0] = rTime.GetSeconds()
				evt_time_nsec[0] = rTime.GetNanoSeconds()

				# Compute Sun Direction
				sun_dir_vector = SunDir(int(rTime.GetDays()), int(rTime.GetSeconds()), int(rTime.GetNanoSeconds()))

				sun_dir[0] = sun_dir_vector.X()
				sun_dir[1] = sun_dir_vector.Y()
				sun_dir[2] = sun_dir_vector.Z()

				# Reconstructed Energy
				fVertex = fResult.GetVertex(0)
				energy[0] = fVertex.GetEnergy()
				print(f'event energy: {energy[0]}')

				# Energy Correction
				positionP3D = P3D(av_id, pos_xyz[0], pos_xyz[1], pos_xyz[2])
				energy_correction = calibrator.CalibrateEnergyRTF(IS_DATA, energy[0], positionP3D, MATERIAL_NAME, CORRECTION_VER)
				energy_corr[0] = energy_correction

				fVertexTime =fVertex.GetTime()
				fit_pos_ev_3d = P3D(psup_id, fPosition.x(), fPosition.y(), fPosition.z())

				
				# Loop over the hits
				calibratedPMTs = rEV.GetCalPMTs()                   # PMTs Brach
				pmtCalStatus = rat.utility().GetPMTCalStatus()
				for iPMT in range(0, calibratedPMTs.GetAllCount()):
					pmtCal = calibratedPMTs.GetAllPMT(iPMT)
					if pmtCalStatus.GetHitStatus(pmtCal) != 0:
						continue

					pmt_point = P3D(psup_id, pmtinfo.GetPosition(pmtCal.GetID()))
					light_path_cal.CalcByPosition(fit_pos_ev_3d, pmt_point)
					inner_av_distance = light_path_cal.GetDistInInnerAV()
					av_distance = light_path_cal.GetDistInAV()
					water_distance = light_path_cal.GetDistInWater()
					transit_time = group_velocity.CalcByDistance(inner_av_distance, av_distance, water_distance)

					pmttime =  pmtCal.GetTime()

					residual = pmttime - transit_time - fVertexTime
					#residual = timeResCalc.CalcTimeResidual(pmtid,pmttime,fit_pos_3d,fVertexTime)
					
					hit_pmtid[0] = pmtCal.GetID()
					hit_residual[0] = residual

					tree.Fill()

	tree.Write()
	fout.Close()


if __name__ == "__main__":

	#read_dir = '/share/neutrino/snoplus/Data/FullFill_2p2/rat_801/bisMSB/Analysis20_bMR/ratds/run_evlist/Analysis20_bMR_r0000358052_s002_p005_runevelist_42.root'
	read_dir = ''
	file_txt_dir = ''
	save_dir = ''

	extract_data(read_dir, file_txt_dir, save_dir)





			















