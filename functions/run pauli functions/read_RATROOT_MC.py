"""
Python script designed to read a RATDS ROOT file from MC data
and extract observables of interest. 
The script  perform the cuts:

	- ScintFit true;
	- Valid Fit result;
	- Vertex Count > 1
	- Valid Vertex Reconstructe Quantites
	- Energy cut above 2.5;
The script employs the pile-up analysis on MC through VertexCount
"""

import rat
import ROOT
from array import array
import glob


def extract_data(read_dir, save_dir):

	"""
	Wrapper function to read and save the PMT info and events info.
	Parameters:

	- read_dir: Directory where the RATDS ROOT files are.
	- save_dir: Directory to save the ouput files
	-
	"""

	print('Extracting RATDS ROOT MC Data')

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

	print('creating TTree Branches')
	tree = ROOT.TTree("T", "output_summary")

	energy = array('d', [0.0])
	tree.Branch('energy', energy, 'energy/D')

	energy_mc = array('d', [0.0])
	tree.Branch('energy_mc', energy_mc, 'energy_mc/D')

	#energy_corr = array('d', [0.0])
	#tree.Branch('energy_corr', energy_corr, 'energy_corr/D')

	evtid = array('i',[0])
	tree.Branch('evtid',evtid,'evtid/I')

	mcid = array('i',[0])
	tree.Branch('mcid',mcid,'mcid/I')

	pos_xyz = array('d',3*[0.0])
	tree.Branch('position',pos_xyz, 'position[3]/D')

	pos_xyz_mc = array('d',3*[0.0])
	tree.Branch('position_mc',pos_xyz, 'position_mc[3]/D')

	sun_dir = array('d',3*[0.0])
	tree.Branch('sun_dir',sun_dir, 'sun_dir[3]/D')

	momentum_mc = array('d',3*[0.0])
	tree.Branch('momentum_mc',momentum_mc, 'momentum_mc[3]/D')

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
	timeResCalc = util.GetTimeResidualCalculator()

	# Point3D with PSUP reference
	P3D = ROOT.RAT.DU.Point3D
	psup_id = P3D.GetSystemId("innerPMT") 

	# Energy Callibrator Definitions 
	#print ('Takinkg Energy Calibrator Tool')
	#calibrator = util.GetReconCalibrator()
	#MATERIAL_NAME = "labppo_2p2_bismsb_2p2_scintillator"
	#CORRECTION_VER = 3   # VER = 2 for bisMSB data/MC
	#IS_DATA = True   # False if MC
	#P3D = ROOT.RAT.DU.Point3D
	#av_id = P3D.GetSystemId("av")

	# ====== Reading file root info ======
	print('Getting in reader ...')
	reader = ROOT.RAT.DU.DSReader(read_dir)
	print(f'reader here: {reader}', flush=True)

	# Loop on entries in root file
	for ievent in range(0, reader.GetEntryCount()):
		print(f'reading root file at entry {ievent}')
		rDS = reader.GetEntry(ievent)

		# ------ Read and Save MC Info ------
		rMC = rDS.GetMC()                                   # MC Branch
	
		mcid[0] = rMC.GetMCID()
		print(f'MCID = {mcid[0]}')

		energy_mc[0] = rMC.GetMCParticle(0).GetKineticEnergy()

		mc_pos = rMC.GetMCParticle(0).GetPosition()
		pos_xyz_mc[0] = mc_pos.x()
		pos_xyz_mc[1] = mc_pos.y()
		pos_xyz_mc[2] = mc_pos.z()

		mc_mom = rMC.GetMCParticle(0).GetMomentum()
		momentum_mc[0] = mc_mom.x()
		momentum_mc[1] = mc_mom.y()
		momentum_mc[2] = mc_mom.z()

		
		# Loop on each triggered event in the entry
		for iev in range(0, rDS.GetEVCount()):
			print(f'reading triggered event at entry {iev}')

			# Select the Getters, Apply Cuts, and Save Data

			rEV = rDS.GetEV(iev)                                # EV Brach
			evtid[0] = rEV.GetGTID()  # GTID

			if not rEV.FitResultExists("scintFitter"):
				continue


			fResult = rEV.GetFitResult("scintFitter")           # Fitter Branch
			if not fResult.GetValid():
				continue

			if fResult.GetVertexCount() < 1:
				continue

			for ivtx in range(0, fResult.GetVertexCount()):

				fVertex = fResult.GetVertex(ivtx)

				if not(fVertex.ContainsPosition() and fVertex.ContainsEnergy() and fVertex.ValidPosition() and fVertex.ValidEnergy()):
					continue

				fPosition = fResult.GetVertex(ivtx).GetPosition()      # Use the ith fitted position
				fEnergy = fResult.GetVertex(ivtx).GetEnergy()          # Grab Reconstructed Energy

				energy_inf_cut = 2.5
				if fEnergy < energy_inf_cut:
					continue


				# ------- Fill Branches -------

				# Reconstructed Position
				pos_xyz[0] = fPosition.x() 
				pos_xyz[1] = fPosition.y()
				pos_xyz[2] = fPosition.z()

				#position = P3D(fPosition.x(), fPosition.y(), fPosition.z())

				# Reconstructed Time
				rTime = rEV.GetUniversalTime()                      # Event Time Branch
				evt_time_day[0] = rTime.GetDays() 
				evt_time_sec[0] = rTime.GetSeconds()
				evt_time_nsec[0] = rTime.GetNanoSeconds()

				# Compute Sun Direction
				sun_dir_vector = SunDir(int(rTime.GetDays()), int(rTime.GetSeconds()), int(rTime.GetNanoSeconds()))

				sun_dir[0] = sun_dir_vector.X()
				sun_dir[1] = sun_dir_vector.Y()
				sun_dir[2] = sun_dir_vector.Z()

				# Reconstructed Energy
				energy[0] = fVertex.GetEnergy()
				#print(f'Reconstucted event energy: {energy[0]} (MeV)')

				fVertexTime =fVertex.GetTime()
				fit_pos_3d = P3D(psup_id, fPosition.x(), fPosition.y(), fPosition.z())
				
				# Loop over the hits
				calibratedPMTs = rEV.GetCalPMTs()                   # PMTs Brach
				for iPMT in range(0, calibratedPMTs.GetAllCount()):
					pmtCal = calibratedPMTs.GetAllPMT(iPMT)
					pmtid = pmtCal.GetID()
					pmttime =  pmtCal.GetTime()

					residual = timeResCalc.CalcTimeResidual(pmtid,pmttime,fit_pos_3d,fVertexTime)
					hit_pmtid[0] = pmtid
					hit_residual[0] = residual

					tree.Fill()

	tree.Write()
	write_pmt_info(fout)
	fout.Close()


if __name__ == "__main__":

	#read_dir = '/share/neutrino/snoplus/Data/FullFill_2p2/rat_801/bisMSB/Analysis20_bMR/ratds/run_evlist/Analysis20_bMR_r0000358052_s002_p005_runevelist_42.root'
	read_dir = '/share/neutrino/snoplus/MonteCarlo/FullFill_bismsb_801/ScintFit_2p2ppo_2p2bismsB8_Solar_NueRun/ratds/ScintFit_2p2ppo_2p2bismsbB8_Solar_NueRun_r354099_s0_p0.root'
	save_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/bisMSB/B8_solar_Nue/proof/x.root'

	extract_data(read_dir, save_dir)
	#flist = glob.glob(read_dir)

	#for i_dx, f_in in enumerate(flist):
		#print(f'reading file {f_in}')
		#file_out_name = f'analysis_ratds_sum_{i_dx}.root'
		#extract_data(f_in, save_dir + file_out_name)





			















