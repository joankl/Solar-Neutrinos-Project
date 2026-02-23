'''
Python Script designed to read data and apply the Daniel's energy callibrator.
For more information look at https://snopl.us/docs/rat/user_manual/html/reconcalibrator_usage.html.
To use the energy calibration we need the reconstructed energy and the 3D position in AV coordinates
and with the AV offset correction 
'''

import ROOT
import rat
import glob
import os
import numpy as np


# ---- Load Data - It should be the reconstructed energy and position (av coordinates) ----
fdir_load = '/lstore/sno/joankl/solar_analysis/real_data/bisMSB/Analysis20_bMR/results/resume_files/'
fdir_save = '/lstore/sno/joankl/solar_analysis/real_data/bisMSB/Analysis20_bMR/results/resume_files/'

energy = np.load(fdir_load + 'energy.npy')
posx = np.load(fdir_load + 'posx.npy')
posy = np.load(fdir_load + 'posy.npy')
posz = np.load(fdir_load + 'posz_av.npy')

# ======= Start Using RAT =======

# Define material and version of the correction
MATERIAL_NAME = "labppo_2p2_bismsb_2p2_scintillator"
CORRECTION_VER = 3   # VER = 2 for bisMSB data/MC
IS_DATA = True   # False if MC

# Load RATDB tables into ReconCalibrator (VERY IMPORTANT!)
du = rat.utility()
du.LoadDBAndBeginRun()

# Load the Energy Calibrator
calibrator = du.GetReconCalibrator()

# offset correctly: it’s easiest if we just give the event positions in AV coords already.

P3D = ROOT.RAT.DU.Point3D
av_id = P3D.GetSystemId("av")

# ---- Loop over the events to apply the energy calibration ----

n_evs = len(energy)
energy_corr = np.empty(n_evs, dtype = np.float64)

for i in range(n_evs):
	position = P3D(av_id, 
		float(posx[i]), 
		float(posy[i]), 
		float(posz[i])
		)

	energy_corr[i] = calibrator.CalibrateEnergyRTF(IS_DATA, 
		float(energy[i]),
		position,
		MATERIAL_NAME,
		CORRECTION_VER
		)


# Save the corrected energy
np.save(fdir_save + 'energy_corrected.npy', energy_corr)
