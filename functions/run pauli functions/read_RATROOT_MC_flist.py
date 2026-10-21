"""
Python script designed to read a RATDS ROOT files from MC data
and extract observables of interest. 
The function should receive a list of input.root files
and will iterae over them to save the info of various
RATDS in the Branches.

The script  perform the cuts:

	- ScintFit true;
	- Valid Fit result;
	- Vertex Count > 1
	- Valid Vertex Reconstructe Quantites
	- Energy cut above 2.5;
The script employs the pile-up analysis on MC through VertexCount.
Te script will generate a summary of the analysis in npz format.
"""

import os
import sys
import numpy as np
import ROOT
import rat

def extract_data(read_dir, file_txt_dir, save_dir):
    print(f"Reading file list from: {file_txt_dir}")
    with open(file_txt_dir, "r", encoding="utf-8") as f:
        flist = [line.strip() for line in f if line.strip()]

    if not flist:
        print(f"No files found in list: {file_txt_dir}")
        return

    # Inicializar herramientas de RAT
    util = rat.utility()
    util.LoadDBAndBeginRun()

    SunDir = rat.RAT.SunDirection
    light_path_cal = util.GetLightPathCalculator()
    group_velocity = util.GetGroupVelocity()
    pmtinfo = util.GetPMTInfo()
    pmtCalStatus = util.GetPMTCalStatus()

    P3D = ROOT.RAT.DU.Point3D
    psup_id = P3D.GetSystemId("innerPMT")
    av_id = P3D.GetSystemId("av")

    calibrator = util.GetReconCalibrator()
    MATERIAL_NAME = "labppo_2p2_scintillator"
    CORRECTION_VER = 3
    IS_DATA = False

    # Guardar metadata de PMTs
    pmt_db = []
    for i_pmt in range(pmtinfo.GetCount()):
        pos = pmtinfo.GetPosition(i_pmt)
        pmt_db.append([i_pmt, pos.x(), pos.y(), pos.z(), pmtinfo.GetType(i_pmt)])
    pmt_db = np.array(pmt_db)

    # Diccionario para almacenar los datos
    data = {
        # A nivel de evento
        'mcid': [],
        'evtid': [],
        'energy_mc': [],
        'parent_energy': [],
        'energy_recon': [],
        'energy_corr': [],
        'pos_mc': [],
        'pos_recon': [],
        'momentum_mc': [],
        'sun_dir': [],
        'evt_time_day': [],
        'evt_time_sec': [],
        'evt_time_nsec': [],
        # A nivel de hit/PMT
        'hit_mcid': [],
        'hit_evtid': [],
        'hit_pmtid': [],
        'hit_time': [],
        'hit_residual': [],
        'hit_qhs': [],
        'hit_qhl': []
    }

    for input_file_i in flist:
        full_path = os.path.join(read_dir, input_file_i)

        print(f"Opening RATDS: {full_path}")
        reader = ROOT.RAT.DU.DSReader(full_path)

        for ievent in range(reader.GetEntryCount()):
            rDS = reader.GetEntry(ievent)
            rMC = rDS.GetMC()

            mc_id = rMC.GetMCID()
            parent_e = rMC.GetMCParent(0).GetKineticEnergy() if rMC.GetMCParentCount() > 0 else 0.0
            particle_0 = rMC.GetMCParticle(0)
            e_mc = particle_0.GetKineticEnergy()
            pos_mc_vec = particle_0.GetPosition()
            pos_mc = [pos_mc_vec.x(), pos_mc_vec.y(), pos_mc_vec.z()]
            mom_mc_vec = particle_0.GetMomentum()
            mom_mc = [mom_mc_vec.x(), mom_mc_vec.y(), mom_mc_vec.z()]

            for iev in range(rDS.GetEVCount()):
                rEV = rDS.GetEV(iev)
                evt_id = rEV.GetGTID()

                if not rEV.FitResultExists("scintFitter"):
                    continue

                fResult = rEV.GetFitResult("scintFitter")
                if not fResult.GetValid() or fResult.GetVertexCount() < 1:
                    continue

                for ivtx in range(fResult.GetVertexCount()):
                    fVertex = fResult.GetVertex(ivtx)
                    if not (fVertex.ContainsPosition() and fVertex.ContainsEnergy() and
                            fVertex.ValidPosition() and fVertex.ValidEnergy()):
                        continue

                    e_recon = fVertex.GetEnergy()
                    if e_recon < 2.5:  # Corte de energía
                        continue

                    fPos = fVertex.GetPosition()
                    pos_recon = [fPos.x(), fPos.y(), fPos.z()]

                    # Corrección de energía RTF
                    pos_p3d_av = P3D(av_id, pos_recon[0], pos_recon[1], pos_recon[2])
                    e_corr = calibrator.CalibrateEnergyRTF(IS_DATA, e_recon, pos_p3d_av, MATERIAL_NAME, CORRECTION_VER)

                    # Tiempos y Dirección del Sol
                    rTime = rEV.GetUniversalTime()
                    t_day, t_sec, t_ns = int(rTime.GetDays()), int(rTime.GetSeconds()), int(rTime.GetNanoSeconds())
                    sun_vec = SunDir(t_day, t_sec, t_ns)
                    sun_dir = [sun_vec.X(), sun_vec.Y(), sun_vec.Z()]

                    # Guardar variables a nivel de evento
                    data['mcid'].append(mc_id)
                    data['evtid'].append(evt_id)
                    data['parent_energy'].append(parent_e)
                    data['energy_mc'].append(e_mc)
                    data['energy_recon'].append(e_recon)
                    data['energy_corr'].append(e_corr)
                    data['pos_mc'].append(pos_mc)
                    data['pos_recon'].append(pos_recon)
                    data['momentum_mc'].append(mom_mc)
                    data['sun_dir'].append(sun_dir)
                    data['evt_time_day'].append(t_day)
                    data['evt_time_sec'].append(t_sec)
                    data['evt_time_nsec'].append(t_ns)

                    fVertexTime = fVertex.GetTime()
                    fit_pos_3d = P3D(psup_id, pos_recon[0], pos_recon[1], pos_recon[2])

                    # Loop de hits PMT
                    calPMTs = rEV.GetCalPMTs()
                    for iPMT in range(calPMTs.GetAllCount()):
                        pmtCal = calPMTs.GetAllPMT(iPMT)
                        if pmtCalStatus.GetHitStatus(pmtCal) != 0:
                            continue

                        pmt_id = pmtCal.GetID()
                        pmt_point = P3D(psup_id, pmtinfo.GetPosition(pmt_id))

                        light_path_cal.CalcByPosition(fit_pos_3d, pmt_point)
                        inner_d = light_path_cal.GetDistInInnerAV()
                        av_d = light_path_cal.GetDistInAV()
                        water_d = light_path_cal.GetDistInWater()
                        t_transit = group_velocity.CalcByDistance(inner_d, av_d, water_d)

                        pmt_time = pmtCal.GetTime()
                        residual = pmt_time - t_transit - fVertexTime

                        # Guardar observables a nivel de hit (incluyendo carga)
                        data['hit_mcid'].append(mc_id)
                        data['hit_evtid'].append(evt_id)
                        data['hit_pmtid'].append(pmt_id)
                        data['hit_time'].append(pmt_time)
                        data['hit_residual'].append(residual)
                        data['hit_qhs'].append(pmtCal.GetQHS())
                        data['hit_qhl'].append(pmtCal.GetQHL())

    # Convertir a arrays de numpy
    for k in data:
        data[k] = np.array(data[k])

    os.makedirs(os.path.dirname(save_dir), exist_ok=True)
    np.savez_compressed(save_dir, pmt_db=pmt_db, **data)
    print(f"Datos guardados exitosamente en: {save_dir}")

if __name__ == "__main__":
    if len(sys.argv) == 4:
        extract_data(sys.argv[1], sys.argv[2], sys.argv[3])




			















