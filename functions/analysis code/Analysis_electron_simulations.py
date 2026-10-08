'''
Functions designed to analyze electron simulations.
The functions will read ROOT files with RATDS structure and will save the
PMT DB information and observables of interest to investigate the separation
of Cherenkov hits from scintilaltion hits by looking at the MC branches of the
photosn recorded by each PMT. Also, the Sun direction will be used to produce
the directional observable which should nor correlated with MC direction of the 
generated events.

Creation: 01/10/2026
'''

import rat
import ROOT
import numpy as np

def extract_cherenkov_scint_data(read_dir, save_dir):
    """
    Read RATDS ROOT MC simulations and extract the data at the hitPMT
    level to differentiate between Cherenkov and Scintillation signal.
    The results will be saved in npz files (numpy format)
    """
    print(f"Extracting RATDS ROOT MC Data from: {read_dir}")

    # Inicializar RAT utility y cargar la base de datos de PMTs
    util = rat.utility()
    util.LoadDBAndBeginRun()
    
    pmtinfo = util.GetPMTInfo()
    pmtCalStatus = util.GetPMTCalStatus()
    light_path = util.GetLightPathCalculator()
    group_velocity = util.GetGroupVelocity()

    Sun_Dir = rat.RAT.SunDirection
    P3D = ROOT.RAT.DU.Point3D
    psup_id = P3D.GetSystemId("innerPMT")

    # Arrays para almacenar datos del evento
    data = {
        'evtid': [], 'mcid': [],
        'energy_mc': [], 'energy_recon': [],
        'pos_mc': [], 'pos_recon': [],
        'dir_mc': [], 'sun_dir': [],
        'hit_pmt_id': [], 'hit_time': [], 'time_residual_mc': [], 'time_residual_recons': [],
        'hit_qhs': [], 'hit_qhl': [],
        'hit_type': [],  # 1 para Cherenkov, 2 para Scintillation, 0 para Otro
        'dir_angle_true_recons': [], 'dir_angle_true_mc': [], 'dir_angle_sun': []
    }
    
    # Array para la base de datos de PMTs
    pmt_db = []
    for i_pmt in range(pmtinfo.GetCount()):
        pos = pmtinfo.GetPosition(i_pmt)
        pmt_db.append([i_pmt, pos.x(), pos.y(), pos.z(), pmtinfo.GetType(i_pmt)])
    pmt_db = np.array(pmt_db)

    # Leer el archivo ROOT
    reader = ROOT.RAT.DU.DSReader(read_dir)
    
    for ievent in range(reader.GetEntryCount()):
        #print(f'reading event #{ievent}')
        rDS = reader.GetEntry(ievent)
        rMC = rDS.GetMC()
        
        # --- Obtener información MC ---
        mcid = rMC.GetMCID()
        mc_particle = rMC.GetMCParticle(0)
        energy_mc = mc_particle.GetKineticEnergy()
        #print(f'event with energy {energy_mc} (MeV)')
        
        mc_pos_vec = mc_particle.GetPosition()
        pos_mc = [mc_pos_vec.x(), mc_pos_vec.y(), mc_pos_vec.z()]
        
        mc_mom_vec = mc_particle.GetMomentum()
        # Normalizar dirección MC
        mom_mag = np.sqrt(mc_mom_vec.x()**2 + mc_mom_vec.y()**2 + mc_mom_vec.z()**2)
        dir_mc = [mc_mom_vec.x()/mom_mag, mc_mom_vec.y()/mom_mag, mc_mom_vec.z()/mom_mag] if mom_mag > 0 else [0,0,0]

        # --- Iterar sobre los eventos reconstruidos (EV) ---
        for iev in range(rDS.GetEVCount()):
            print(f'reading event #{iev}')
            rEV = rDS.GetEV(iev)
            evtid = rEV.GetGTID()
            
            # Requerir un ajuste válido
            if not rEV.FitResultExists("scintFitter"): continue
            fResult = rEV.GetFitResult("scintFitter")
            if not fResult.GetValid() or fResult.GetVertexCount() < 1: continue

            print('passed the fitters')
            
            fVertex = fResult.GetVertex(0) # Usar el primer vértice
            if not (fVertex.ContainsPosition() and fVertex.ContainsEnergy() and fVertex.ValidPosition() and fVertex.ValidEnergy()): continue
            
            energy_recon = fVertex.GetEnergy()
            fPosition = fVertex.GetPosition()
            pos_recon = [fPosition.x(), fPosition.y(), fPosition.z()]
            fVertexTime = fVertex.GetTime()
            
            # Obtener dirección del sol para el tiempo del evento
            rTime = rEV.GetUniversalTime()
            sun_dir_vec = Sun_Dir(int(rTime.GetDays()), int(rTime.GetSeconds()), int(rTime.GetNanoSeconds()))
            sun_dir = [sun_dir_vec.X(), sun_dir_vec.Y(), sun_dir_vec.Z()]
            
            # Preparar posición 3D para el cálculo del time residual
            fit_pos_3d = P3D(psup_id, pos_recon[0], pos_recon[1], pos_recon[2])
            mc_pos_3d = P3D(psup_id, pos_mc[0], pos_mc[1], pos_mc[2])
            
            # --- Iterar sobre los PMTs activados (Hits) ---
            calibratedPMTs = rEV.GetCalPMTs()
            
            # Crear un diccionario del MC PMT para clasificar rápido los fotones
            mc_pmt_dict = {}
            for i_mc_pmt in range(rMC.GetMCPMTCount()):
                mc_pmt = rMC.GetMCPMT(i_mc_pmt)
                pmt_id = mc_pmt.GetID()
                
                hit_type = 0 # Default: Otro
                
                # Revisar los PE (Photoelectrons) en este PMT para determinar su historia
                if mc_pmt.GetMCPECount() > 0:
                    pe = mc_pmt.GetMCPE(0) # Usamos el primer fotón detectado en el PMT
                    
                    # 1 = Cherenkov (bit 1 según documentación, revisa el index real si falla)
                    if pe.GetFromHistory(1): 
                        hit_type = 1
                    # 2 = Scintillation (bit 2 según documentación)
                    elif pe.GetFromHistory(2):
                        hit_type = 2
                
                mc_pmt_dict[pmt_id] = hit_type
            
            # iteration over the callibrated hit PMTs
            for i_pmt in range(calibratedPMTs.GetAllCount()):
                pmtCal = calibratedPMTs.GetAllPMT(i_pmt)
                pmt_id = pmtCal.GetID()

                if pmtCalStatus.GetHitStatus(pmtCal) != 0: continue

                # ---- Time Residual calculation using fitted and reconstructed quantities ----

                pmt_point = P3D(psup_id, pmtinfo.GetPosition(pmtCal.GetID()))

                # Reconstructed quantities
                light_path.CalcByPosition(fit_pos_3d, pmt_point)
                inner_av_distance_recons = light_path.GetDistInInnerAV()
                av_distance_recons = light_path.GetDistInAV()
                water_distance_recons = light_path.GetDistInWater()
                transit_time_recons = group_velocity.CalcByDistance(inner_av_distance_recons, av_distance_recons, water_distance_recons)

                # Simulated quantities
                light_path.CalcByPosition(mc_pos_3d, pmt_point)
                inner_av_distance_mc = light_path.GetDistInInnerAV()
                av_distance_mc = light_path.GetDistInAV()
                water_distance_mc = light_path.GetDistInWater()
                transit_time_mc = group_velocity.CalcByDistance(inner_av_distance_mc, av_distance_mc, water_distance_mc)

                pmt_time = pmtCal.GetTime()

                residual_recons = pmt_time - transit_time_recons - fVertexTime
                residual_mc = pmt_time - transit_time_mc - fVertexTime

                # Obtener la posición del PMT de la base de datos
                pmt_pos = pmtinfo.GetPosition(pmt_id)
                pmt_pos_vec = np.array([pmt_pos.x(), pmt_pos.y(), pmt_pos.z()])
                
                # Vector desde el evento reconstruido y simulado hacia el PMT
                hit_vec_recons = pmt_pos_vec - np.array(pos_recon)
                hit_vec_mc = pmt_pos_vec - np.array(pos_mc)

                hit_mag_recons = np.linalg.norm(hit_vec_recons)
                hit_mag_mc = np.linalg.norm(hit_vec_mc)

                hit_dir_recons = hit_vec_recons / hit_mag_recons if hit_mag_recons > 0 else np.zeros(3)
                hit_dir_mc = hit_vec_mc / hit_mag_mc if hit_mag_mc > 0 else np.zeros(3)
                
                # Calcular el coseno del ángulo con la dirección MC y con el Sol
                cos_theta_recons = np.dot(hit_dir_recons, dir_mc)  # based on the recons. position
                cos_theta_mc = np.dot(hit_dir_mc, dir_mc)          # based on the mc position
                cos_theta_sun = np.dot(hit_dir_recons, sun_dir)
                
                # Buscar el tipo de hit en el diccionario MC
                hit_type = mc_pmt_dict.get(pmt_id, 0)
                
                # Guardar los datos del Hit
                data['evtid'].append(evtid)
                data['mcid'].append(mcid)
                data['energy_mc'].append(energy_mc)
                data['energy_recon'].append(energy_recon)
                data['pos_mc'].append(pos_mc)
                data['pos_recon'].append(pos_recon)
                data['dir_mc'].append(dir_mc)
                data['sun_dir'].append(sun_dir)
                data['hit_pmt_id'].append(pmt_id)
                data['hit_time'].append(pmt_time)
                data['time_residual_recons'].append(residual_recons)
                data['time_residual_mc'].append(residual_mc)
                data['hit_qhs'].append(pmtCal.GetQHS())
                data['hit_qhl'].append(pmtCal.GetQHL())
                data['hit_type'].append(hit_type)
                data['dir_angle_true_recons'].append(cos_theta_recons)
                data['dir_angle_true_mc'].append(cos_theta_mc)
                data['dir_angle_sun'].append(cos_theta_sun)
                
    # Convertir a numpy arrays
    for key in data:
        data[key] = np.array(data[key])
        
    # Guardar en archivo .npz
    np.savez(save_dir, pmt_db=pmt_db, **data)
    print(f"Data saved to {save_dir}")


if __name__ == "__main__":
    read_dir = "/lstore/sno/joankl/solar_analysis/mc_data/2p2_ppo/electrons/output_root/e_2.5MeV_0.root"
    save_dir = "/lstore/sno/joankl/solar_analysis/mc_data/2p2_ppo/electrons/results_npz/e_2.5MeV_0.npz"
    extract_cherenkov_scint_data(read_dir, save_dir)