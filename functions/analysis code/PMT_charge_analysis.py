
import numpy as np

def calculate_integrated_charge_per_event(npz_file, t_min=-4.0, t_max=1.0, save_output=None):
    """
    Compute the total charge acumulated by each PMT hit type per event, i.e.,
    the charge of all PMTs is the summed for each hit_type.
    
    Compute the charge acumulated by each hit PMT id of each hit_type per event.
    The hit types of interest are hit_type = 2 (pure scintillator) and hit_type = 3 
    (Cherenkov plus scintillation light).
    
    The PMTs are selected only through a time residual cut.

    Parameters*:
    -npz_file: Directory where are the files to be analyzed in npz format.
    .t_min, t_max: min and max time residual window limits to select hit PMTs
    -save_output: defaul to None. Directory where to save the analysis results.

    Returns: results (npz array with the analysis results)
    """
    print(f"Loading Files: {npz_file}")
    data = np.load(npz_file)

    # Extract the data
    mcid = data['mcid']
    evtid = data['evtid']
    hit_pmt_id = data['hit_pmt_id']
    t_res = data['time_residual_recons'] 
    qhs = data['hit_qhs']
    hit_type = data['hit_type']

    # 1. Perform the tres cut
    time_mask = (t_res >= t_min) & (t_res <= t_max)

    mcid_cut = mcid[time_mask]
    evtid_cut = evtid[time_mask]
    hit_pmt_id_cut = hit_pmt_id[time_mask]
    qhs_cut = qhs[time_mask]
    hit_type_cut = hit_type[time_mask]

    # Identify the events: Use mcid and evtid to avoid mixing the events
    unique_events, event_indices = np.unique(
        np.column_stack((mcid_cut, evtid_cut)), axis=0, return_index=True
    )

    # Structure of the output/results
    results = {
        'mcid': [],
        'evtid': [],
        'q_tot_scint_only': [],   # Carga integrada de todos los hit_type == 2 por evento
        'q_tot_cher_scint': [],   # Carga integrada de todos los hit_type == 3 por evento
        'q_per_pmt_scint_only': [],       # Carga integrada de PMTs hit_type == 2 por evento
        'q_per_pmt_cher_scint': [],       # Carga integrada de PMTs hit_type == 3 por evento
        'norm_q_scint_only': [],  # Carga integrada de PMTs hit_type == 2 por evento normalizado por nhits
        'norm_q_cher_scint': [],  # Carga integrada de PMTs hit_type == 3 por evento normalizado por nhits
        'nhits_scint_only': [],   # Número de hits/PMTs de tipo 2 en la ventana
        'nhits_cher_scint': []    # Número de hits/PMTs de tipo 3 en la ventana
    }

    # 2. Event by event analysis
    for evt_mc, evt_gtid in unique_events:
        # Mask to select event-by-event through mcid and eventID entries
        mask_evt = (mcid_cut == evt_mc) & (evtid_cut == evt_gtid)

        hit_pmt_id_evt = hit_pmt_id_cut[mask_evt]
        q_evt = qhs_cut[mask_evt]
        ht_evt = hit_type_cut[mask_evt]

        # hit_type masks
        mask_type2 = (ht_evt == 2) 
        mask_type3 = (ht_evt == 3) #| (ht_evt == 1)

        # NHits of each type within the event
        nhits_type2 = np.sum(mask_type2)
        nhits_type3 = np.sum(mask_type3)

        # Charge Integration of all PMTs wihtin and event for each hit type
        q_sum_type2 = np.sum(q_evt[mask_type2]) if np.any(mask_type2) else 0.0
        q_sum_type3 = np.sum(q_evt[mask_type3]) if np.any(mask_type3) else 0.0

        #Charge integration for each PMT hit id within an event foe each hit_type
        hit_pmtid_type2 = hit_pmt_id_evt[mask_type2]
        hit_pmtid_type3 = hit_pmt_id_evt[mask_type3]

        unique_pmtid_type2 = np.unique(hit_pmtid_type2, axis = 0)
        unique_pmtid_type3 = np.unique(hit_pmtid_type3, axis = 0)

        q_evt_type2 = q_evt[mask_type2]
        q_evt_type3 = q_evt[mask_type3]

        # Iteration over hit PMT ID for each hit_type:
        for hit_id_i in unique_pmtid_type2:
            # Mask to select only this hit PMT ID
            mask_pmtid = (hit_pmtid_type2 == hit_id_i)

            q_evt_pmtid = q_evt_type2[mask_pmtid]
            q_sum_pmt = np.sum(q_evt_pmtid)

            # Integrate the charge observed by this PMT
            results['q_per_pmt_scint_only'].append(q_sum_pmt)
            results['norm_q_scint_only'].append(q_sum_pmt/nhits_type2)

        for hit_id_i in unique_pmtid_type3:
            # Mask to select only this hit PMT ID
            mask_pmtid = (hit_pmtid_type3 == hit_id_i)

            q_evt_pmtid = q_evt_type3[mask_pmtid]
            q_sum_pmt = np.sum(q_evt_pmtid)
            
            # Save the Integrated and normalized charge observed by this PMT
            results['q_per_pmt_cher_scint'].append(q_sum_pmt)
            results['norm_q_cher_scint'].append(q_sum_pmt/nhits_type3)

        # save the results
        results['mcid'].append(evt_mc)
        results['evtid'].append(evt_gtid)
        results['q_tot_scint_only'].append(q_sum_type2)
        results['q_tot_cher_scint'].append(q_sum_type3)
        results['nhits_scint_only'].append(nhits_type2)
        results['nhits_cher_scint'].append(nhits_type3)

    # Transform to numpy arrays
    for k in results:
        results[k] = np.array(results[k])

    print(f"Events Analyzed: {len(results['mcid'])}")

    # 3. Guardar el resultado procesado si se especifica ruta de guardado
    if save_output:
        np.savez(save_output, **results)
        print(f"Cargas integradas guardadas en: {save_output}")

    return results