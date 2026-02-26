'''
Function designed to read ntuples of real data and perform:
1) A Hotspot and atmospheric vetoing.
2) A coincidence analysis between a prompt and a delay events.

The vetoing algorithm will save the number of found events and followers in a vetoing
dictionary. Also, it will save the runID and the GTID.
The coincidence analysis will save the GTID of the event as well as the energy, position,
and time of the found prompt and delay.
'''

import uproot
import numpy as np
import glob
import re
import os
import pickle


# =================== Auxiliar Functions ===================

def read_files_txt(file_txt_dir):
    '''Read and file.txt and returns a list with the names of the files to be readen'''
    with open(file_txt_dir, "r", encoding="utf-8") as f:
        file_names = [line.strip() for line in f if line.strip()]
    return file_names

def natural_order(file):
    """Function to arange the file names by order if run and subrun"""
    return [int(text) if text.isdigit() else text.lower() for text in re.split('(\d+)', file)]

# ==========================================================

def veto_and_coincidence_analysis(read_dir, file_txt_dir, save_dir):

    '''
    Parameters:

    -read_dir: Directory where the real_data.ntuple.root files are.
    -file_txt_dit: Name of the file to be readen. The function can be implemented
                   within a loop of txt files which contain a list of the file_names
                   to be readen.
    -save_dir: Directory to save the coincidence analysis results.
    '''

    #Ensure that save directory exist. If not, then create it
    os.makedirs(save_dir, exist_ok=True)

     # ---- List of variables to read -----
    var_read_name_list = ['eventID', 'runID', 'scintFit', 'fitValid', 'dcFlagged',
                           'nhits', 'energy', 'clockCount50', 'posr_av', 'posx', 'posy', 'posz_av']


    # ---- List of variables to save ----
    # The names on the list not necessarily correspond to the ntuples entries.
    # A dictionary for prompt and delayed candidates will be created and those entries
    # should correspond to the var_save_name_list.
    data_var_save_name_list = ['runID', 'eventID', 'energy', 'posx', 'posy', 'posz_av', 'nhits', 'time']
    coincidence_var_save_name_list = ['runID', 'eventID', 'energy', 'posx', 'posy', 'posz', 'time']

    # ----- Data Dictionary -----
    # This dictionary will save the data extracted from the loop over the files.
    # It will be used to perform the coincidence analysis
    data_dict = {var: [] for var in data_var_save_name_list}

    # ----- HS and atmospheric Dictionary -----
    hs_dict = { 'counter': [],
                'eventID': [],
                'runID': []}

    atm_dict = {'counter': [],
                'eventID': [],
                'runID': []}

    # ----- Prompt and delayed data dictionaties -----
    prompt_coinc_dict = {var: [] for var in coincidence_var_save_name_list}
    delay_coinc_dict = {var: [] for var in coincidence_var_save_name_list}


    # ------ Files reader ------
    fname_list = read_files_txt(file_txt_dir)
    fname_list.sort(key=natural_order)

    for i_dx, fname_i in enumerate(fname_list):
        full_fdir_i = os.path.join(read_dir, fname_i)

        try:
            print(f'Reading file: {full_fdir_i}')
            file_i = uproot.open(full_fdir_i)
            output = file_i['output;1']
        except Exception as e:
            print(f"Error while trying to read {fname_i}: {e}")
            continue


        # Create a temporal dictionary with the useful ntuple variables to perform the analysis
        temp_vars = {}

        for var in var_read_name_list:
            if var == 'dcFlagged':
                temp_vars[var] = np.array(output[var]).astype(np.uint64)
            else:
                temp_vars[var] = np.array(output[var])

        # ============ Cut Conditions ============

        # ---- General Cut Conditions overall dataset ----
        valid_condition = (temp_vars['scintFit'] & temp_vars['fitValid'])

        nhits_min = 20
        nhits_condition = (temp_vars['nhits'] >= nhits_min)

        energy_cut = 1.0 #MeV
        energy_condition = (temp_vars['energy'] >= energy_cut)

        posr_cut = 5500.0 # mm
        posr_condition = (temp_vars['posr_av'] <= posr_cut)

        mask_cut = 0xD82100000162C6
        dcflag_condition = ((int(mask_cut) & temp_vars['dcFlagged']) == int(mask_cut))

        general_condition = valid_condition & nhits_condition & energy_condition & posr_condition & dcflag_condition

        # Apply general cuts
        for var in var_read_name_list:
            temp_vars[var] = temp_vars[var][general_condition]

        # Transform ClockCount50 into time
        time = (temp_vars['clockCount50']*20)/1000 # in microseconds
        temp_vars['time'] = time

        # Save the filtered result on the data_dict out of the file loop:
        for var in data_var_save_name_list:
            data_dict[var].append(temp_vars[var].tolist())

    # Flat the sublists of the data_dict
    for var in data_var_save_name_list:
        data_dict[var] = [x for sublist in data_dict[var] for x in sublist]

    # Extract the observables of interest from the data_dict
    energy = np.array(data_dict['energy'])
    time = np.array(data_dict['time'])
    nhits = np.array(data_dict['nhits'])
    posx, posy, posz = np.array(data_dict['posx']), np.array(data_dict['posy']), np.array(data_dict['posz_av'])
    runID = np.array(data_dict['runID'])
    eventID = np.array(data_dict['eventID'])

    #print(f'nhits data = {nhits}')

    N_evs = energy.shape[0]

    # =========== Vetoing of Hotspots and Atmospherics ===========
    print('Performin Vetoing of Hotspots and Atmospherics')

    hs_counter = 0
    atm_counter = 0

    # ---- Cut Conditions ----
    nhits_hs_number = 1500  #nhits above this value is a hotspot suspected
    nhits_atm_number = 3000 #nhits above this value is an atmospheric suspected
    dt_hs_remove = 1e6      # Remove 1s of data after a hs
    dt_atm_remove = 20*1e6  # Remove 20s of data after an atmospheric

    nhits_hs_condition = (nhits > nhits_hs_number)
    nhits_atm_condition = (nhits > nhits_atm_number)

    #position condition for HS events
    pos_condition_1 = (posx > 500) & (posx < 900) & (posy > -900) & (posy < -500) & (posz > 5300)
    pos_condition_2 = (posx > 0) & (posx < 500) & (posy > -1200) & (posy < -800) & (posz > 5300)
    pos_condition = (pos_condition_1) | (pos_condition_2)

    hs_condition = nhits_hs_condition & pos_condition

    # Index of the found prompt and followers hs and atm event to be removed from the dataset
    # before performing the coincidence analysis.
    hs_index_to_remove = []
    atm_index_to_remove = []

    hs_ev_index = np.where(hs_condition)[0]
    atm_ev_index = np.where(nhits_atm_condition)[0]

    #Evaluates if there are repeated index of atm_ev_index in hs_ev_index.
    if len(hs_ev_index) >= len(atm_ev_index):
        same_index_condition = np.in1d(atm_ev_index, hs_ev_index)
        index_to_del = np.where(same_index_condition)[0]
        atm_ev_index = np.delete(atm_ev_index, index_to_del)

    if len(atm_ev_index) > len(hs_ev_index):
        same_index_condition = np.in1d(hs_ev_index, atm_ev_index)
        index_to_del = np.where(same_index_condition)[0]
        atm_ev_index = np.delete(atm_ev_index, index_to_del)


    # Save the corresponding GTID and runID of the suspected HS and atm
    hs_dict['eventID'].append(eventID[hs_ev_index].tolist())
    hs_dict['runID'].append(runID[hs_ev_index].tolist())

    atm_dict['eventID'].append(eventID[atm_ev_index].tolist())
    atm_dict['runID'].append(runID[atm_ev_index].tolist())

    # ------ Loop on events to save the followers ------
    # (1) HS:
    # First verify if there are Hotspots events before saving followers
    if len(hs_ev_index) > 0:
        print('Extracting HS followers')
        hs_counter += len(hs_ev_index)
        hs_dict['counter'].append(hs_counter)

        # Iteration over HS followers
        for i_dx in hs_ev_index:
            hs_idx = i_dx
            hs_next_idx = hs_idx + 1

            hs_index_to_remove.append(hs_idx)

            time_hs = time[hs_idx]

            try:
                runID_next_hs = runID[hs_next_idx]
                eventID_next_hs = eventID[hs_next_idx]
                time_next_hs = time[hs_next_idx]
                dt = time_next_hs - time_hs

            except IndexError:
                continue

            while (dt > 0) and (dt <= dt_hs_remove):

                # Save the followers eventID and runID into dictionary
                hs_dict['runID'].append(runID_next_hs)
                hs_dict['eventID'].append(eventID_next_hs)

                hs_index_to_remove.append(hs_next_idx)

                #Refresh the follower index
                hs_next_idx += 1

                try:
                    runID_next_hs = runID[hs_next_idx]
                    eventID_next_hs = eventID[hs_next_idx]
                    time_next_hs = time[hs_next_idx]
                    dt = time_next_hs - time_hs

                except IndexError:
                    continue
    #else:
        #continue

    # (2) atm:
    # First verify if there are atmospheric events before saving followers
    if len(atm_ev_index) > 0:
        print('Extracting atm followers')
        atm_counter += len(atm_ev_index)
        atm_dict['counter'].append(atm_counter)

        # Iteration over HS followers
        for i_dx in atm_ev_index:
            atm_idx = i_dx
            atm_next_idx = atm_idx + 1

            atm_index_to_remove.append(atm_idx)

            time_atm = time[atm_idx]

            try:
                runID_next_atm = runID[atm_next_idx]
                eventID_next_atm = eventID[atm_next_idx]

                time_next_atm = time[atm_next_idx]
                dt = time_next_atm - time_atm

            except IndexError:
                continue

            while (dt > 0) and (dt <= dt_atm_remove):

                # Save the followers eventID and runID into dictionary
                atm_dict['runID'].append(runID_next_atm)
                atm_dict['eventID'].append(eventID_next_atm)

                atm_index_to_remove.append(atm_next_idx)

                #Refresh the follower index
                atm_next_idx += 1

                try:
                    runID_next_atm = runID[atm_next_idx]
                    eventID_next_atm = eventID[atm_next_idx]

                    time_next_atm = time[atm_next_idx]
                    dt = time_next_atm - time_atm

                except IndexError:
                    continue
    #else:
        #continue

    print(f'Number of HS events = {hs_counter}')
    print(f'Number of atm events = {atm_counter}')

    # =========== Coincidence finder ===========
    print('Performing coincidence analysis')

    # ---- Cut the Hotspots and Atmospheric prompt and followers ----
    full_index_to_remove = np.array(hs_index_to_remove + atm_index_to_remove, dtype = int)
    full_index_to_remove = np.unique(full_index_to_remove)

    #print(f'index to remove = {full_index_to_remove}')

    energy = np.delete(energy, full_index_to_remove)
    time = np.delete(time, full_index_to_remove)
    nhits = np.delete(nhits, full_index_to_remove)
    posx = np.delete(posx, full_index_to_remove)
    posy = np.delete(posy, full_index_to_remove)
    posz = np.delete(posz, full_index_to_remove)
    runID = np.delete(runID, full_index_to_remove)
    eventID = np.delete(eventID, full_index_to_remove)

    # ---- Coincidence Cuts ----
    energy_prompt_inf_cut = 1.0
    energy_prompt_sup_cut = 8.0

    energy_delay_inf_cut = 1.0
    energy_delay_sup_cut = 4.0

    dt_inf_cut = 0.5
    dt_sup_cut = 1000

    dr_inf_cut = 0
    dr_sup_cut = 2500

    # ---- Loop on events ----

    N_evs = energy.shape[0]

    #print(f'Number of events = {N_evs}')

    for ev_i in range(N_evs):

        prompt_idx = ev_i
        delay_idx = prompt_idx + 1

        # Suspected Prompt observables
        runID1 = runID[prompt_idx]
        eventID1 = eventID[prompt_idx]
        energy1 = energy[prompt_idx]
        posx1 = posx[prompt_idx]
        posy1 = posy[prompt_idx]
        posz1 = posz[prompt_idx]
        t1 = time[prompt_idx]

        # Suspected delay observables
        try:
            runID2 = runID[delay_idx]
            eventID2 = eventID[delay_idx]
            energy2 = energy[delay_idx]
            posx2 = posx[delay_idx]
            posy2 = posy[delay_idx]
            posz2 = posz[delay_idx]
            t2 = time[delay_idx]

            # Evaluate dt and dr between events
            dx = posx2 - posx1
            dy = posy2 - posy1
            dz = posz2 - posz1

            dr = np.sqrt(dx**2 + dy**2 + dz**2)
            dt = t2 - t1

        except IndexError:
            continue

        # Look for the coincidence
        while (dt > dt_inf_cut) and (dt <= dt_sup_cut):


            if (dr >= dr_inf_cut) and (dr <= dr_sup_cut) and (dt > dt_inf_cut) and (dt <= dt_sup_cut) and (energy2 >= energy_delay_inf_cut) and (energy2 <= energy_delay_sup_cut):
                print(f'Pair found with dr = {dr} and dt = {dt}')

                prompt_obs_lst = [runID1, eventID1, energy1, posx1, posy1, posz1, t1]
                delay_obs_lst = [runID2, eventID2, energy2, posx2, posy2, posz2, t2]

                # Save the coincidences in the coincidence dictionaries
                for idy, var_name in enumerate(coincidence_var_save_name_list):
                    prompt_coinc_dict[var_name].append(prompt_obs_lst[idy])
                    delay_coinc_dict[var_name].append(delay_obs_lst[idy])

                break
            delay_idx += 1

            # Look for the next suspected delay
            try:
                runID2 = runID[delay_idx]
                eventID2 = eventID[delay_idx]
                energy2 = energy[delay_idx]
                posx2 = posx[delay_idx]
                posy2 = posy[delay_idx]
                posz2 = posz[delay_idx]
                t2 = time[delay_idx]

                # Evaluate dt and dr between events
                dx = posx2 - posx1
                dy = posy2 - posy1
                dz = posz2 - posz1

                dr = np.sqrt(dx**2 + dy**2 + dz**2)
                dt = t2 - t1

            except IndexError:
                break

    # ====== Save the Dictionaries ======
    dict_save_list = [hs_dict, atm_dict, prompt_coinc_dict, delay_coinc_dict]
    dict_name_save = ['hs_dict', 'atm_dict', 'prompt_coinc_dict', 'delay_coinc_dict']

    for dict_i, dict_name in zip(dict_save_list, dict_name_save):
        with open(save_dir + dict_name + '.pkl', 'wb') as f:
            pickle.dump(dict_i, f)

    return print('Analysis Done!')


if __name__ == '__main__':

    read_dir = '/share/neutrino/snoplus/Data/FullFill_2p2/rat_801/bisMSB/Analysis15/'
    file_txt_dir = '/lstore/sno/joankl/solar_analysis/real_data/bisMSB/Analysis15/file_name_list/sublist_18.txt'
    save_dir = '/lstore/sno/joankl/solar_analysis/real_data/bisMSB/Analysis15/proof/'

    veto_and_coincidence_analysis(read_dir, file_txt_dir, save_dir)
