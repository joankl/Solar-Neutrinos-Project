'''                                                                                                                                                               
Function to read the ntuples.root of MC data. It will extract the variables of 
interest and save it in separeted numpy arrays.
This script makes use of the RAT ROOT Tools to perform:

- Energy-Position dependent Corrections

Then, the script will only work for RAT above v8

Edit Dates:
- 23/01/2026 Include RAT code to implement the energy correction
- 22/05/2026 Include ITR Cut above 0.2
'''

import numpy as np
import uproot
import glob
import re
import os
import rat
import ROOT


def extract_data(read_dir, file_txt_dir, save_dir):

    '''
    Function designed to read ntuples.root files. It will take a list of the ntuples.root
    file names and will save the analysis result in numpy files.

    Parameters:
    - read_dir: directory where the files with data are.
    - file_txt_dir:  directory where the files.txt with the ntuples file names are
    - save_dir: directory where the analysis files will be saved
    '''
    
    #Ensure that save directory exiss by creating it
    os.makedirs(save_dir, exist_ok=True)

    # Load RATDB tables into ReconCalibrator (VERY IMPORTANT!)
    du = rat.utility()
    du.LoadDBAndBeginRun()

    # ---------------- Subfunciones auxiliares ----------------
    def read_files_txt(file_txt_dir):
        """Lee un archivo .txt y devuelve una lista de nombres de archivo."""
        with open(file_txt_dir, "r", encoding="utf-8") as f:
            file_names = [line.strip() for line in f if line.strip()]
        return file_names

    def orden_natural(archivo):
        """Función para ordenar archivos naturalmente por número de run/subrun."""
        return [int(texto) if texto.isdigit() else texto.lower() for texto in re.split('(\d+)', archivo)]

    def extract_subrunID(filename_list):
        """Extrae el subrunID de cada archivo en la lista."""
        subrunID_list = np.array([], dtype=np.int16)
        for filename_i in filename_list:
            base_filename = os.path.basename(filename_i)
            match = re.search(r'_s(\d+)_', base_filename)
            if match:
                subrun_number = int(match.group(1))
                subrunID_list = np.append(subrunID_list, subrun_number)
        return subrunID_list

    def posr_cal(x, y, z):
        """Recalcula posr usando la corrección en z."""
        return np.sqrt(x**2 + y**2 + z**2)

    def ReconCalibrator(posx, posy, posz, energy):
        '''
        Function designed to return the energy corrected from
        Daniel Coockman energy callibrator

        parameters:
        - posx,y,x: numpy array coordinates of the reconstructed event position
        - energy: numpy array reconstructed event energy

        return:
        numpy array with the corrected energy of the events
        '''

        # ==== Define material and version of the correction ====
        MATERIAL_NAME = "labppo_2p2_bismsb_2p2_scintillator"
        CORRECTION_VER = 3   # VER = 3 for bisMSB data/MC
        IS_DATA = False   # False if MC

        du = rat.utility()

        # Load the Energy Calibrator
        calibrator = du.GetReconCalibrator()

        # offset correctly: it’s easiest if we just give the event positions in AV coords already.
        P3D = ROOT.RAT.DU.Point3D
        av_id = P3D.GetSystemId("av")

        # --- Loop over the events to apply the energy calibration ---
        n_evs = len(energy)
        energy_corr = np.empty(n_evs, dtype = np.float64)

        for i in range(n_evs):
            position = P3D(av_id, 
                float(posx[i]), 
                float(posy[i]), 
                float(posz[i]))

            energy_corr[i] = calibrator.CalibrateEnergyRTF(IS_DATA, 
                float(energy[i]),
                position,
                MATERIAL_NAME,
                CORRECTION_VER)

        return energy_corr

    # ============= List of parameters to read =============

    # ------- List of the variables to filter the data -------
    var_read_name_list = ['evIndex', 'scintFit', 'fitValid','nhits', 
                          'energy', 'posr_av', 'posx', 'posy', 'posz', 'posz_av', 
                          'parentKE1', 'mcke1', 'itr']

    # ------- List of the variables to save -------
    var_save_name_list = ['energy', 'posr_av','posx', 'posy', 'posz_av', 'parentKE1', 'mcke1']

    # Diccionary to save the acummulated  data
    data_dict = {var: np.array([]) for var in var_save_name_list + ['n_init_evs']}

    # ============= Files reader =============
    file_name_list = read_files_txt(file_txt_dir)
    file_name_list.sort(key=orden_natural)

    print(f'file name list {file_name_list}')

    for i_dx, file_name_i in enumerate(file_name_list):
        full_dir_file_i = os.path.join(read_dir, file_name_i)
        print(f'Reading file: {full_dir_file_i}')

        try:
            file_i = uproot.open(full_dir_file_i)
            output = file_i['output;1']
        except Exception as e:
            print(f"Error while trying to read {file_name_i}: {e}")
            continue

        temp_vars = {} # Temporal variables dictionary for the file in loop

        # Extract the Observables from the ntuples
        for var in var_read_name_list:
            temp_vars[var] = np.array(output[var])

        # simulated events condition
        gen_ev_condition = (temp_vars['evIndex'] <= 0)

        # Apply gen_ev_condition to only account for events of interest (generated solar events)
        for var in var_read_name_list:
            temp_vars[var] = temp_vars[var][gen_ev_condition] 

        temp_vars['n_init_evs'] = len(temp_vars['energy'])  # Save the number of the total simulated events

        # ============ Analysis Conditions ============
        print(' ------- Performing Analysis cuts ------- ')

        valid_condition = (temp_vars['scintFit'] & temp_vars['fitValid'])

        nhits_min = 20
        nhits_condition = (temp_vars['nhits'] >= nhits_min)
        
        energy_cut = 0 #MeV
        energy_condition = (temp_vars['energy'] >= energy_cut)

        posr_cut = 5700.0 # mm
        posr_condition = (temp_vars['posr_av'] <= posr_cut)

        itr_cut = 0.2
        itr_condition = (temp_vars['itr'] >= itr_cut)

        #mask_cut = 0xD82100000162C6
        #dcflag_condition = ((int(mask_cut) & temp_vars['dcFlagged']) == int(mask_cut))
        
        general_condition = valid_condition & nhits_condition & energy_condition & posr_condition & itr_condition

        # ========== Apply cut conditions and save the data ==========

        #subrunID_sv = np.append(subrunID_sv, subrunID_array[general_condition].astype(np.int16))
        #np.save(os.path.join(save_dir, 'subrunID'), subrunID_sv)
       
        # Apply Cuts        
        for var in var_save_name_list:
            data_filtered = np.extract(general_condition, temp_vars[var])
            data_dict[var] = np.append(data_dict[var], data_filtered)

        print(' ------- Analysis Cuts Applied -------')

        # Compute the corrected energy
        energy = data_dict['energy']
        posx = data_dict['posx']
        posy = data_dict['posy']
        posz = data_dict['posz_av']

        #data_dict['energy_corrected'] = ReconCalibrator(posx, posy, posz, energy)

        # Save the number of Initial Events
        data_dict['n_init_evs'] = np.append(data_dict['n_init_evs'], temp_vars['n_init_evs'])

        n_init_evs = np.sum(data_dict['n_init_evs']) # Value refreshed on each iteration as an acummulation of the Nº of events
        print(f'number of MC events (cumulative for all files!) {n_init_evs}')

        # ====== Save Data ======

        for var_sv in var_save_name_list:
            np.save(os.path.join(save_dir, var_sv), data_dict[var_sv])
        #np.save(os.path.join(save_dir, 'energy_corrected'), data_dict['energy_corrected'])

        # Save the generated number of events
        np.save(os.path.join(save_dir, 'n_init_evs'), np.array([n_init_evs]))


    return print('extraction concluded!')

'''
if __name__ == "__main__":

    read_dir = '/share/neutrino/snoplus/MonteCarlo/FullFill_bismsb_801/ScintFit_2p2ppo_2p2bismsB8_Solar_NueRun/'
    file_txt_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/2p2_ppo_bismsb/B8_solar_Nue/file_name_list/sublist_0.txt'
    save_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/2p2_ppo_bismsb/B8_solar_Nue/proof/'

    extract_data(read_dir, file_txt_dir, save_dir)
'''
