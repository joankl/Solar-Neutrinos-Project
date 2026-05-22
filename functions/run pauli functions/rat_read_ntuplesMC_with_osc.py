'''                                                                                                                                                               
Function to read the ntuples.root of real data. It will extract the variables of 
interest and save it in separeted numpy arrays.
This script makes use of the RAT ROOT Tools to perform:

- Energy-Position dependent Corrections

Then, the script will only work for RAT above v8

The script needs:

- The ntuple.root file 
- PSelmaa data files

The script will save:
- The ntuple data of interest to analyze
- The Pee values for each saved ntuple event.

Edit Dates:
- 23/01/2026 Include RAT code to implement the energy correction
- 22/05/2026 Include ITR Cut above 0.2
'''

import numpy as np
from scipy.interpolate import interp1d
import uproot
import glob
import re
import os
import rat
import ROOT


def extract_data(read_data_dir, file_txt_dir, save_dir, pselmaa_data_dir):

    '''
    Function designed to read ntuples.root files. It will take a list of the ntuples.root
    file names and will save the analysis result in numpy files.

    Parameters:
    - read_data_dir: directory where the files with data are.
    - file_txt_dir:  directory where the files.txt with the ntuples file names are
    - save_dir: directory where the analysis files will be saved
    - pselmaa_data_dir: directory where the PSelmaa data is.
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
    var_read_name_list = ['evIndex', 'scintFit', 'fitValid', 'nhits',
                          'energy', 'posr_av', 'posx', 'posy', 'posz_av', 'parentKE1', 'itr']

    # ------- List of the variables to save -------
    var_save_name_list = ['energy', 'posr_av','posx', 'posy', 'posz_av', 'parentKE1']

    # Diccionary to save the acummulated  data
    data_dict = {var: np.array([]) for var in var_save_name_list  + ['n_init_evs']}

    # ============= Load PSelmaa Data =============

    PSelmaa_data = np.loadtxt(pselmaa_data_dir, skiprows=1)
    Pee = PSelmaa_data[:,1]
    Pee_energy = PSelmaa_data[:,0]

    # ============= Files reader =============
    file_name_list = read_files_txt(file_txt_dir)
    file_name_list.sort(key=orden_natural)

    print(f'file name list {file_name_list}')

    for i_dx, file_name_i in enumerate(file_name_list):
        full_dir_file_i = os.path.join(read_data_dir, file_name_i)
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

        n_init_evs_bfr_osc = len(temp_vars['energy'])  #N of events before rejection sampling for osc.
        print(f'Nº of initial events for this file before rejection sampling = {n_init_evs_bfr_osc}')

        # --------- Compute the Survival Probability ---------

        print(' ------ Computing Survival Probability ------ ')

        E_nu = temp_vars['parentKE1']  # Simulated neutrino energy

        # Interpolate Pee with the MC data energy
        Pee_f = interp1d(Pee_energy, Pee, kind='linear', bounds_error=False, fill_value=0)  #function that describes the interpolation
        Pee_int = Pee_f(E_nu) # Interpolation Values

        print(f'Interpolated Pee with values: \n {Pee_int}')

        # Random vector generation
        N_ev = len(E_nu) # Number of simulated event
        rng = np.random.default_rng(seed=42) # Seed for Reproducibility
        r_vec = rng.random(N_ev) # Random vector

        survival_mask = (r_vec < Pee_int)  # Binary mask. If False the event has 'oscillated' and will be removed

        # Event selection
        for var in var_read_name_list:
            temp_vars[var] = temp_vars[var][survival_mask] 
        #Pee_int = Pee_int[survival_mask] 

        print('------ Oscillated Data Obtained ------ ')

        # Number of initially simulated events after applying oscillations
        n_init_evs_aft_osc = len(temp_vars['energy']) #N of events after rejection sampling for osc.
        temp_vars['n_init_evs'] = n_init_evs_aft_osc

        print(f'Nº of initial events for this file after rejection sampling = {n_init_evs_aft_osc}')

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

        #Pee_int = Pee_int[general_condition]  # Apply cuts in Pee

        print(' ------- Analysis Cuts Applied -------')

        # Compute the corrected energy
        energy = data_dict['energy']
        posx = data_dict['posx']
        posy = data_dict['posy']
        posz = data_dict['posz_av']

        data_dict['energy_corrected'] = ReconCalibrator(posx, posy, posz, energy)

        #Save number of Initial Events
        n_init_evs = np.sum(data_dict['n_init_evs']) # Value refreshed on each iteration as an acummulation of the Nº of events


        # Save the number of Initial Events
        data_dict['n_init_evs'] = np.append(data_dict['n_init_evs'], temp_vars['n_init_evs'])

        # ====== Save Data ======

        for var_sv in var_save_name_list:
            np.save(os.path.join(save_dir, var_sv), data_dict[var_sv])
        np.save(os.path.join(save_dir, 'energy_corrected'), data_dict['energy_corrected'])

        # Save the generated number of events
        np.save(os.path.join(save_dir, 'n_init_evs'), np.array([n_init_evs]))

        # Save the Pee
        np.save(os.path.join(save_dir, 'Pee'), Pee_int)


    return print('extraction concluded!')

'''
if __name__ == "__main__":

    read_data_dir = '/share/neutrino/snoplus/MonteCarlo/FullFill_bismsb_801/ScintFit_2p2ppo_2p2bismsB8_Solar_NueRun/'
    file_txt_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/2p2_ppo_bismsb/B8_solar_Nue/file_name_list/sublist_0.txt'
    save_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/2p2_ppo_bismsb/B8_solar_Nue/proof/'
    pselmaa_data_dir = '/lstore/sno/joankl/solar_analysis/flux_prediction/PSelmaa/output_files/pselmaa_test_sun_pee_B16_GS98_b8.txt'
    extract_data(read_data_dir, file_txt_dir, save_dir, pselmaa_data_dir)
'''

