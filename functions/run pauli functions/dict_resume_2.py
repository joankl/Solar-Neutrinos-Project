'''
Python script (v2 of dict_resume) to read and save the resume of the
hotspot, atmospheric and IBD coincidence analysis.
'''

import glob
import os
import pickle
import numpy as np

main_dir = '/lstore/sno/joankl/solar_analysis/real_data/bisMSB/Analysis15_bMR/ntuple_output/bkg_candidates/output_*/'
save_dir = '/lstore/sno/joankl/solar_analysis/real_data/bisMSB/Analysis15_bMR/ntuple_output/resume_files/'

os.makedirs(save_dir, exist_ok=True)

dict_list = ['atm_dict', 'hs_dict'] 
dict_keys = ['counter', 'eventID', 'runID']


def flatten_to_list(item, target_list):
    """
    Agrega elementos a target_list de forma recursiva sin importar
    si son escalares, listas, tuplas o numpy arrays.
    """
    if isinstance(item, (list, tuple)):
        for sub in item:
            flatten_to_list(sub, target_list)
    elif isinstance(item, np.ndarray):
        for sub in item.ravel():
            target_list.append(sub.item() if hasattr(sub, 'item') else sub)
    else:
        # Es un escalar (int, float, np.int32, np.int64, etc.)
        target_list.append(item.item() if hasattr(item, 'item') else item)


for dict_i in dict_list:
    print(f'Reading dictionary {dict_i}.pkl')

    dict_to_save = {var_i: [] for var_i in dict_keys}
    full_dict_dir_list = glob.glob(main_dir + dict_i + '.pkl')

    for dir_i in full_dict_dir_list:
        print(f"Trying to read: {dir_i}")

        with open(dir_i, 'rb') as f:
            dict_file = pickle.load(f)

        for var_i in dict_keys:
            val = dict_file.get(var_i, [])
            flatten_to_list(val, dict_to_save[var_i])

        del dict_file

    out_path = os.path.join(save_dir, f'{dict_i}.pkl')
    with open(out_path, 'wb') as f:
        pickle.dump(dict_to_save, f, protocol=pickle.HIGHEST_PROTOCOL)

    del dict_to_save
    print(f'Dictionary {dict_i}.pkl saved successfully!')

print('Done! :)')