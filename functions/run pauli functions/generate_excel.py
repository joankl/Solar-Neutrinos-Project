'''
Python module to generate an unified excel sheet from
a set of observables in numpy format. You may write
in the observable list the correct name of the numpy
files to be readen.
'''

import numpy as np 
import pandas as pd

# ======= Observable List =======

obs_list = ['eventID', 'runID']

# ===============================


obs_dict = {obs_i: [] for obs_i in obs_list} # Create empty dictionary to save the observables.

type = 'Analysis20_bMR'                      # Type of document according ot it real data generation name. To be included in the excel file name


for obs_i in obs_dict:

	obs = np.load(obs_i + '.npy')             #Load the numpy file
	obs = obs.tolist()

	obs_dict[obs_i].append(obs)


# Transform to pandas dataframe
df = pd.DataFrame(obs_list)

# Save as xlsx file
df.to_excel(f"filtered_solar_{type}.xlsx", index=False)
