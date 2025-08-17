'''
Analysis code to read the .root files and extract the variables of interest and perform calculations on the angle between
a hit PMT and the Sun direction. It will return the separated np.arrays of the observables energy, posr, time_res and cos_alpha for each event record
'''

import numpy as np
import uproot

#Function to Compute the radial position of events
def magnitude(vector): 
    x = vector[:,0]
    y = vector[:,1]
    z = vector[:,2]

    r = np.sqrt(x**2 + y**2 + z**2)
    r = r.astype(np.float32)
    return r

#Analysis code settings  ------------------------------
#Read files:
input_file_dir = '' #directory and name to read the input file
out_file_dir = '' #directory to save the output files

#Data Cuts:
energy_inf_cut = 2.5
energy_sup_cut = 12

time_res_inf_cut = -5.0
time_res_sup_cut = 7.0

posr_cut = 5500

dcflag_cut = int(0x2100000042C2)

#Nº of splitten loops to analyze the data: 
#N_split = 20 # the more split, the less memory expensive is the code
#------------------------------------------------------

#Load the Data:
load_data = uproot.open(input_file_dir)

print('trees: ', load_data.keys())

#select the tree of event data and PMT info
event_data = load_data[]
pmt_data = load_data[]

#event info to be used:
var_event_list = ['evtid', 'energy', 'position', 'hit_pmtid', 
                  'hit_residual', 'ev_time_day', 
                  'ev_time_sec', 'ev_time_nanosec', 'dc_flag']  #list the name of the varibles to be extracted and used for the solarnu analysis.

#pmt info to be used
var_pmt_list = ['pmt_id', 'pmt_pos_xyz', 'pmt_type']

#Observables to save
var_name_save_list = ['energy', 'posr', 'hit_residual']
multi_cos_alpha = np.array([]) #create the empty list of the cos_alpha for the multiple PMTs record

# Extract the variables with the name of the var_event_list in numpy.array from the .root file
for var_name_i in var_event_list:
    locals()[var_name_i] = np.array(event_data[var_name_i])
    print(f'Loaded Observable {var_name_i} and values: ', locals()[var_name_i])

# Compute the distance from the reconst. vertx to the center of the detector
posr = magnitude(position)

# Extract the pmt info
for var_name_i in var_pmt_list:
    locals()[var_name_i] = np.array(pmt_data[var_name_i])
    print(f'Loaded PMT infor {var_name_i} and values: ', locals()[var_name_i])

# Select the valid PMT information through pmt_type selection
pmt_type_condition = (pmt_type == 1)

#Selecting the pmt_id_valid is enough to filter all the PMT info (as the pmt_pos):
#for example, selecting the valid hit_pmtid will automatically choose the valid pmt_pos using the pmtid as the index of the data
pmt_id_valid = pmt_id[pmt_type_condition] #list with the valid pmt_id

#Extract the valid components of the event data through the selection of the pmt_id_valid:
hit_pmt_id_condition = np.in1d(hit_pmtid, pmt_id_valid)
#Also include  general cuts on energy and event position, and a dcFlag cut
energy_condition = (energy >= energy_inf_cut) & (energy <= energy_sup_cut)
time_res_condition = (hit_residual >= time_res_inf_cut) & (hit_residual <= time_res_sup_cut)
posr_condition = (posr <= 5500 )
dcFlag_condition = ((dcflag_cut & dc_flag) == int(mask_cut))

general_condition = hit_pmt_id_condition & energy_condition & time_res_condition & posr_condition & dcFlag_condition

#Loop to extract the samples which verify general condition
for var_name_i in var_event_list:
    #print('before pmt_id selection', locals()[var_name_i].shape)
    locals()[var_name_i] = locals()[var_name_i][general_condition]  #selection of the observable entries that verify the general condition
    #print('after pmt_id selection', locals()[var_name_i].shape)
posr = posr[general_condition]

print('UTDays: ', ev_time_day)

#Now, save the observables of interest from var_event_list after general condition cut
for var_name_i in var_name_save_list:
    np.save(out_file_dir + var_name_i, locals()[var_name_i])

#Loop over PMTs records to compute the angle between each hit PMT and the Sun direction:
N_samples = len(hit_residual) #Number of samples to be analyzed
#sample_idx_split = np.slipt(np.array(range(N_samples)), N_split) # A list splitten in N_sample arrays which contains the sample index

for sample_idx in range(N_samples):

    #Extract universal time for this event:
    Utime_day = ev_time_day[sample_idx]
    Utime_sec = ev_time_sec[sample_idx]
    Utime_nsec = ev_time_nanosec[sample_idx]

    #Compute the Sun's outward direction
    sun_dir = sun_direction(UTdays = Utime_day, UTsecs = Utime_sec, UTnsecs = Utime_nsec)*(-1)
    print(f'sun direction: {sun_dir}')

    #Extract the hit PMT coordinates:
    pmt_hit_id = hit_pmtid[sample_idx]  # first take the id of the hit PMT to look for its coordinates. 
    pmt_hit_xyz = pmt_pos_xyz[pmt_hit_id] # The ID is equivalent to the index of the pmt_pos.
    print(f'pmt direction: {pmt_hit_xyz}')

    #Compute the cos_alpha of the angles between hit_pmt with the sun direction by the scalar product definition: The Sun vector verifies norm = 1
    #cos(alpha) = dot_prod/norm1*norm2  
    norm2 = np.linalg.norm(pmt_hit_xyz)
    pmt_hit_xyz = pmt_hit_xyz/norm2  #normalized vector
    
    dot_prod = np.dot(sun_dir, pmt_hit_xyz)
    
    cos_alpha = dot_prod
    print(f'cos of angle = {cos_alpha}')
    
    multi_cos_alpha = np.append(multi_cos_alpha, cos_alpha)

    #save the cos_alpha computation
    np.save(out_file_dir + 'cos_alpha', multi_cos_alpha)