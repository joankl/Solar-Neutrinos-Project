'''
Main Code used to perform analysis and construction of observables for solar Nue analysis. This code must be complemented with the function SunDirection.py to obtain the
solar direction when analyzing real data. If the data is from MC then just extract the mc_direction.
Also the files to be readen must be defined. The code defines directories of the root files.
'''

#Function to Compute the radial position of events
def magnitude(vector): 
    x = vector[:,0]
    y = vector[:,1]
    z = vector[:,2]

    r = np.sqrt(x**2 + y**2 + z**2)
    r = r.astype(np.float32)
    return r

#Analysis code settings  ------------------------------
out_file_dir = 'analysis_out/en_up_2_pt_5_MeV/'

#Data Cuts:
energy_inf_cut = 2.5
energy_sup_cut = 12

time_res_inf_cut = -5.0
time_res_sup_cut = 7.0

posr_cut = 5000
# -------------------------------------------------------

for i_dx, input_file_dir in enumerate(flist):
    print(f'In file {input_file_dir}')
    
    #Load the Data:
    load_data = uproot.open(input_file_dir)
    
    #select the tree of event data and PMT info
    TTree_data_name = load_data.keys()[0]
    TTree_pmt_info_name = load_data.keys()[-1]
    
    event_data = load_data[TTree_data_name]
    pmt_data = load_data[TTree_pmt_info_name]
    
    #event info to be used:
    var_event_list = ['evtid', 'evtid_bi214', 'energy', 'position', 'hit_pmtid', 
                      'hit_residual', 'dc_flag', 'ev_time_day', 
                      'ev_time_sec', 'ev_time_nanosec', 'clockCount50']  #list the name of the varibles to be extracted and used for the solarnu analysis.
    
    #pmt info to be used
    var_pmt_list = ['pmt_id', 'pmt_pos_xyz', 'pmt_type']
    
    #Observables to save
    var_name_save_list = ['evtid', 'evtid_bi214', 'energy', 'posr', 'hit_residual', 'clockCount50']
    multi_cos_alpha = np.array([]) #create the empty list of the cos_alpha for the multiple PMTs record
    
    # Extract the variables with the name of the var_event_list in numpy.array from the .root file
    observables = {}
    
    # Cargar variables del árbol de eventos
    for var_name_i in var_event_list:
        observables[var_name_i] = np.array(event_data[var_name_i])
    
        # Transformar clockCount50
        if var_name_i == 'clockCount50':
            observables['clockCount50'] = observables['clockCount50'] * 20
    
    # Calcular posr
    observables['posr'] = magnitude(observables['position'])
    
    # Cargar variables del árbol de PMTs
    for var_name_i in var_pmt_list:
        observables[var_name_i] = np.array(pmt_data[var_name_i])
    
    # Filtrado de pmt_type válido
    pmt_type_condition = (observables['pmt_type'] == 1)
    pmt_id_valid = observables['pmt_id'][pmt_type_condition]
    
    # Condiciones
    hit_pmt_id_condition = np.in1d(observables['hit_pmtid'], pmt_id_valid)
    energy_condition = (observables['energy'] >= energy_inf_cut) & (observables['energy'] <= energy_sup_cut)
    time_res_condition = (observables['hit_residual'] >= time_res_inf_cut) & (observables['hit_residual'] <= time_res_sup_cut)
    posr_condition = (observables['posr'] <= posr_cut)
    
    general_condition = hit_pmt_id_condition & energy_condition & time_res_condition & posr_condition
    
    # Aplicar condición general a las variables del evento
    for var_name_i in var_event_list:
        observables[var_name_i] = observables[var_name_i][general_condition]
    observables['posr'] = observables['posr'][general_condition]

    print(f'selected energies form {input_file_dir}: {observables['energy']} with shape {observables['energy'].shape}')
    
    # Guardar variables deseadas
    for var_name_i in var_name_save_list:
        np.save(out_file_dir + var_name_i + f'_{i_dx}.npy', observables[var_name_i])
    
    N_samples = len(observables['hit_residual'])
    
    for sample_idx in range(N_samples):
        Utime_day = observables['ev_time_day'][sample_idx]
        Utime_sec = observables['ev_time_sec'][sample_idx]
        Utime_nsec = observables['ev_time_nanosec'][sample_idx]
    
        sun_dir = sun_direction(UTdays=Utime_day, UTsecs=Utime_sec, UTnsecs=Utime_nsec) * (-1)
    
        pmt_hit_id = observables['hit_pmtid'][sample_idx]
        pmt_hit_xyz = observables['pmt_pos_xyz'][pmt_hit_id]
    
        norm2 = np.linalg.norm(pmt_hit_xyz)
        pmt_hit_xyz = pmt_hit_xyz / norm2
    
        dot_prod = np.dot(sun_dir, pmt_hit_xyz)
        cos_alpha = dot_prod
    
        multi_cos_alpha = np.append(multi_cos_alpha, cos_alpha)
    
    # Guardar los cos(alpha)
    np.save(out_file_dir + 'cos_alpha' + f'_{i_dx}.npy', multi_cos_alpha)