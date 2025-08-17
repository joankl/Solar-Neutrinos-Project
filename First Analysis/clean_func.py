def unique_ID(file, evID, subset):

    #function which extracts info. of an event given the event ID.
    #It function also clean the repeated eventID in PMT info using subsets division
    #It function return a Pandas Dataframe with the information

    # -Parameters:

    #file: file with information after passing through uproot
    #evID is the ID number
    #subset -> 2x1 array of index: [i_init_subset(n) , i_init_subset(n+1)]. exm: [0,1], [1,2], [2,3], ...
    #-------------------------------------------------------
    import numpy as np
    from numpy import array, where, shape, reshape
    import pandas as pd

    #Load file
    #EV info
    data1 = file['T;7']
    
    evtid = np.array(data1['evtid'])
    mc_position = np.array(data1['mc_position'])
    position = np.array(data1['position'])
    hit_pmtid = np.array(data1['hit_pmtid'])
    hit_pmttime = np.array(data1['hit_pmttime'])
    hit_residual = np.array(data1['hit_residual'])
    hit_type = np.array(data1['hit_type'])


    #Exclusion of repeated eventID of different simulations: Subset division

    index = where(evtid == evID)[0]
    break_data = np.array([index[0]], dtype = np.int64)  # #indices where data Breaks - Indices where subsets start

    for i in range(len(index)-1):
        if index[i+1] != index[i]+1:
            break_data = np.append(break_data,index[i+1])
    break_data = np.append(break_data,index[-1]+1)

    init =  break_data[subset[0]]
    final_ = break_data[subset[1]]

    #Extract new data: 

    new_evtid = evtid[init:final_]
    new_mc_position = mc_position[init:final_]
    new_position = position[init:final_]
    new_hit_pmtid = hit_pmtid[init:final_]
    new_hit_pmttime = hit_pmttime[init:final_]
    new_hit_residual = hit_residual[init:final_]
    new_hit_type = hit_type[init:final_]

    #Construct DataFrame:

    data = {'ID': new_evtid,
            'mc_position': new_mc_position.tolist(),
            'position': new_position.tolist(),
            'hit_pmtid': new_hit_pmtid,
            'hit_residual': new_hit_residual,
            'hit_type': new_hit_type,
            }

    data_frame = pd.DataFrame(data)

    return data_frame
#--------------------------------------------------------------------------------------------------------------------------

def EV_info(file, evID,):

    #function which extracts info. of an event given the event ID.
    #It function also clean the repeated eventID in PMT info using subsets division
    
    #import uproot
    import numpy as np
    from numpy import array, where, shape, reshape
    import pandas as pd

    #Load file
    #EV info
    data1 = file['T;7']
    
    evtid = np.array(data1['evtid'])
    mc_position = np.array(data1['mc_position'])
    position = np.array(data1['position'])
    hit_pmtid = np.array(data1['hit_pmtid'])
    hit_residual = np.array(data1['hit_residual'])
    hit_type = np.array(data1['hit_type'])
    
    #Extract the event given the evtid----------------------
    condition = (evID == evtid)

    #if condition satisfied, extract the information in the index where the condition is True
    mc_position_ev = np.extract(condition, mc_position)
    position_ev = np.extract(condition, position)
    residual_ev = np.extract(condition, hit_residual)
    type_e = np.extract(condition, hit_type)
    
    #Indices for which elements are True in eventID then extract the hit_pmtID and save in pmtif_ev
    pmtid_ev = np.extract(condition, hit_pmtid) 
    
    return pmtid_ev

#-------------------------------------------------------------------------------------------

def clean_pmt_type(file, type_, evID):
    
    import numpy as np
    from numpy import array, where, shape, reshape
    
    #upload file
    #file = uproot.open("/sno/py_out1.root")
    
    #pmt_info
    pmt_info = file['pmt;1']
    
    pmt_id = array(pmt_info['pmt_id'])
    # pmt_pos_xyz = array(pmt_info['pmt_pos_xyz'])
    # pmt_pos_sph = array(pmt_info['pmt_pos_sph'])
    pmt_type = array(pmt_info['pmt_type'])
    
    #selection of PMT type-------------------------------------
    clean_pmt_id = np.array([],np.int32)                   #pmtID of the desired type
    valid_index_pmt_type = np.array([],np.int32)           #index of the desired PMTs for the type required
    
    for type_i in type_:
        p = where(pmt_type == type_i)[0]  # p get the indices where the PMT_type are for each type_ desired
        valid_index_pmt_type = np.concatenate((valid_index_pmt_type,p))
        
    for valid_i in valid_index_pmt_type:
        clean_pmt_id = np.append(clean_pmt_id,pmt_id[valid_i])
        

    # Extract the indices and ID of the clean PMT given the event-------------------
    
    pmtid_ev = EV_info(file, evID)                                    #Get the array of pmtID which verify a hit
    
    pmt_clean_id_i = where(np.in1d(clean_pmt_id, pmtid_ev))[0]        # Get the index on the PMT info when the PMT cleaned list is in the PMT ID event
    
    return pmt_clean_id_i, clean_pmt_id

#-----------------------------------------------------------------

def GetCoordinateshitPMT(file, type_, evID):
    
    import numpy as np
    from numpy import array, where, shape, reshape

    #Function that receives the type of PMT and event ID of interest and return the PMT coordinates
    #ruturn PMT coordinates in [(x,y,z), (tetha,phi,r)]
    
    #Load file
    #file = uproot.open("/sno/py_out1.root")
    pmt_info = file['pmt;1']
    
    pmt_id = array(pmt_info['pmt_id'])
    pmt_pos_xyz = array(pmt_info['pmt_pos_xyz'])
    pmt_pos_sph = array(pmt_info['pmt_pos_sph'])
    
    hitpmt_coordxyz = []
    hitpmt_coordsphe = []
    
    #Look fot the clean PMT IDs from clean_pmt_type function
    
    pmt_clean_id_i ,clean_pmt_id = clean_pmt_type(file, type_, evID)
    
    for i in pmt_clean_id_i:
        hitpmt_coordxyz.append(pmt_pos_xyz[clean_pmt_id[i]])
        hitpmt_coordsphe.append(pmt_pos_sph[clean_pmt_id[i]])
        

        
    return array(hitpmt_coordxyz), array(hitpmt_coordsphe)