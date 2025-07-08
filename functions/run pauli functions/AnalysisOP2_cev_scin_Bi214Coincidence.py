# Imports that are relevant for the analysis

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

import ROOT
import rat
import glob
import os
import math
from array import array


'''
*Function Parameters:

# location of the input files
read_dir = "/share/neutrino/snoplus/Data/FullFill_2p2/ratds/"
file_txt_dir = ?
save_dir = '/lstore/sno/joankl/solar_analysis/real_data/data/'  
file_out_name = 'real_data_solar_analysis_i.root' #where the index i should be iterated for each sublist.txt!

'''

def extract_data(read_dir, file_txt_dir, save_dir, file_out_name):

  #Code settings parameters: --------------------------------------

  dcflag_cut = int(0x2100000042C2)  #Data cleanning Mask
  energy_inf_cut = 2.5    #Lower cut on the energy (MeV) of the events for the analysis
  Bi_214_max_energy = 3.5 #end-point energy of the Bi214 (MeV) (value from MC)
  Po_214_min_energy = 0.5 #int-point energy of the Po214 (MeV) (value from MC)
  Po_214_max_energy = 1.5 #end-point energy of the Po214 (MeV) (value from MC)
  dt_window = 1.0         #range of the time window to try to find a delayed Po214 for a prompt Bi214 in miliseconds
  dr_window = 1000.0      # upper limit for the distance between the prompt and the delayed reconstructed positions in milimeters

  #----------------------------------------------------------------

  def read_files_txt(file_txt_dir):
      """Lee una sublist.txt y devuelve una lista con solo los nombres de archivos en sublist.txt"""
      with open(file_txt_dir, "r", encoding="utf-8") as f:
          file_names = [line.strip() for line in f if line.strip()]
      return file_names

  def write_pmt_info(fout,reload=False):
    du =  rat.utility()
    if reload:
      du.LoadDBAndBeginRun()
    #du.LoadDBAndBeginRun()
    pmtinfo = du.GetPMTInfo()
    fout.cd()
    tpmt = ROOT.TTree("pmt","PMT information")

    pmt_id = array('i',[0])
    tpmt.Branch('pmt_id',pmt_id,'pmt_id/I')
      
    pmt_pos_xyz = array('d',3*[0.0])
    tpmt.Branch('pmt_pos_xyz',pmt_pos_xyz,'pmt_pos_xyz[3]/D')
    pmt_pos_sph = array('d',3*[0.0])
    tpmt.Branch('pmt_pos_sph',pmt_pos_sph,'pmt_pos_sph[3]/D')

    pmt_type = array('i',[0])
    tpmt.Branch('pmt_type',pmt_type,'pmt_type/I')

    for i_pmt in range(0,pmtinfo.GetCount()):
      pos = pmtinfo.GetPosition(i_pmt)
      pmt_pos_xyz[0] = pos.x()
      pmt_pos_xyz[1] = pos.y()
      pmt_pos_xyz[2] = pos.z()
      pmt_pos_sph[0] = pos.Theta()
      pmt_pos_sph[1] = pos.Phi()
      pmt_pos_sph[2] = pos.Mag()
      pmt_id[0] = i_pmt
      pmt_type[0] = pmtinfo.GetType(i_pmt)
      tpmt.Fill()
    #my_pmt_db[i_pmt] = {'type':pmtinfo.GetType(ipmt), 'position_xyz': [pos.x(),pos.y(),pos.z()], 'position_sph':[pos.Theta(),pos.Phi(),pos.Mag()]}
    tpmt.Write()

  #do_analysis()
  # Total number of processed events
  events = 0
  fitter = 0
  fout = ROOT.TFile.Open(save_dir + file_out_name,"RECREATE")
  print('type: ',type(fout))
  fout.cd()
  print('fout.cd() :', fout.cd())

  #TTree and Branches Creation ----------------------------------------------------

  tree = ROOT.TTree("T","output summary")
  #mcID = array('i',[0])         #To save montecarlo ID
  #tree.Branch('mcID', mcID, 'mcID/I')

  evt_en = array('d', [0.0])    #To save reconst. energy of event
  tree.Branch('energy', evt_en, 'energy/D')

  #evt_pos_en_err = array('d', [0.0])    #To save reconst. energy postive error of event
  #tree.Branch('positiveEnergyError', evt_en, 'positiveEnergyError/D')

  #evt_neg_en_err = array('d', [0.0])    #To save reconst. energy postive error of event
  #tree.Branch('negativeEnergyError', evt_en, 'negativeEnergyError/D')

  evtid = array('i',[0])
  tree.Branch('evtid',evtid,'evtid/I')

  #Brach which contains the GTID of the suspected Bi214
  evtid_bi214 = array('i', [0])
  tree.Branch('evtid_bi214',evtid_bi214,'evtid_bi214/I')
 
  #Save event Universal time
  ev_time_day = array('d', [0.0]) 
  tree.Branch('ev_time_day',ev_time_day, 'ev_time_day/D')
  ev_time_sec = array('d', [0.0])
  tree.Branch('ev_time_sec',ev_time_sec, 'ev_time_sec/D')
  ev_time_nanosec = array('d', [0.0])
  tree.Branch('ev_time_nanosec',ev_time_nanosec, 'ev_time_nanosec/D')

  #mc_pos = array('d',3*[0.0])
  #tree.Branch('mc_position',mc_pos,'mc_position[3]/D')
  #mc_mom = array('d',3*[0.0])
  #tree.Branch('mc_momentum',mc_mom,'mc_momentum[3]/D')

  evt_pos = array('d',3*[0.0])
  tree.Branch('position',evt_pos,'position[3]/D')

  evt_pos_Pos_err = array('d',3*[0.0])
  tree.Branch('pos_position_error',evt_pos,'pos_position_error[3]/D')

  evt_neg_Pos_err = array('d',3*[0.0])
  tree.Branch('neg_position_error',evt_pos,'neg_position_error[3]/D')

  #cleanning info: dcFlag
  dc_flag = array('l', [0])
  tree.Branch('dc_flag', dc_flag, 'dc_flag/L')

  #hit_info
  neck_nhit_pmt = array('i', [0])
  tree.Branch('neck_nhit_pmt', neck_nhit_pmt, 'neck_nhit_pmt/I')

  hit_pmtid = array('i',[0])
  tree.Branch('hit_pmtid',hit_pmtid,'hit_pmtid/I')

  hit_pmttime = array('d',[0.])
  tree.Branch('hit_pmttime',hit_pmttime,'hit_pmttime/D')

  #clockCount50:
  clockCount50 = array('l',[0])
  tree.Branch('clockCount50',clockCount50,'clockCount50/L')

  hit_residual = array('d',[0.])
  tree.Branch('hit_residual',hit_residual,'hit_residual/D')

  # Reading the Files -------------------------------------------------

  flist = read_files_txt(file_txt_dir)
 # print(type(flist))
  assert len(flist) != 0, f"Unexpected number of input files. Got {0}".format(len(flist))


  # Bringing RAT Tools:
  util = rat.utility()
  # time residual calculator
  timeResCalc = util.GetTimeResidualCalculator()
  stop = False

  for input_file_i in flist:
    print('file:', input_file_i)
    ev_counter = 0
    full_dir_flist = os.path.join(read_dir, input_file_i)  #Construct the full directory of the input files
    reader = ROOT.RAT.DU.DSReader(full_dir_flist)
   # print('Number OF FILES{}'.format(reader.GetEntryCount()))
    for ievent in range(0,reader.GetEntryCount()):
      rDS = reader.GetEntry(ievent)

      if stop:
        break
 
      for iev in range(0, rDS.GetEVCount()):
        events = events + 1
        
        ev = None
        if (events == ev):
          stop = True
          break
        
        #Grab the EV Branch
        rEV = rDS.GetEV(iev)

        #Get the flag data cleaning branch
        dataCleanFlag = rEV.GetDataCleaningFlags()
        flag_index = dataCleanFlag.GetLatestPass()
        #print('flag:', dataCleanFlag.GetFlags(flag_index).GetULong64_t(0))

        #DC Flag Cut Condition: 
        dcflag_value = dataCleanFlag.GetFlags(flag_index).GetULong64_t(0)
        dcflag_condition = ((dcflag_cut & dcflag_value) == dcflag_cut)

        if not (dcflag_condition):
          continue

        #Get the event time Brach
        rTime = rEV.GetUniversalTime()

        #Get the eventID
        evtid[0] = rEV.GetGTID()
        #print('evID', evtid[0])

        #Get the uncalibrated PMT branch
        uncalPMTs = rEV.GetUncalPMTs()
        #print('neckHits:', uncalPMTs.GetNeckCount())

        '''
        There is a subtlety here that was not considered yet
        A single event can actually lead to multiple vertices
        For example, pile-up could give that. Or very long events like muons that
        could spill over into the next time window
        For now assume just the first vertex
        '''
        # Check if this ev has a fit result
        n_vtx = 0
        if not rEV.FitResultExists("scintFitter"):
          continue
        fResult = rEV.GetFitResult("scintFitter")
        #print('scintf',fResult.GetValid())
        #print('fitvalid',rEV.FitResultExists("fitValid"))
        #f_valid = rEV.FitResultExists("fitValid")
        #f_scint = rEV.FitResultExists("scintFitter")
        
        if not fResult.GetValid():
          continue
   
        # Check that there is at least one reconstructed vertex. If there are nor valid reconstructed quantities, step to the next event
        if fResult.GetVertexCount() < 1:
          continue
        for ivtx in range(0,fResult.GetVertexCount()):

          fVertex = fResult.GetVertex(ivtx)

          if not(fVertex.ContainsPosition() and fVertex.ContainsEnergy() and fVertex.ValidPosition() and fVertex.ValidEnergy()):
            continue
          else:
            n_vtx = n_vtx + 1

        # Grab a few quantities of interest:
        # Print bool quantities
          #print('contains position?',fVertex.ContainsPosition())
          #print('contains energy?',fVertex.ContainsEnergy())
          #print('ValidPos?',fVertex.ValidPosition())
          #print('ValidEn?',fVertex.ValidEnergy())

          # Grab Reconstructed position and errors
          fPosition = fResult.GetVertex(ivtx).GetPosition()
          fPosPositionErr = fVertex.GetPositivePositionError()
          fNegPositionErr = fVertex.GetNegativePositionError()

          # Grab Reconstructed Energy 
          fEnergy = fResult.GetVertex(ivtx).GetEnergy()

          # Grab Reconstructed time
          fVertexTime =  fVertex.GetTime()

          #Energy condition to enter to the analysis coincidence due to a suspected Bi214:
          if (fEnergy >= energy_inf_cut) and (fEnergy <= Bi_214_max_energy):

            print('Performing BiPo214 Coincidence Analysis ')

            prompt_iev = iev
            print(f'prompt ev index = {prompt_iev}')

            rEV_prompt = rDS.GetEV(prompt_iev)
            
            #extract the quantities for the prompt
            t1 = (rEV_prompt.GetClockCount50()*20)*10**(-6)  #Transform from ns to ms
            x1 = fPosition.x()
            y1 = fPosition.y()
            z1 = fPosition.z()

            #initialize the quantities of the delayed as zeroes
            dt = 0
            dr = 0

            delay_iev = iev + 1
            while dt >= 0 and dt <= dt_window:

              print('In while loop dt window')
              print(f'delay ev index = {delay_iev}')

              rDS_delay = reader.GetEntry(delay_iev)

              #if delay_iev < rDS_delay.GetEVCount():
              #    rEV_delay = rDS_delay.GetEV(delay_iev)
              #else:
              #    print(f"Índice inválido: delay_iev={delay_iev}, EVCount={rDS_delay.GetEVCount()}")
              #    continue  # o return, o lo que tenga sentido en tu flujo

              for iev_delay in range(0, rDS_delay.GetEVCount()):

                rEV_delay = rDS_delay.GetEV(iev_delay)
                print(f'GTID of Prompt {rEV_prompt.GetGTID()}')
                print(f'GTID of Delay {rEV_delay.GetGTID()}')

                # Ensure the validity of the Delay event:
                # Correct Data-Cleanning Flag
                dataCleanFlag = rEV_delay.GetDataCleaningFlags()
                flag_index = dataCleanFlag.GetLatestPass()
                dcflag_value = dataCleanFlag.GetFlags(flag_index).GetULong64_t(0)
                dcflag_condition = ((dcflag_cut & dcflag_value) == dcflag_cut)

                #validity of existing quantities and scintFitter result on the delay event:
                print('testing the validity of the data ...')
                if not (dcflag_condition):
                  delay_iev += 1
                  continue

                n_vtx = 0
                if not(rEV_delay.FitResultExists("scintFitter")):
                  delay_iev += 1
                  continue

                fResult_delay = rEV_delay.GetFitResult("scintFitter")

                if not fResult_delay.GetValid():
                  delay_iev += 1
                  continue

                # Check that there is at least one reconstructed vertex
                if fResult_delay.GetVertexCount() < 1:
                  delay_iev += 1
                  continue

                for ivtx in range(0,fResult_delay.GetVertexCount()):
                  fVertex_delay = fResult_delay.GetVertex(ivtx)

                  if not(fVertex_delay.ContainsPosition() and fVertex_delay.ContainsEnergy() and fVertex_delay.ValidPosition() and fVertex_delay.ValidEnergy()):
                    delay_iev += 1
                    continue
                  else:
                    n_vtx = n_vtx + 1

                  print('Validity passed!')

                  #Now that all validity conditions are fulfilled, compute the observables of interest for the delay
                  fEnergy_delay = fResult_delay.GetVertex(ivtx).GetEnergy()
                  fPosition_delay = fResult_delay.GetVertex(ivtx).GetPosition()

                  #Cut on the energy of the delay event
                  print(f'Energy of the delay {fEnergy_delay} (MeV)')
                  if (fEnergy_delay < Po_214_min_energy) or (fEnergy_delay > Po_214_max_energy):
                    print('Delay energy out of range.')
                    delay_iev += 1
                    continue

                  print('performing computation of dt_delay and dr_delay')
                  try:
                    t2 = (rEV_delay.GetClockCount50()*20)*10**(-6)
                    print(f't1 = {t1}')
                    print(f't2 = {t2}')

                    x2 = fPosition_delay.x()
                    y2 = fPosition_delay.y()
                    z2 = fPosition_delay.z()

                    #Re-evaluate dt and dr
                    dt = t2 - t1

                    dx = x2 - x1
                    dy = y2 - y1
                    dz = z2 - z1 

                    dr = math.sqrt(dx**2 + dy**2 + dz**2)

                  except IndexError:
                    continue

                  print(f'dt = {dt} (ms) and dr = {dr} (mm)')
                  #Condition to save the suspected GTID of the prompt and get out of while
                  if (dt > 0 and dt <= dt_window) and (dr >= 0 and dr <= dr_window):
                    print(f'Pair of Coincidences found with dt = {dt} (ms) and dr = {dr} (mm)')
                    evtid_bi214[0] = rEV_prompt.GetGTID()
                    break

                  else:
                    print(f'Incompatible Delay with {dt} (ms) and {dr} (mm). Going for the next one.')
                    delay_iev += 1

          #print('Out of BiPo214 Coincidence Analysis.')

          #Energy Cut Condition to save all the events: Reject a lot of Background!
          if fEnergy < energy_inf_cut:
            continue

          evt_pos[0] = fPosition.x()
          evt_pos[1] = fPosition.y()
          evt_pos[2] = fPosition.z()

          print('ev_position:', evt_pos)

          evt_pos_Pos_err[0] = fPosPositionErr.x()
          evt_pos_Pos_err[1] = fPosPositionErr.y()
          evt_pos_Pos_err[2] = fPosPositionErr.z()
          #print('evt_pos_Pos_err:', evt_pos_Pos_err)

          evt_neg_Pos_err[0] = fNegPositionErr.x()
          evt_neg_Pos_err[1] = fNegPositionErr.y()
          evt_neg_Pos_err[2] = fNegPositionErr.z()
          #print('evt_neg_Pos_err:', evt_neg_Pos_err)

          #print('file:', input_file_i)
          print('event Energy', fEnergy)
          print('GTID', evtid[0])
       
          evt_en[0] = fEnergy
          #evt_pos_en_err[0] = fVertex.GetPositiveEnergyError()
          #evt_neg_en_err[0] = fVertex.GetNegativeEnergyError()

          #print('evt_pos_en_err ', fVertex.GetPositiveEnergyError())
          #print('evt_neg_en_err ', fVertex.GetPositiveEnergyError())


          #Reconstructed Time
          ev_time_day[0] = rTime.GetDays()
          ev_time_sec[0] = rTime.GetSeconds()
          ev_time_nanosec[0] = rTime.GetNanoSeconds()
          #print('ev_time',rTime.GetDays())

          #Data Cleanning dcFlag
          dc_flag[0] = dataCleanFlag.GetFlags(flag_index).GetULong64_t(0)
          print('dcFlagged:', dc_flag[0])

          #Number of neck hits        
          neck_nhit_pmt[0] = uncalPMTs.GetNeckCount()

          #clockCount50
          print(f'clockCount50 of event {rEV.GetClockCount50()}')
          clockCount50[0] = rEV.GetClockCount50()

          ev_counter+= 1

          print('ev_counter: ', ev_counter)

          #print('neckHit:', necknhit_pmt[0])
        # normalized position vector
        #R = fPosition.Unit()
        #if fVertex.ValidDirection():
          #fDirection = fResult.GetVertex(ivtx).GetDirection()
          #evt_mom[0] = fDirection.x()
          #evt_mom[1] = fDirection.y()
          #evt_mom[2] = fDirection.z()
        #else:
          #evt_mom[0] = 0.
          #evt_mom[1] = 0.
          #evt_mom[2] = 0.
          
        #fitter = fitter + 1
        
          
        # now loop over the hits
          calibratedPMTs = rEV.GetCalPMTs() 
          for iPMT in range(0, calibratedPMTs.GetAllCount()):
            #print('PRINT HERE', iPMT)
            pmtCal = calibratedPMTs.GetAllPMT( iPMT )
            pmttime =  pmtCal.GetTime()
            pmtid = pmtCal.GetID()
            #fp = ROOT.RAT.DU.Point3D(0,fPosition.x(),fPosition.y(),fPosition.z())
            #residimu_fileal = timeResCalc.CalcTimeResidual( pmtid, pmttime, fp,fVertex.GetTime())
            residual = timeResCalc.CalcTimeResidual( pmtid, pmttime,fPosition, fVertexTime)          

            hit_pmtid[0] = pmtid
            hit_pmttime[0] = pmttime
            hit_residual[0] = residual
          
          # variable to hold the PMT hit type
          #htype = 0
          
          # PMT specific information
          #fPMTPosition = pmtInfo.GetPosition( pmtid )
          #PMTx = fPMTPosition.x();
          #PMTy = fPMTPosition.y();
          #PMTz = fPMTPosition.z();
          #PMTr = math.sqrt( PMTx*PMTx + PMTy*PMTy + PMTz*PMTz );
          
          # Now check the MC branch for the history of this specific PMT hit
         
          #for imcpmt in range(0,rMC.GetMCPMTCount()):
            #print('PRINT HERE', impcpmt)
            #rMCPMT = rMC.GetMCPMT(imcpmt)
            #if rMCPMT.GetID() == pmtid:
              # Found the right PMT
              # Of course, now we could actually have multiple PEs
              
              #Usually one would need to match through the track ID
              
             # rMCPE = rMCPMT.GetMCPE(0)
              #for i_mcpe in range(0,rMCPMT.GetMCPECount()):
              #  rMCPE = rMCPMT.GetMCPE(i_mcpe)
              #  # Now try to compare the track IDs
              #  #if rMCPE.GetPhotonTrackID() == 
              #if rMCPE.GetFromHistory(ROOT.RAT.DS.MCPE.hCherenkov):
                # this is a hit from Cherenkov
               # hit_type[0] = ROOT.RAT.DS.MCPE.hCherenkov  
              #elif rMCPE.GetFromHistory(ROOT.RAT.DS.MCPE.hScintillation):
               # hit_type[0] = ROOT.RAT.DS.MCPE.hScintillation
              #else:
               # continue
              #hit_type[0] = hit_type
              # Get out of the MCPMT loop
              #break
          # end of MCPMT loop
          #tree.Print()
          #print(f"data : {evtid} {hit_pmtid} {hit_residual} {hit_pmttime} {hit_type}")
            tree.Fill()
          # Fill in some of the histograms
        # end of calpmt loop
      # end of GetEVCount loop
  tree.Write()

  write_pmt_info(fout)

  fout.Close()

if __name__ == "__main__":
  read_dir = '/share/neutrino/snoplus/Data/FullFill_2p2/ratds/'
  file_txt_dir = '/lstore/sno/joankl/solar_analysis/real_data/file_name_list/sublist_0.txt'
  save_dir = '/lstore/sno/joankl/solar_analysis/real_data/proof/'
  file_out_name = 'xReal_BiPoCoincidence.root'
  extract_data(read_dir, file_txt_dir, save_dir, file_out_name)
