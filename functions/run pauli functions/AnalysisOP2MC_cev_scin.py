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
  print('In extract_data function')

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


#def do_analysis(ev = None):
  # Total number of processed events
  events = 0
  fitter = 0
  fout = ROOT.TFile.Open(save_dir + file_out_name,"RECREATE")
  fout.cd()
  tree = ROOT.TTree("T","output summary")
  mcID = array('i',[0])         #To save montecarlo ID
  evt_en = array('d', [0.0])    #To save energy of event
  evtid = array('i',[0])
  tree.Branch('evtid',evtid,'evtid/I')
  tree.Branch('mcID', mcID, 'mcID/I')
  tree.Branch('energy', evt_en, 'energy/D')
 
  mc_pos = array('d',3*[0.0])
  tree.Branch('mc_position',mc_pos,'mc_position[3]/D')
  mc_mom = array('d',3*[0.0])
  tree.Branch('mc_momentum',mc_mom,'mc_momentum[3]/D')

  evt_pos = array('d',3*[0.0])
  tree.Branch('position',evt_pos,'position[3]/D')
  evt_mom = array('d',3*[0.0])
  tree.Branch('momentum',evt_mom,'momentum[3]/D')
  #evt_mom = array('f',3*[0.0])
  #tree.Branch('momentum',evt_mom,'momentum[3]/D')

  hit_pmtid = array('i',[0])
  tree.Branch('hit_pmtid',hit_pmtid,'hit_pmtid/I')
  hit_pmttime = array('d',[0.])
  tree.Branch('hit_pmttime',hit_pmttime,'hit_pmttime/D')
  hit_residual = array('d',[0.])
  tree.Branch('hit_residual',hit_residual,'hit_residual/D')
  hit_type = array('i',[0])
  tree.Branch('hit_type',hit_type,'hit_type/I')

  #clockCount50:
  clockCount50 = array('l',[0])
  tree.Branch('clockCount50',clockCount50,'clockCount50/L')
  
  flist = read_files_txt(file_txt_dir)
 # print(type(flist))
  assert len(flist) != 0, f"Unexpected number of input files. Got {0}".format(len(flist))
 # assert len(flist) != 0, f"Unexpected number of input files. Got {len(flist)}"

  util = rat.utility()
  # time residual calculator
  timeResCalc = util.GetTimeResidualCalculator()
  #pmtInfo = util.GetPMTInfo()
  stop = False
 # print('NUMBER OF FILES{}'.format(len(flist)))
  for input_file_i in flist:
    full_dir_flist = os.path.join(read_dir, input_file_i)  #Construct the full directory of the input files
    reader = ROOT.RAT.DU.DSReader(full_dir_flist)
   # print('Number OF FILES{}'.format(reader.GetEntryCount()))
    for ievent in range(0,reader.GetEntryCount()):
      rDS = reader.GetEntry(ievent)
  
      #for rDS, _ in  rat.dsreader(flist):
      
      if stop:
        break
      # Grab the MC branch firsfirst
      rMC = rDS.GetMC()
      
      mcID[0] = rMC.GetMCID()
      print(mcID[0],'mcID')
      mc_p = rMC.GetMCParticle(0).GetPosition()
      mc_pos[0] = mc_p.x()
      mc_pos[1] = mc_p.y()
      mc_pos[2] = mc_p.z()
      
      mc_m = rMC.GetMCParticle(0).GetMomentum()
      mc_mom[0] = mc_m.x()
      mc_mom[1] = mc_m.y()
      mc_mom[2] = mc_m.z()

    
      for iev in range(0, rDS.GetEVCount()):
        events = events + 1
       # print("Event {}".format(events))

        ev = None
        if (events == ev):
          stop = True
          break
        
        rEV = rDS.GetEV(iev)
        evtid[0] = rEV.GetGTID()
        
        
        '''
        There is a subtlety here that was not considered yet
        A single event can actually lead to multiple vertices
        For example, pile-up could give that. Or very long events like muons that
        could spill over into the next time window
        For now assume just the first vertex
        '''
        # Check if this ev has a fit result
        n_vtx = 0
        if (not rEV.FitResultExists("scintFitter")):
          continue
        fResult = rEV.GetFitResult("scintFitter")
        if not fResult.GetValid():
          continue
        # Check that there is at least one reconstructed vertex
        if fResult.GetVertexCount() < 1:
          continue
        for ivtx in range(0,fResult.GetVertexCount()):
          fVertex = fResult.GetVertex(ivtx)
          evt_en[0] = fVertex.GetEnergy()
          print('event Energy', evt_en[0])

          #clockCount50
          clockCount50[0] = rEV.GetClockCount50()
          print('clockCount50', clockCount50[0])

          if not (fVertex.ContainsPosition() and fVertex.ValidPosition() and fVertex.ValidEnergy()):
            continue
          else:
            n_vtx = n_vtx + 1
          # Grab a few quantities of interest
          # Reconstructed position
        fPosition = fResult.GetVertex(ivtx).GetPosition()
        evt_pos[0] = fPosition.x()
        evt_pos[1] = fPosition.y()
        evt_pos[2] = fPosition.z()
        
        # normalized position vector
        #R = fPosition.Unit()
        if fVertex.ValidDirection():
          fDirection = fResult.GetVertex(ivtx).GetDirection()
          evt_mom[0] = fDirection.x()
          evt_mom[1] = fDirection.y()
          evt_mom[2] = fDirection.z()
        else:
          evt_mom[0] = 0.
          evt_mom[1] = 0.
          evt_mom[2] = 0.
          
        fitter = fitter + 1
        
          
        # now loop over the hits
        calibratedPMTs = rEV.GetCalPMTs() 
        for iPMT in range(0, calibratedPMTs.GetAllCount()):
         # print('PRINT HERE', iPMT)
          pmtCal = calibratedPMTs.GetAllPMT( iPMT )
          pmttime =  pmtCal.GetTime()
          pmtid = pmtCal.GetID()
          #fp = ROOT.RAT.DU.Point3D(0,fPosition.x(),fPosition.y(),fPosition.z())
          #residual = timeResCalc.CalcTimeResidual( pmtid, pmttime, fp,fVertex.GetTime())
          residual = timeResCalc.CalcTimeResidual( pmtid, pmttime,fPosition,fVertex.GetTime())          

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
          
          for imcpmt in range(0,rMC.GetMCPMTCount()):
            #print('PRINT HERE', impcpmt)
            rMCPMT = rMC.GetMCPMT(imcpmt)
            if rMCPMT.GetID() == pmtid:
              # Found the right PMT
              # Of course, now we could actually have multiple PEs
              '''
              Usually one would need to match through the track ID
              '''
              rMCPE = rMCPMT.GetMCPE(0)
              #for i_mcpe in range(0,rMCPMT.GetMCPECount()):
              #  rMCPE = rMCPMT.GetMCPE(i_mcpe)
              #  # Now try to compare the track IDs
              #  #if rMCPE.GetPhotonTrackID() == 
              if rMCPE.GetFromHistory(ROOT.RAT.DS.MCPE.hCherenkov):
                # this is a hit from Cherenkov
                hit_type[0] = ROOT.RAT.DS.MCPE.hCherenkov  
              elif rMCPE.GetFromHistory(ROOT.RAT.DS.MCPE.hScintillation):
                hit_type[0] = ROOT.RAT.DS.MCPE.hScintillation
              else:
                continue
              #hit_type[0] = hit_type
              # Get out of the MCPMT loop
              break
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
  read_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/solar_Nue/Nue_simulation_data/'
  file_txt_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/solar_Nue/file_name_list/sublist_0.txt'
  save_dir = '/lstore/sno/joankl/solar_analysis/mc_data/main_simulations/solar_Nue/proof/'
  file_out_name = 'xMC.root'
  extract_data(read_dir, file_txt_dir, save_dir, file_out_name)
