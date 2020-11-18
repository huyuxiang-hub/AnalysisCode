
#include "NormalAnaMgr.hh"
//  for event
#include <sstream>
#include <cassert>
#include "junoHit_PMT.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4EventManager.hh"
#include "G4TrackingManager.hh"
#include "G4OpticalPhoton.hh"
#include "G4PhysicalVolumeStore.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/SniperDataPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperLog.h"
#include "RootWriter/RootWriter.h"

#include "NormalTrackInfo.hh"
#include "EvtNavigator/NavBuffer.h"
#include "BufferMemMgr/IDataMemMgr.h"
#include "Event/SimHeader.h"

DECLARE_TOOL(NormalAnaMgr);

NormalAnaMgr::NormalAnaMgr(const std::string& name) 
    : ToolBase(name)
{
    declProp("EnableNtuple", m_flag_ntuple=true);
    m_evt_tree = 0;
    m_step_no = 0;
    declProp("EnableHitInfo",m_flag_hitinfo=false);


}

NormalAnaMgr::~NormalAnaMgr()
{

}

void
NormalAnaMgr::BeginOfRunAction(const G4Run* /*aRun*/) {
    if (not m_flag_ntuple) {
        return;
    }
    // check the RootWriter is Valid.

    SniperPtr<RootWriter> svc(*getParent(), "RootWriter");
    if (svc.invalid()) {
        LogError << "Can't Locate RootWriter. If you want to use it, please "
                 << "enalbe it in your job option file."
                 << std::endl;
        return;
    }
    m_evt_tree = svc->bookTree("SIMEVT/evt", "evt");
    m_evt_tree->Branch("evtID", &m_eventID, "evtID/I");
    m_evt_tree->Branch("edep", &m_edep, "edep/F");
    m_evt_tree->Branch("edepX", &m_edep_x, "edepX/F");
    m_evt_tree->Branch("edepY", &m_edep_y, "edepY/F");
    m_evt_tree->Branch("edepZ", &m_edep_z, "edepZ/F");
    m_evt_tree->Branch("nPhotons", &m_nPhotons, "nPhotons/I");
    m_evt_tree->Branch("totalPE", &m_totalPE, "totalPE/I");     

 
   if(m_flag_hitinfo==true)
   {
    m_evt_tree->Branch("nPE", &m_nPE);
    m_evt_tree->Branch("energy", &m_energy);
    m_evt_tree->Branch("hitTime",&m_hitTime);
    m_evt_tree->Branch("pmtID", &m_pmtID);
    m_evt_tree->Branch("PETrackID", &m_peTrackID);

    m_evt_tree->Branch("isCerenkov", &m_isCerenkov);
    m_evt_tree->Branch("isReemission", &m_isReemission);
    m_evt_tree->Branch("isOriginalOP", &m_isOriginalOP);
    m_evt_tree->Branch("OriginalOPTime", &m_OriginalOPTime);

    // PMT
    m_evt_tree->Branch("nPMTs", &m_npmts_byPMT, "nPMTs/I");
    m_evt_tree->Branch("nPE_byPMT", &m_nPE_byPMT);
    m_evt_tree->Branch("PMTID_byPMT",&m_PMTID_byPMT);
    // - 2015.10.10 Tao Lin <lintao@ihep.ac.cn>
    //   Hit's position
    m_evt_tree->Branch("LocalPosX",&m_localpos_x);
    m_evt_tree->Branch("LocalPosY",&m_localpos_y);
    m_evt_tree->Branch("LocalPosZ",&m_localpos_z);
    // - 2016.04.17 Tao Lin <lintao@ihep.ac.cn>
    //   Hit's direction
    m_evt_tree->Branch("LocalDirX",&m_localdir_x);
    m_evt_tree->Branch("LocalDirY",&m_localdir_y);
    m_evt_tree->Branch("LocalDirZ",&m_localdir_z);

    // - 2017.03.01 Tao Lin <lintao@ihep.ac.cn>
    //   Hit's Global Position
    m_evt_tree->Branch("GlobalPosX",&m_globalpos_x);
    m_evt_tree->Branch("GlobalPosY",&m_globalpos_y);
    m_evt_tree->Branch("GlobalPosZ",&m_globalpos_z);

    m_evt_tree->Branch("BoundaryPosX", &m_boundarypos_x);
    m_evt_tree->Branch("BoundaryPosY", &m_boundarypos_y);
    m_evt_tree->Branch("BoundaryPosZ", &m_boundarypos_z);
   

    m_step_no = new TH1I("stepno", "step number of optical photons", 1000, 0, 1000);
    svc->attach("SIMEVT", m_step_no);
    } 

}


void
NormalAnaMgr::EndOfRunAction(const G4Run* /*aRun*/) {

}

void
NormalAnaMgr::BeginOfEventAction(const G4Event* evt) {
    // initialize the evt tree
 
   m_eventID = evt->GetEventID();
    m_edep = 0.;
    m_edep_x = 0.;
    m_edep_y = 0.;
    m_edep_z = 0.;
    m_nPhotons = 0;
    m_totalPE = 0;
 if(m_flag_hitinfo==true)
   {
    m_nPE            .clear()                ;
    m_energy         .clear()                ;
    m_hitTime        .clear()                ;
    m_pmtID	     .clear()   	     ;
    m_peTrackID	     .clear()                ;
    m_isCerenkov     .clear()                ;
    m_isReemission   .clear()                ;
    m_isOriginalOP   .clear()                ;
    m_OriginalOPTime .clear()                ;    

      m_localpos_x	.clear();    
      m_localpos_y	.clear();
      m_localpos_z	.clear();

      m_localdir_x	.clear();
      m_localdir_y	.clear();
      m_localdir_z	.clear();

      m_boundarypos_x	.clear();
      m_boundarypos_y	.clear();
      m_boundarypos_z	.clear();    
      
      m_globalpos_x      .clear(); 
      m_globalpos_y      .clear(); 
      m_globalpos_z      .clear(); 
          
      m_PMTID_byPMT     .clear();
      m_nPE_byPMT    .clear();   
      m_cache_bypmt.clear();
     }

}

void
NormalAnaMgr::EndOfEventAction(const G4Event* evt) {


if(m_flag_hitinfo==false)
{
    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    G4int CollID = SDman->GetCollectionID("hitCollection");

    junoHit_PMT_Collection* col = 0;
    G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
    if (!HCE or CollID<0) {
        LogError << "No hits collection found." << std::endl;
    } else {
        col = (junoHit_PMT_Collection*)(HCE->GetHC(CollID));
    }

    int totPE = 0;
    if (col) {
        int n_hit = col->entries();
        m_nPhotons = n_hit;
        if (n_hit > 2000000) { m_nPhotons = 2000000; }
        for (int i = 0; i < n_hit; ++i) {
            totPE += (*col)[i]->GetCount();
        }  
         m_totalPE = totPE;
     }
    

}






if(m_flag_hitinfo==true)
 {
    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    G4int CollID = SDman->GetCollectionID("hitCollection");

    junoHit_PMT_Collection* col = 0; 
    G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
    if (!HCE or CollID<0) {
        LogError << "No hits collection found." << std::endl;
    } else {
        col = (junoHit_PMT_Collection*)(HCE->GetHC(CollID));
    }
    
    // fill evt data
    int totPE = 0;
    if (col) {
        int n_hit = col->entries();
        m_nPhotons = n_hit;
        // FIXME: Make sure not overflow
        if (n_hit > 2000000) { m_nPhotons = 2000000; }

        for (int i = 0; i < n_hit; ++i) {
            totPE += (*col)[i]->GetCount(); 
            // if overflow, don't save anything into the array.
            // but still count the totalPE.
            if (i >= 2000000) { continue; }
            m_energy		.push_back((*col)[i]->GetKineticEnergy());
            m_nPE		.push_back((*col)[i]->GetCount());
            m_hitTime		.push_back((*col)[i]->GetTime());
            m_pmtID		.push_back((*col)[i]->GetPMTID());

            m_cache_bypmt[m_pmtID[i]] += m_nPE[i];

            if ((*col)[i]->IsFromCerenkov()) {
                LogDebug << "+++++ from cerenkov" << std::endl;
                m_isCerenkov.push_back(1);
            }
            else
            {
               m_isCerenkov.push_back(0);
             }
            if ((*col)[i]->IsReemission()) {
                LogDebug << "+++++ reemission" << std::endl;
                m_isReemission.push_back(1);
            }
            else
             {
               m_isReemission.push_back(0);
             }

            m_isOriginalOP	.push_back( (*col)[i]->IsOriginalOP());
            m_OriginalOPTime	.push_back((*col)[i]->GetOriginalOPStartT());
            m_peTrackID		.push_back((*col)[i]->GetProducerID());

            G4ThreeVector local_pos = (*col)[i]->GetPosition();
            m_localpos_x	.push_back(local_pos.x());
            m_localpos_y	.push_back(local_pos.y());
            m_localpos_z	.push_back(local_pos.z());

            G4ThreeVector local_dir = (*col)[i]->GetMomentum();
            m_localdir_x	.push_back(local_dir.x());
            m_localdir_y	.push_back(local_dir.y());
            m_localdir_z	.push_back(local_dir.z());

            G4ThreeVector global_pos = (*col)[i]->GetGlobalPosition();
            m_globalpos_x	.push_back(global_pos.x());
            m_globalpos_y	.push_back(global_pos.y());
            m_globalpos_z	.push_back(global_pos.z());

            G4ThreeVector boundary_pos = (*col)[i]->GetBoundaryPosition();
            m_boundarypos_x	.push_back(boundary_pos.x());
            m_boundarypos_y	.push_back(boundary_pos.y());
            m_boundarypos_z	.push_back(boundary_pos.z());           


        }

    }

    m_npmts_byPMT = 0;
    for (std::map<int,int>::iterator it = m_cache_bypmt.begin();
            it != m_cache_bypmt.end(); ++it) {
        m_PMTID_byPMT	.push_back( it->first);
        m_nPE_byPMT	.push_back(it->second);
        ++m_npmts_byPMT;
    }

    m_totalPE = totPE;
}
  
      if (m_edep>0) {
        m_edep_x /= m_edep;
        m_edep_y /= m_edep;
        m_edep_z /= m_edep;
    }

    if (m_flag_ntuple and m_evt_tree) {
        m_evt_tree -> Fill();
    }
    save_into_data_model();

}


void
NormalAnaMgr::PreUserTrackingAction(const G4Track* aTrack) {

  if(aTrack->GetParentID()==0 && aTrack->GetUserInformation()==0)
    {
        NormalTrackInfo* anInfo = new NormalTrackInfo(aTrack);
        G4Track* theTrack = (G4Track*)aTrack;
        theTrack->SetUserInformation(anInfo);
    }
    NormalTrackInfo* info = (NormalTrackInfo*)(aTrack->GetUserInformation());

    if (!info) {
         return;
    }

    if (aTrack->GetDefinition() == G4OpticalPhoton::Definition()
            and aTrack->GetCreatorProcess()) {
        LogDebug << "###+++ "<< aTrack ->GetCreatorProcess()->GetProcessName() <<std::endl; 
    }

    // original OP
    // set the info 
if(m_flag_hitinfo==true)
  {
    if (aTrack->GetDefinition() == G4OpticalPhoton::Definition()
            and info->isOriginalOP()
            and info->getOriginalOPStartTime() == 0.0) {
        // make sure this track info is not changed before.
        assert(info->getOriginalOPStartTime() == 0.0);
        LogDebug << "------ original OP" << std::endl;
        info->setOriginalOPStartTime(aTrack->GetGlobalTime());
    }

    // An example: Get the parent particle name from track info
    if (info && aTrack->GetDefinition()!=G4OpticalPhoton::Definition()) {
        // const G4String& parent_name = info->getParentName();
        // G4cout << " The parent of " << aTrack->GetDefinition()->GetParticleName()
        //        << " is " << parent_name
        //        << G4endl;
    }

  }

}

void
NormalAnaMgr::PostUserTrackingAction(const G4Track* aTrack) {
 
   if (aTrack->GetParentID() == 0) {
        // this is the primary particle
        const G4ThreeVector& pos = aTrack->GetPosition();
        LogDebug << "!!!Primary Track " << aTrack->GetTrackID() << ": ";
        LogDebug << "+ " << pos.x() << " " << pos.y() << " " << pos.z() << std::endl;
        LogDebug << "+ " << aTrack->GetKineticEnergy() << std::endl;
    }
    G4TrackingManager* tm = G4EventManager::GetEventManager() 
                                            -> GetTrackingManager();
    G4TrackVector* secondaries = tm->GimmeSecondaries();
    if(secondaries)
    {
        NormalTrackInfo* info = (NormalTrackInfo*)(aTrack->GetUserInformation());

        if (!info) {
             return;
        }

        size_t nSeco = secondaries->size();
        if(nSeco>0)
        {
            for(size_t i=0;i<nSeco;i++)
            { 
             
                // make sure the secondaries' track info is empty
                // if already attached track info, skip it.
                if ((*secondaries)[i]->GetUserInformation()) {
                    LogDebug << "The secondary already has user track info. skip creating new one" << std::endl;
                    continue;
                }
                NormalTrackInfo* infoNew = new NormalTrackInfo(info);
               if(m_flag_hitinfo==true)
               {
                // cerekov tag
                if ((*secondaries)[i]->GetCreatorProcess() 
                    and (*secondaries)[i]->GetCreatorProcess()->GetProcessName() == "Cerenkov") {
                    infoNew->setFromCerenkov();
                    LogDebug << "### infoNew->setFromCerenkov()" << std::endl;
                }
                // reemission tag
                // + parent track is an OP
                // + secondary is also an OP
                // + the creator process is Scintillation
                if (aTrack->GetDefinition() == G4OpticalPhoton::Definition()
                    and (*secondaries)[i]->GetDefinition() == G4OpticalPhoton::Definition()
                    and (*secondaries)[i]->GetCreatorProcess()->GetProcessName() == "Scintillation") {
                    infoNew->setReemission();
                }
                // original optical photon tag
                if (aTrack->GetDefinition() != G4OpticalPhoton::Definition() 
                    and (*secondaries)[i]->GetDefinition() == G4OpticalPhoton::Definition()
                    ) {
                    LogDebug << "------ original OP" << std::endl;
                    infoNew->setOriginalOP();
                }
                }
                // save parent track info
                infoNew->setParentName(aTrack->GetDefinition()->GetParticleName());

                (*secondaries)[i]->SetUserInformation(infoNew);
            }
        }
    }
  

}

void
NormalAnaMgr::UserSteppingAction(const G4Step* step) {
    
    G4Track* track = step->GetTrack();
    G4double edep = step->GetTotalEnergyDeposit();

    // cache the material pointer, avoid comparison using string.
    static G4Material* s_mat_LS = nullptr;
    if (!s_mat_LS) {
        s_mat_LS = G4Material::GetMaterial("LS");
    }

    if (edep > 0 and track->GetDefinition()!= G4OpticalPhoton::Definition()
                 and track->GetMaterial() == s_mat_LS) {
        m_edep += edep;
        G4ThreeVector pos = step -> GetPreStepPoint() -> GetPosition();
        m_edep_x += edep * pos.x();
        m_edep_y += edep * pos.y();
        m_edep_z += edep * pos.z();

    }
   
    // if the step number of optical photon bigger than X, mark it as killed
    if (track->GetDefinition() == G4OpticalPhoton::Definition()) {
        G4int stepno = track->GetCurrentStepNumber();
       
      if(m_flag_hitinfo==true)
       {
        if (track->GetTrackStatus() == fStopAndKill) {
            // if the opticalphoton is killed, save the step no
            m_step_no->Fill(stepno);
        }
       }
        if (stepno >= 1000) {
            G4String phyname;
            if (track->GetVolume()) { phyname = track->GetVolume()->GetName(); }
            const G4ThreeVector& tmppos = track->GetPosition();
            LogWarn << "opticalphoton [" << track->GetTrackID() << "]"
                    << "@[" << phyname << "]"
                    << " (" << tmppos.x() << ", "
                    << tmppos.y() << ", "
                    << tmppos.z() << ", "
                    << track->GetGlobalTime() << ") "
                    << " step number >= " << 1000
                    << std::endl;
            track->SetTrackStatus(fStopAndKill);
        }

        // update the last hit acrylic surface
     if(m_flag_hitinfo==true)
      {
        G4StepPoint* prepoint = step->GetPreStepPoint();
        G4StepPoint* postpoint = step->GetPostStepPoint();
        // if(postpoint->GetStepStatus()==fGeomBoundary) {
        //     LogInfo << " * "
        //             << prepoint->GetPhysicalVolume()->GetName() << "/"
        //             << postpoint->GetPhysicalVolume()->GetName() << " "
        //             << std::endl;
        // }

        static G4PhysicalVolumeStore* s_physvol_store = nullptr;
        static G4VPhysicalVolume* s_physvol_pAcrylic = nullptr;
        static G4VPhysicalVolume* s_physvol_pInnerWater = nullptr;
        if (!s_physvol_store) {
            s_physvol_store = G4PhysicalVolumeStore::GetInstance();
            s_physvol_pAcrylic = s_physvol_store->GetVolume("pAcrylic");
            s_physvol_pInnerWater = s_physvol_store->GetVolume("pInnerWater");
        }

        if(postpoint->GetStepStatus()==fGeomBoundary 
        && prepoint->GetPhysicalVolume() == s_physvol_pAcrylic
        && postpoint->GetPhysicalVolume() == s_physvol_pInnerWater
            ) {
            G4ThreeVector prepos = prepoint->GetPosition();
            G4ThreeVector postpos = postpoint->GetPosition();
            // LogInfo << " * "
            //         << prepoint->GetPhysicalVolume()->GetName() << "/"
            //         << postpoint->GetPhysicalVolume()->GetName() << " "
            //         << "(" << prepos.x() << ", " << prepos.y() << ", " << prepos.z() << "), "
            //         << "(" << postpos.x() << ", " << postpos.y() << ", " << postpos.z() << "), "
            //         << std::endl;
            // TODO use postpoint as the boundary between acrylic and water
            NormalTrackInfo* info = (NormalTrackInfo*)(track->GetUserInformation());
            if (info) {
                info->setBoundaryPos(postpos);
            }

        }
    }
  }

}

bool NormalAnaMgr::save_into_data_model() {
    SniperDataPtr<JM::NavBuffer>  navBuf(*getParent(), "/Event");
    if (navBuf.invalid()) {
        return false;
    }
    LogDebug << "navBuf: " << navBuf.data() << std::endl;
    JM::EvtNavigator* evt_nav = navBuf->curEvt();
    LogDebug << "evt_nav: " << evt_nav << std::endl;
    if (not evt_nav) {
        return false;
    }
    JM::SimHeader* m_simheader = dynamic_cast<JM::SimHeader*>(evt_nav->getHeader("/Event/Sim"));
    LogDebug << "simhdr: " << m_simheader << std::endl;
    if (not m_simheader) {
        return false;
    }
    JM::SimEvent* m_simevent = dynamic_cast<JM::SimEvent*>(m_simheader->event());
    LogDebug << "simevt: " << m_simevent << std::endl;
    if (not m_simevent) {
        return false;
    }

    // DO NOTHING
    return true;
}
