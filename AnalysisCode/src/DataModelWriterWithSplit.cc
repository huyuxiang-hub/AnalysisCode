#include "DataModelWriterWithSplit.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperLog.h"
#include "SniperKernel/Incident.h"
#include "SniperKernel/SniperException.h"
#include "SniperKernel/SniperDataPtr.h"


#include "DataRegistritionSvc/DataRegistritionSvc.h"

#include "BufferMemMgr/IDataMemMgr.h"
#include "EvtNavigator/NavBuffer.h"
#include "Event/SimHeader.h"
#include "Event/SimEvent.h"

#include "G4Event.hh"
#include "G4Track.hh"
#include "G4SDManager.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"

#include "junoHit_PMT.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"

#include "DetSimAlg/ISimTrackSvc.h"

DECLARE_TOOL(DataModelWriterWithSplit);

DataModelWriterWithSplit::DataModelWriterWithSplit(const std::string& name)
    : ToolBase(name)
{
    iotaskname = "detsimiotask";
    iotask = 0;
    declProp("HitsMax", m_hitcol_max = 10000); // max hits in one sub event
}

DataModelWriterWithSplit::~DataModelWriterWithSplit()
{

}

void
DataModelWriterWithSplit::BeginOfRunAction(const G4Run* /*aRun*/) {
    // Here, the scope is another I/O Task 
    Task* toptask = getRoot();
    iotask = dynamic_cast<Task*>(toptask->find(iotaskname));
    if (iotask == 0) {
        LogError << "Can't find the task for I/O." << std::endl;
        throw SniperException("Make sure the IO task is created");
    }
    // check the BufferMgr in iotask
    SniperPtr<IDataMemMgr> mMgr(*iotask, "BufferMemMgr");
    if ( mMgr.invalid() ) {
        LogError << "Failed to get BufferMemMgr!" << std::endl;
        throw SniperException("Make sure you have load the BufferMemMgr.");
    }
    
    SniperPtr<ISimTrackSvc> simtracksvc_ptr(getParent(), "SimTrackSvc");
    if (simtracksvc_ptr.invalid()) {
      simtracksvc = dynamic_cast<ISimTrackSvc*>(getParent()->createSvc("SimTrackSvc"));
    } else {
        simtracksvc = simtracksvc_ptr.data();
    }
    // SniperPtr<DataRegistritionSvc> drsSvc(iotask, "DataRegistritionSvc");
    // if ( drsSvc.invalid() ) {
    //     LogError << "Failed to get DataRegistritionSvc!" << std::endl;
    //     throw SniperException("Make sure you have load the DataRegistritionSvc.");
    // }
    // FIXME: Why we need register Data???
    // drsSvc->registerData("JM::SimEvent", "/Event/Sim");
}

void
DataModelWriterWithSplit::EndOfRunAction(const G4Run* /*aRun*/) {

}

void
DataModelWriterWithSplit::BeginOfEventAction(const G4Event* /*evt*/) {

}

void
DataModelWriterWithSplit::EndOfEventAction(const G4Event* evt) {
    // == retrieve the geant4 collection ==
    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    G4int CollID = SDman->GetCollectionID("hitCollection");

    G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
    junoHit_PMT_Collection* col = (junoHit_PMT_Collection*)(HCE->GetHC(CollID));
    if (not col) {
        LogWarn << "Not Hit found" << std::endl;
        return;
    }
    int n_hit = col->entries();
    if(n_hit==0)
       {col=0;}
   
    CollID = SDman->GetCollectionID("hitCollectionMuon");
    junoHit_PMT_muon_Collection* col_muon = (junoHit_PMT_muon_Collection*)(HCE->GetHC(CollID));
    if (col_muon) {
        LogDebug << "size of muon hit collection: " << col_muon->entries()
                 << std::endl;
        n_hit += col_muon->entries();
    }
   
    if(col){col_muon=0;}

    m_start_idx = 0;
    ///=======================//
    int m_nPhotons = 0;
    double m_timewindow = 0.0;
    bool minmax_initialized = false;
    double max_CDLPMT_hittime = 0;
    double min_CDLPMT_hittime = 0;

    if(col_muon){
    
    for (int i = 0; i < n_hit; ++i) {
    int copyno = (*col_muon)[i]->GetPMTID();
     if(copyno < 30000) 
                {
                    if (!minmax_initialized) {
                        max_CDLPMT_hittime = ((*col_muon)[i]->GetTime());
                        min_CDLPMT_hittime = ((*col_muon)[i]->GetTime());
                        minmax_initialized = true;
                    } else { 
                        if ((*col_muon)[i]->GetTime() < min_CDLPMT_hittime) {
                            min_CDLPMT_hittime = (*col_muon)[i]->GetTime();
                        }
                        if ((*col_muon)[i]->GetTime() > max_CDLPMT_hittime) {
                            max_CDLPMT_hittime = (*col_muon)[i]->GetTime();
                        }
                    }
                }
      if( copyno < 30000 || copyno >= 300000)
        { m_nPhotons += (*col_muon)[i]->GetCount();}
         
        }
       m_timewindow = max_CDLPMT_hittime - min_CDLPMT_hittime;
     }
   
    if(col)
    {
       for (int i = 0; i < n_hit; ++i) {
       int copyno = (*col)[i]->GetPMTID();
       if(copyno < 30000)
                {
                    if (!minmax_initialized) {
                        max_CDLPMT_hittime = ((*col)[i]->GetTime());
                        min_CDLPMT_hittime = ((*col)[i]->GetTime());
                        minmax_initialized = true;
                    } else {
                        if ((*col)[i]->GetTime() < min_CDLPMT_hittime) {
                            min_CDLPMT_hittime = (*col)[i]->GetTime();
                        }
                        if ((*col)[i]->GetTime() > max_CDLPMT_hittime) {
                            max_CDLPMT_hittime = (*col)[i]->GetTime();
                        }
                    }
                }

        

        }
       m_nPhotons=col->entries();
       m_timewindow = max_CDLPMT_hittime - min_CDLPMT_hittime;

    
     }
    
    SniperDataPtr<JM::NavBuffer>  navBuf(*getParent(), "/Event");
    if (navBuf.invalid()) {
        LogError << "Can't find the NavBuffer." << std::endl;
        return;
    }
    JM::EvtNavigator* evt_nav = navBuf->curEvt();
    if (not evt_nav) {
        LogError << "Can't find the event navigator." << std::endl;
        return;
    }

    // == get the BufferMemMgr of I/O task ==
    SniperPtr<IDataMemMgr> mMgr(*iotask, "BufferMemMgr");
    // == begin magic ==
    LogInfo << "writing events with split begin. " << TTimeStamp() << std::endl;
    while (m_start_idx<n_hit) {
        LogDebug << "Event idx: " << evt->GetEventID() << std::endl;
        LogDebug << "start idx: " << m_start_idx << std::endl;
        JM::EvtNavigator* nav = new JM::EvtNavigator();
        if(m_start_idx==0){
         TTimeStamp ts = evt_nav->TimeStamp();
         nav->setTimeStamp(ts);

        }else{
        TTimeStamp ts;
        nav->setTimeStamp(ts);
        }
        mMgr->adopt(nav, "/Event");
        
        // == create header ==
        
        JM::SimHeader* sim_header = new JM::SimHeader;
        if (m_start_idx==0){
           sim_header->setCDLPMTtotalHits(m_nPhotons);
           sim_header->setCDLPMTtimeWindow(m_timewindow);
         }
        // == create event ==
        JM::SimEvent* sim_event = new JM::SimEvent(evt->GetEventID());
        // == fill tracks ==
        if (m_start_idx==0 ){     
        collect_primary_track(evt);
        fill_tracks(sim_event);
        }
        // == fill hits ==
        fill_hits(sim_event, evt);
        // == add header ==
        sim_header->setEvent(sim_event);
        nav->addHeader("/Event/Sim", sim_header);
        // == trigger the io event ==
        Incident::fire(*getParent(), iotaskname);
    }
    LogInfo << "writing events with split end. " << TTimeStamp() << std::endl;
    // == end magic ==

}

void 
DataModelWriterWithSplit::fill_hits(JM::SimEvent* dst, const G4Event* evt)
{

    LogDebug << "Begin Fill Hits" << std::endl;
    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    G4int CollID = SDman->GetCollectionID("hitCollection");

    G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
    junoHit_PMT_Collection* col = (junoHit_PMT_Collection*)(HCE->GetHC(CollID));
    // muon hit type
    CollID = SDman->GetCollectionID("hitCollectionMuon");
    junoHit_PMT_muon_Collection* col_muon = (junoHit_PMT_muon_Collection*)(HCE->GetHC(CollID));
    if (col_muon) {
        fill_hits_tmpl(col_muon, dst);
    }
    
    // fill evt data
    // int totPE = 0;
    if (col) {
        fill_hits_tmpl(col, dst);
    }
    LogDebug << "End Fill Hits" << std::endl;

}

void 
DataModelWriterWithSplit::collect_primary_track(const G4Event* evt)
{
    LogDebug << "Begin Fill Tracks" << std::endl;
    std::vector<JM::SimTrack*>& m_tracks=simtracksvc->all();     

    G4int nVertex = evt -> GetNumberOfPrimaryVertex();
    for (G4int index=0; index < nVertex; ++index) {
        G4PrimaryVertex* vtx = evt->GetPrimaryVertex( index );

        // Vertex info
        double x = vtx->GetX0();
        double y = vtx->GetY0();
        double z = vtx->GetZ0();
        double t = vtx->GetT0();

        // Loop Over Particle
        G4PrimaryParticle* pp = vtx -> GetPrimary();

        while (pp) {

            int trkid = pp -> GetTrackID();
            int pdgid = pp -> GetPDGcode();
            double px = pp -> GetPx();
            double py = pp -> GetPy();
            double pz = pp -> GetPz();
            double mass = pp -> GetMass();

            // new track
             for (auto trk: m_tracks) {
             if (trkid == trk->getTrackID()) {
                  trk->setPDGID(pdgid);
                  trk->setInitPx(px);
                  trk->setInitPy(py);
                  trk->setInitPz(pz);
                  trk->setInitMass(mass);
                  trk->setInitX(x);
                  trk->setInitY(y);
                  trk->setInitZ(z);
                  trk->setInitT(t);
                  
                  break;
               }
             }

            pp = pp->GetNext();
        }
    }
    LogDebug << "End Fill Tracks" << std::endl;
}

void
DataModelWriterWithSplit::fill_tracks(JM::SimEvent* dst)
{
    if (!simtracksvc) {
        LogWarn << "SimTrackSvc is not available to save additional tracks" << std::endl;
        return;
    }
    std::vector<JM::SimTrack*>& alltracks = simtracksvc->all();
    for (auto track: alltracks) {
        JM::SimTrack* jm_trk = dst->addTrack();

        jm_trk->setPDGID   (track->getPDGID());
        jm_trk->setTrackID (track->getTrackID());
        jm_trk->setInitMass(track->getInitMass());

        jm_trk->setInitPx  (track->getInitPx());
        jm_trk->setInitPy  (track->getInitPy());
        jm_trk->setInitPz  (track->getInitPz());
        jm_trk->setInitX   (track->getInitX());
        jm_trk->setInitY   (track->getInitY());
        jm_trk->setInitZ   (track->getInitZ());
        jm_trk->setInitT   (track->getInitT());

        jm_trk->setExitPx  (track->getExitPx());
        jm_trk->setExitPy  (track->getExitPy());
        jm_trk->setExitPz  (track->getExitPz());
        jm_trk->setExitX   (track->getExitX());
        jm_trk->setExitY   (track->getExitY());
        jm_trk->setExitZ   (track->getExitZ());
        jm_trk->setExitT   (track->getExitT());

        jm_trk->setTrackLength(track->getTrackLength());

        jm_trk->setEdep    (track->getEdep());
        jm_trk->setEdepX   (track->getEdepX());
        jm_trk->setEdepY   (track->getEdepY());
        jm_trk->setEdepZ   (track->getEdepZ());
        
        jm_trk->setQEdep   (track->getQEdep());
        jm_trk->setQEdepX  (track->getQEdepX());
        jm_trk->setQEdepY  (track->getQEdepY());
        jm_trk->setQEdepZ  (track->getQEdepZ());

        jm_trk->setEdepNotInLS(track->getEdepNotInLS());
    }

    simtracksvc->reset();
}

