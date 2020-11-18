#include "DataModelWriterWithSplit.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperLog.h"
#include "SniperKernel/Incident.h"
#include "SniperKernel/SniperException.h"

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


DECLARE_TOOL(DataModelWriterWithSplit);

DataModelWriterWithSplit::DataModelWriterWithSplit(const std::string& name)
    : ToolBase(name)
{
    iotaskname = "detsimiotask";
    iotask = 0;
    declProp("HitsMax", m_hitcol_max = 100); // max hits in one sub event
}

DataModelWriterWithSplit::~DataModelWriterWithSplit()
{

}

void
DataModelWriterWithSplit::BeginOfRunAction(const G4Run* /*aRun*/) {
    // Here, the scope is another I/O Task 
    Task* toptask = getParent()->getParent();
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
    CollID = SDman->GetCollectionID("hitCollectionMuon");
    junoHit_PMT_muon_Collection* col_muon = (junoHit_PMT_muon_Collection*)(HCE->GetHC(CollID));
    if (col_muon) {
        LogDebug << "size of muon hit collection: " << col_muon->entries()
                 << std::endl;
        n_hit += col_muon->entries();
    }
    m_start_idx = 0;
    // == get the BufferMemMgr of I/O task ==
    SniperPtr<IDataMemMgr> mMgr(*iotask, "BufferMemMgr");
    // == begin magic ==
    LogInfo << "writing events with split begin. " << TTimeStamp() << std::endl;
    while (m_start_idx<n_hit) {
        LogDebug << "Event idx: " << evt->GetEventID() << std::endl;
        LogDebug << "start idx: " << m_start_idx << std::endl;
        JM::EvtNavigator* nav = new JM::EvtNavigator();
        TTimeStamp ts;
        nav->setTimeStamp(ts);
        mMgr->adopt(nav, "/Event");
        // == create header ==
        JM::SimHeader* sim_header = new JM::SimHeader;
        // == create event ==
        JM::SimEvent* sim_event = new JM::SimEvent(evt->GetEventID());
        // == fill hits ==
        fill_hits(sim_event, evt);
        // == fill tracks ==
        fill_tracks(sim_event, evt);
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
DataModelWriterWithSplit::fill_tracks(JM::SimEvent* dst, const G4Event* evt)
{
    LogDebug << "Begin Fill Tracks" << std::endl;

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
            JM::SimTrack* jm_trk = dst->addTrack();
            jm_trk->setPDGID(pdgid);
            jm_trk->setTrackID(trkid);
            jm_trk->setInitPx(px);
            jm_trk->setInitPy(py);
            jm_trk->setInitPz(pz);
            jm_trk->setInitMass(mass);
            jm_trk->setInitX(x);
            jm_trk->setInitY(y);
            jm_trk->setInitZ(z);
            jm_trk->setInitT(t);

            pp = pp->GetNext();
        }
    }
    LogDebug << "End Fill Tracks" << std::endl;
}
