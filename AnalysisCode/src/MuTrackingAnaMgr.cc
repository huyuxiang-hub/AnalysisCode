#include "MuTrackingAnaMgr.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperLog.h"
#include "RootWriter/RootWriter.h"

#include <G4Event.hh>
#include <G4SteppingManager.hh>
#include <G4SystemOfUnits.hh>

#include <TROOT.h>
#include <TString.h>

DECLARE_TOOL(MuTrackingAnaMgr);

MuTrackingAnaMgr::MuTrackingAnaMgr(const std::string& name)
    : ToolBase(name)
{
    declProp("DeviationThreshold", deviation_threshold = 1e-3); // rad
    // Note: 1e-3 rad \sim 0.06deg
    declProp("MomentumThreshold",  momentum_threshold  = 0*MeV );

    m_MuPosx = new std::vector<float>;
    m_MuPosy = new std::vector<float>;
    m_MuPosz = new std::vector<float>;
    m_MuPx   = new std::vector<float>;
    m_MuPy   = new std::vector<float>;
    m_MuPz   = new std::vector<float>;
}

MuTrackingAnaMgr::~MuTrackingAnaMgr(){
    delete m_MuPosx;
    delete m_MuPosy;
    delete m_MuPosz;
    delete m_MuPx  ;
    delete m_MuPy  ;
    delete m_MuPz  ;
}

void MuTrackingAnaMgr::BeginOfRunAction(const G4Run* /*aRun*/) {
    SniperPtr<RootWriter> svc(*getParent(), "RootWriter");
    gROOT->ProcessLine("#include <vector>");
    m_muon_track_tree = svc->bookTree("SIMEVT/mu_tracking", "muon tracking");
    m_muon_track_tree->Branch("evtID"      , &m_eventID    , "evtID/I"      );
    m_muon_track_tree->Branch("MuTrackID"  , &m_MuTrackID  , "MuTrackID/I"  );
    m_muon_track_tree->Branch("MuParentID" , &m_MuParentID , "MuParentID/I" );
    m_muon_track_tree->Branch("Mu_Posx"    , &m_MuPosx);
    m_muon_track_tree->Branch("Mu_Posy"    , &m_MuPosy);
    m_muon_track_tree->Branch("Mu_Posz"    , &m_MuPosz);
    m_muon_track_tree->Branch("Mu_Px"      , &m_MuPx  );
    m_muon_track_tree->Branch("Mu_Py"      , &m_MuPy  );
    m_muon_track_tree->Branch("Mu_Pz"      , &m_MuPz  );
}

void MuTrackingAnaMgr::EndOfRunAction(const G4Run* /*aRun*/) {
}

void MuTrackingAnaMgr::BeginOfEventAction(const G4Event* evt) {
    m_eventID = evt->GetEventID();

    for(unsigned int iMu = 0; iMu < m_v_MuTrackID.size(); ++iMu){
        m_MuTrackID   = m_v_MuTrackID  [iMu];
        delete m_v_MuPosx[m_MuTrackID];
        delete m_v_MuPosy[m_MuTrackID];
        delete m_v_MuPosz[m_MuTrackID];
        delete m_v_MuPx  [m_MuTrackID];
        delete m_v_MuPy  [m_MuTrackID];
        delete m_v_MuPz  [m_MuTrackID];
    }
    m_v_MuTrackID  .clear();
    m_v_MuParentID .clear();
    m_v_MuPosx.clear();
    m_v_MuPosy.clear();
    m_v_MuPosz.clear();
    m_v_MuPx  .clear();
    m_v_MuPy  .clear();
    m_v_MuPz  .clear();
    m_v_MuExitPosx.clear();
    m_v_MuExitPosy.clear();
    m_v_MuExitPosz.clear();
    m_v_MuExitPx  .clear();
    m_v_MuExitPy  .clear();
    m_v_MuExitPz  .clear();

    m_trackID_lowP.clear();
}

void MuTrackingAnaMgr::EndOfEventAction(const G4Event* /*evt*/) {
    for(unsigned int iMu = 0; iMu < m_v_MuTrackID.size(); ++iMu){
        m_MuTrackID   = m_v_MuTrackID  [iMu];
        m_MuParentID  = m_v_MuParentID [iMu];
        *m_MuPosx     = *(m_v_MuPosx[m_MuTrackID]);
        *m_MuPosy     = *(m_v_MuPosy[m_MuTrackID]);
        *m_MuPosz     = *(m_v_MuPosz[m_MuTrackID]);
        *m_MuPx       = *(m_v_MuPx  [m_MuTrackID]);
        *m_MuPy       = *(m_v_MuPy  [m_MuTrackID]);
        *m_MuPz       = *(m_v_MuPz  [m_MuTrackID]);
        // Add exit point to segment
        m_MuPosx->push_back(m_v_MuExitPosx[m_MuTrackID]);
        m_MuPosy->push_back(m_v_MuExitPosy[m_MuTrackID]);
        m_MuPosz->push_back(m_v_MuExitPosz[m_MuTrackID]);
        m_MuPx  ->push_back(m_v_MuExitPx  [m_MuTrackID]);
        m_MuPy  ->push_back(m_v_MuExitPy  [m_MuTrackID]);
        m_MuPz  ->push_back(m_v_MuExitPz  [m_MuTrackID]);

        m_muon_track_tree->Fill();
    }

    if(m_trackID_lowP.size() > 0){
        LogInfo << "Number of muons below threshold: " << m_trackID_lowP.size() << std::endl;
    }
}

void MuTrackingAnaMgr::UserSteppingAction(const G4Step* fStep) {
    G4Track* fTrack = fStep->GetTrack();
    G4int trackId   = fTrack->GetTrackID();
    G4int StepNo    = fTrack->GetCurrentStepNumber();

    G4int trackPDG  = fTrack->GetDefinition()->GetPDGEncoding();

    // only deals with muons
    if(fabs(trackPDG)!=13){
        return;
    }

    G4StepPoint* thePrePoint  = fStep->GetPreStepPoint() ;
    G4StepPoint* thePostPoint = fStep->GetPostStepPoint();

    G4int parId = fTrack->GetParentID();

    if(StepNo==1)   { // first step of this muon
        if(parId != 0 && thePrePoint->GetMomentum().mag() < momentum_threshold){
            // Do not save secondary muons with low momentum
            m_trackID_lowP.insert(trackId);
            return;
        }

        m_v_MuTrackID .push_back(trackId);
        m_v_MuParentID.push_back(parId);

        m_v_MuPosx[trackId] = new std::vector<float>;
        m_v_MuPosy[trackId] = new std::vector<float>;
        m_v_MuPosz[trackId] = new std::vector<float>;
        m_v_MuPx  [trackId] = new std::vector<float>;
        m_v_MuPy  [trackId] = new std::vector<float>;
        m_v_MuPz  [trackId] = new std::vector<float>;

        m_v_MuPosx[trackId]->push_back(thePrePoint->GetPosition().x()/mm);
        m_v_MuPosy[trackId]->push_back(thePrePoint->GetPosition().y()/mm);
        m_v_MuPosz[trackId]->push_back(thePrePoint->GetPosition().z()/mm);
        m_v_MuPx  [trackId]->push_back(thePrePoint->GetMomentum().x()   );
        m_v_MuPy  [trackId]->push_back(thePrePoint->GetMomentum().y()   );
        m_v_MuPz  [trackId]->push_back(thePrePoint->GetMomentum().z()   );
    }
    else {
        if(m_trackID_lowP.find(trackId) != m_trackID_lowP.end()){
            // Track with low momentum => skip
            return;
        }
        if(m_v_MuPosx.find(trackId) == m_v_MuPosx.end()){
            LogError << "Could not find muon with track ID " << trackId << "." <<
                "This is not first step (step " << StepNo << ")." << std::endl;
            throw;
        }
    }

    m_v_MuExitPosx[trackId] = thePostPoint->GetPosition().x()/mm;
    m_v_MuExitPosy[trackId] = thePostPoint->GetPosition().y()/mm;
    m_v_MuExitPosz[trackId] = thePostPoint->GetPosition().z()/mm;
    m_v_MuExitPx  [trackId] = thePostPoint->GetMomentum().x()   ;
    m_v_MuExitPy  [trackId] = thePostPoint->GetMomentum().y()   ;
    m_v_MuExitPz  [trackId] = thePostPoint->GetMomentum().z()   ;

    CLHEP::Hep3Vector last_dir(m_v_MuPx[trackId]->back(),
                               m_v_MuPy[trackId]->back(),
                               m_v_MuPz[trackId]->back() );

    last_dir = last_dir.unit();
    if(last_dir.mag() > 0){
        float deviation = acos(last_dir.dot(thePostPoint->GetMomentum().unit()));
        if(deviation > deviation_threshold){
            m_v_MuPosx[trackId]->push_back(m_v_MuExitPosx[trackId]);
            m_v_MuPosy[trackId]->push_back(m_v_MuExitPosy[trackId]);
            m_v_MuPosz[trackId]->push_back(m_v_MuExitPosz[trackId]);
            m_v_MuPx  [trackId]->push_back(m_v_MuExitPx  [trackId]);
            m_v_MuPy  [trackId]->push_back(m_v_MuExitPy  [trackId]);
            m_v_MuPz  [trackId]->push_back(m_v_MuExitPz  [trackId]);
        }
    }
}

// vim: et:sts=4:ts=4:sw=4
