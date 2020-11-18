#include "RadioAnaMgr.hh"

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4EventManager.hh"
#include "G4TrackingManager.hh"
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/SniperDataPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperLog.h"
#include "RootWriter/RootWriter.h"

#include "NormalTrackInfo.hh"

#include "TROOT.h"


DECLARE_TOOL(RadioAnaMgr);

RadioAnaMgr::RadioAnaMgr(const std::string& name)
    : ToolBase(name)
{
    declProp("StopAtPa234m", m_StopAtPa234m=true);
}

RadioAnaMgr::~RadioAnaMgr() {

}

void RadioAnaMgr::BeginOfRunAction(const G4Run*) {
    SniperPtr<RootWriter> svc(*getParent(), "RootWriter");
    if (svc.invalid()) {
        LogError << "Can't Locate RootWriter. If you want to use it, please "
                 << "enalbe it in your job option file."
                 << std::endl;
        return;
    }
    gROOT->ProcessLine("#include <vector>");

    m_radio_tree = svc->bookTree("SIMEVT/radio", "radioactivity decay");
    m_radio_tree->Branch("evtID", &m_evtid);
    m_radio_tree->Branch("radioidx", &m_radioidx);
    m_radio_tree->Branch("trkid", &m_trkid);
    m_radio_tree->Branch("parentid", &m_parentid);
    m_radio_tree->Branch("pdgid", &m_pdgid);
    m_radio_tree->Branch("name", &m_name);
    m_radio_tree->Branch("time", &m_time);
    m_radio_tree->Branch("Ek", &m_Ek);
}

void RadioAnaMgr::EndOfRunAction(const G4Run*) {

}

void RadioAnaMgr::BeginOfEventAction(const G4Event* anEvent) {
    // reset all variables
    m_evtid = anEvent->GetEventID();
    m_radioidx.clear();
    m_trkid   .clear();
    m_parentid.clear();
    m_pdgid   .clear();
    m_name    .clear();
    m_time    .clear();
    m_Ek      .clear();

    current_radioidx = -1;

}

void RadioAnaMgr::EndOfEventAction(const G4Event*) {
    m_radio_tree->Fill();

}

void RadioAnaMgr::PreUserTrackingAction(const G4Track* aTrack) {
    if (aTrack->GetDefinition() == G4OpticalPhoton::Definition()) {
        return;
    }

    if (aTrack->GetDefinition()->IsGeneralIon()) {
        LogDebug << "- Track: particle name " << (aTrack->GetDefinition()->GetParticleName())
                 << " pdgcode: " << (aTrack->GetDefinition()->GetPDGEncoding())
                 << std::endl;
    }

    const G4VProcess* proc = aTrack->GetCreatorProcess();
    if (!proc) {
        return;
    }

    // TODO: save Radio related
}

void RadioAnaMgr::PostUserTrackingAction(const G4Track* aTrack) {

}

void RadioAnaMgr::UserSteppingAction(const G4Step* aStep) {
    const G4Track* aTrack = aStep->GetTrack();
    if (aTrack->GetDefinition() == G4OpticalPhoton::Definition()) {
        return;
    }

    G4TrackingManager* tm = G4EventManager::GetEventManager() 
                                            -> GetTrackingManager();
    G4TrackVector* secondaries = tm->GimmeSecondaries();
    if(not secondaries) {
        return;
    }

    size_t nSeco = secondaries->size();

    bool is_radio_decay = false;
    // we don't want Geant4 simulate decay of Pa234[73.920X]
    bool kill_rad_secondaries = false;

    // The first case is when Pa234m is not a primary particle,
    // it will stop decay. 
    // This assumes we only want to simulation Th234 -> Pa234m.

    if (aTrack->GetParentID() != 0 
        && aTrack->GetDefinition()->GetParticleName() == "Pa234[73.920X]"
        && m_StopAtPa234m) {
        kill_rad_secondaries = true;
    }


    for (size_t i = 0; i < nSeco; ++i) {
        G4Track* sectrk = (*secondaries)[i];
        const G4VProcess* creatorProcess = sectrk->GetCreatorProcess();

        if (!creatorProcess) {
            continue;
        }

        if (creatorProcess->GetProcessName() == "RadioactiveDecay") {
            is_radio_decay = true;

            // kill 
            if (kill_rad_secondaries) {
                sectrk->SetTrackStatus(fStopAndKill);
            }
        }
        // LogInfo << "creatorProcess: " << creatorProcess->GetProcessName() << std::endl;
    }

    // if there is no radioactivedecay, just return.
    if (!is_radio_decay) {
        return;
    }
    // if the secondaries are killed, just return.
    if (kill_rad_secondaries) {
        return;
    }
        
    // Save all the radioactivity decay related. 
    // Note, even though sometimes we don't want Pa234m,
    // but it will be still simulated. So it will be also recorded.
    // Maybe we could check the status of track is fStopAndKill.

    ++current_radioidx;

    // parent
    m_radioidx.push_back(current_radioidx);
    m_trkid.push_back(aTrack->GetTrackID());
    m_parentid.push_back(aTrack->GetParentID());
    m_pdgid.push_back(aTrack->GetDefinition()->GetPDGEncoding());
    m_name.push_back(aTrack->GetDefinition()->GetParticleName());
    m_time.push_back(aTrack->GetGlobalTime());
    m_Ek.push_back(aTrack->GetKineticEnergy());
    LogDebug << "--> Parent pdgcode " << m_pdgid.back()
             << " time: " << m_time.back()
             << std::endl;

    // daughters
    for (size_t i = 0; i < nSeco; ++i) {
        G4Track* sectrk = (*secondaries)[i];
        const G4VProcess* creatorProcess = sectrk->GetCreatorProcess();

        if (!creatorProcess) {
            continue;
        }

        if (creatorProcess->GetProcessName() == "RadioactiveDecay") {

            m_radioidx.push_back(current_radioidx);
            m_trkid.push_back(sectrk->GetTrackID());
            m_parentid.push_back(sectrk->GetParentID());
            m_pdgid.push_back(sectrk->GetDefinition()->GetPDGEncoding());
            m_name.push_back(sectrk->GetDefinition()->GetParticleName());
            m_time.push_back(sectrk->GetGlobalTime());
            m_Ek.push_back(sectrk->GetKineticEnergy());
            LogDebug << "--> Daughter pdgcode " << m_pdgid.back()
                     << " time: " << m_time.back()
                     << " Ek: " << m_Ek.back()
                     << std::endl;

        }
        // LogInfo << "creatorProcess: " << creatorProcess->GetProcessName() << std::endl;
    }

}
