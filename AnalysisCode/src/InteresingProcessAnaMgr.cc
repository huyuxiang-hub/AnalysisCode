
#include "InteresingProcessAnaMgr.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperLog.h"
#include "RootWriter/RootWriter.h"

#include "G4SDManager.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4EventManager.hh"
#include "G4TrackingManager.hh"
#include "G4SteppingManager.hh"
#include "G4VProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4Positronium.hh"
#include "G4Positron.hh"
#include "G4Electron.hh"
#include "G4AntiNeutrinoE.hh"
#include "TMath.h"

#include "NormalTrackInfo.hh"

#include "TROOT.h"

#include "DetSimAlg/ISimTrackSvc.h"
#include "Event/SimTrack.h"

DECLARE_TOOL(InteresingProcessAnaMgr);

InteresingProcessAnaMgr::InteresingProcessAnaMgr(const std::string& name)
    : ToolBase(name), simtracksvc(NULL)
{

}

InteresingProcessAnaMgr::~InteresingProcessAnaMgr()
{

}

// == Run Action
void
InteresingProcessAnaMgr::BeginOfRunAction(const G4Run* /*aRun*/) {
    SniperPtr<RootWriter> svc(*getParent(), "RootWriter");

    gROOT->ProcessLine("#include <vector>");

    m_evt_tree = svc->bookTree("SIMEVT/michael", "evt");
    m_evt_tree -> Branch("evtID", &m_eventID, "evtID/I");
    m_evt_tree -> Branch("pdgid", &michael_pdgid);
    m_evt_tree -> Branch("kine", &michael_kine);
    m_evt_tree -> Branch("px", &michael_px);
    m_evt_tree -> Branch("py", &michael_py);
    m_evt_tree -> Branch("pz", &michael_pz);
    m_evt_tree -> Branch("x", &michael_pos_x);
    m_evt_tree -> Branch("y", &michael_pos_y);
    m_evt_tree -> Branch("z", &michael_pos_z);
    m_evt_tree -> Branch("t", &michael_t);

    m_neutron_tree = svc->bookTree("SIMEVT/nCapture", "Neutron Capture");
    m_neutron_tree -> Branch("evtID", &m_eventID, "evtID/I");
    m_neutron_tree -> Branch("NeutronTrkid", &neutron_only_trkid);
    m_neutron_tree -> Branch("NeutronKine",  &neutron_only_kine );
    m_neutron_tree -> Branch("NeutronCaptureT", &neutron_only_t );
    m_neutron_tree -> Branch("NCStartX", &neutron_only_start_x );
    m_neutron_tree -> Branch("NCStartY", &neutron_only_start_y );
    m_neutron_tree -> Branch("NCStartZ", &neutron_only_start_z );
    m_neutron_tree -> Branch("NCStopX", &neutron_only_stop_x );
    m_neutron_tree -> Branch("NCStopY", &neutron_only_stop_y );
    m_neutron_tree -> Branch("NCStopZ", &neutron_only_stop_z );
    m_neutron_tree -> Branch("NCTrackLength", &neutron_only_track_length );

    m_neutron_tree -> Branch("trkid", &neutron_trkid);
    m_neutron_tree -> Branch("pdgid", &neutron_pdgid);
    m_neutron_tree -> Branch("kine", &neutron_kine);
    m_neutron_tree -> Branch("px", &neutron_px);
    m_neutron_tree -> Branch("py", &neutron_py);
    m_neutron_tree -> Branch("pz", &neutron_pz);
    m_neutron_tree -> Branch("x", &neutron_pos_x);
    m_neutron_tree -> Branch("y", &neutron_pos_y);
    m_neutron_tree -> Branch("z", &neutron_pos_z);
    m_neutron_tree -> Branch("t", &neutron_t);


    // SimTrackSvc
    SniperPtr<ISimTrackSvc> simtracksvc_ptr(getParent(), "SimTrackSvc");
    if (simtracksvc_ptr.invalid()) {
        // create the service
        simtracksvc = dynamic_cast<ISimTrackSvc*>(getParent()->createSvc("SimTrackSvc"));
    } else {
        simtracksvc = simtracksvc_ptr.data();
    }

}

void
InteresingProcessAnaMgr::EndOfRunAction(const G4Run* /*aRun*/) {

}

// == Event Action
void
InteresingProcessAnaMgr::BeginOfEventAction(const G4Event* evt) {
    m_eventID = evt->GetEventID();
    michael_pdgid.clear();
    michael_kine .clear();
    michael_px   .clear();
    michael_py   .clear();
    michael_pz   .clear();
    michael_pos_x.clear();
    michael_pos_y.clear();
    michael_pos_z.clear();
    michael_t    .clear();

    neutron_only_trkid.clear();
    neutron_only_kine .clear();
    neutron_only_t    .clear();
    neutron_only_start_x.clear();
    neutron_only_start_y.clear();
    neutron_only_start_z.clear();
    neutron_only_stop_x.clear();
    neutron_only_stop_y.clear();
    neutron_only_stop_z.clear();
    neutron_only_track_length.clear();

    neutron_trkid.clear();
    neutron_pdgid.clear();
    neutron_kine .clear();
    neutron_px   .clear();
    neutron_py   .clear();
    neutron_pz   .clear();
    neutron_pos_x.clear();
    neutron_pos_y.clear();
    neutron_pos_z.clear();
    neutron_t    .clear();

}

void
InteresingProcessAnaMgr::EndOfEventAction(const G4Event* /*evt*/) {

    m_evt_tree -> Fill();
    m_neutron_tree -> Fill();
}

// == Tracking Action
void
InteresingProcessAnaMgr::PreUserTrackingAction(const G4Track* aTrack) {
    NormalTrackInfo* info = (NormalTrackInfo*)(aTrack->GetUserInformation());
    // LogDebug << "NormalTrackInfo: " << info << std::endl;
    if (info) {
        // select michael electron
        if (info->isMichaelElectron()) {
            saveMichael(aTrack);
        }
        if (info->isNeutronCapture()) {
            saveNeutronCapture(aTrack);
        }

        if (info->getDecayMode()) {
            // if it is the secondaries of g4decay
            saveDecayInit(aTrack);
        }

    }
}

void
InteresingProcessAnaMgr::PostUserTrackingAction(const G4Track* aTrack) {
    NormalTrackInfo* info = (NormalTrackInfo*)(aTrack->GetUserInformation());

    selectMichael(aTrack);
    selectNeutronCapture(aTrack);
    selectDecay(aTrack);

    if (info) {
        // as secondaries, store the exit info
        saveNeutronCaptureSecondary(aTrack);

        if (info->getDecayMode()) {
            // if it is the secondaries of g4decay
            saveDecayExit(aTrack);
        }
    }
}

// == Stepping Action
void
InteresingProcessAnaMgr::UserSteppingAction(const G4Step* /*step*/) {

}

// Helper
// == Select Michael
void 
InteresingProcessAnaMgr::selectMichael(const G4Track* aTrack) {
    G4ParticleDefinition* particle = aTrack->GetDefinition();
    G4int pdgid = particle->GetPDGEncoding();
    // only select the mu- and mu+
    if (pdgid != 13 and pdgid != -13) {
        return;
    }
    // get the secondaries
    G4TrackingManager* tm = G4EventManager::GetEventManager() 
                                            -> GetTrackingManager();
    G4TrackVector* secondaries = tm->GimmeSecondaries();
    if(not secondaries) {
        return;
    }
    size_t nSeco = secondaries->size();

    bool find_michael = false;
    int idx_michael = -1;
    G4Track* trk_michael = 0;
    // what particle is to select
    // -- mu-
    // --- Capture at Rest
    if (pdgid == 13) {
        int idx_anti_nu_e = -1;
        double kine_michael_electron = -1;
        for (size_t i = 0; i < nSeco; ++i) {
            G4Track* sectrk = (*secondaries)[i];
            G4ParticleDefinition* secparticle = sectrk->GetDefinition();
            // const G4String& sec_name = secparticle->GetParticleName();
            const G4VProcess* creatorProcess = sectrk->GetCreatorProcess();
            G4double sec_kine = sectrk->GetKineticEnergy();

            if (creatorProcess->GetProcessName()!="muMinusCaptureAtRest") {
                continue;
            }
            if (secparticle == G4AntiNeutrinoE::Definition() ) {
                idx_anti_nu_e = i;
            } else if (secparticle == G4Electron::Definition() and sec_kine > kine_michael_electron) {
                idx_michael = i;
                trk_michael = sectrk;
            }
        }
        // find the michael electron
        if (idx_michael != -1 and idx_anti_nu_e != -1 and trk_michael) {
            find_michael = true;
        }

    } else if (pdgid == -13) {
        for (size_t i = 0; i < nSeco; ++i) {
            G4Track* sectrk = (*secondaries)[i];
            G4ParticleDefinition* secparticle = sectrk->GetDefinition();
            // const G4String& sec_name = secparticle->GetParticleName();
            const G4VProcess* creatorProcess = sectrk->GetCreatorProcess();

            if (creatorProcess->GetProcessName()!="Decay") {
                continue;
            }
            if (secparticle == G4Positron::Definition()) {
                idx_michael = i;
                trk_michael = sectrk;
            }

        }
        // find the michael positron
        if (idx_michael != -1 and trk_michael) {
            find_michael = true;
        }

    }

    // Fill the info
    if (not find_michael) {
        return;
    }
    NormalTrackInfo* info = (NormalTrackInfo*)(trk_michael->GetUserInformation());
    if (not info) {
        LogWarn << "--- update the track. create new Track Info" << std::endl;
        info = new NormalTrackInfo(aTrack);
        trk_michael->SetUserInformation(info);
    }
    info -> setMichaelElectron();

}

// == save Michael Electron
void
InteresingProcessAnaMgr::saveMichael(const G4Track* aTrack) {
    NormalTrackInfo* info = (NormalTrackInfo*)(aTrack->GetUserInformation());
    if ((not info) or (not info->isMichaelElectron())) {
        return;
    }
    LogDebug << "!-- Michael Electron: " << std::endl;
    // Please note, in Post Tracking Action, we can't get the right
    // track ID of the secondaries.
    LogDebug << "!--- TrkID: " << aTrack->GetTrackID() << std::endl;
    LogDebug << "!--- Kine : " << aTrack->GetKineticEnergy() << std::endl;

    michael_pdgid.push_back(aTrack->GetDefinition()->GetPDGEncoding());
    michael_kine .push_back(aTrack->GetKineticEnergy());
    const G4ThreeVector& mom = aTrack->GetMomentum();
    michael_px   .push_back(mom.x());
    michael_py   .push_back(mom.y());
    michael_pz   .push_back(mom.z());
    const G4ThreeVector& pos = aTrack->GetPosition();
    michael_pos_x.push_back(pos.x());
    michael_pos_y.push_back(pos.y());
    michael_pos_z.push_back(pos.z());
    michael_t    .push_back(aTrack->GetGlobalTime());

    LogDebug << "!---- Kine : " << michael_kine .back() << std::endl;
    LogDebug << "!---- Mome : " << michael_px.back()
                        << ", " << michael_py.back()
                        << ", " << michael_pz.back() << std::endl;
    LogDebug << "!---- Posi : " << michael_pos_x.back()
                        << ", " << michael_pos_y.back()
                        << ", " << michael_pos_z.back() << std::endl;
    LogDebug << "!---- Time : " << michael_t    .back() << std::endl;

}
// == Select Neutron Capture
void
InteresingProcessAnaMgr::selectNeutronCapture(const G4Track* aTrack) {
    G4ParticleDefinition* particle = aTrack->GetDefinition();
    G4int pdgid = particle->GetPDGEncoding();
    // only select neutron
    if (pdgid != 2112) {
        return;
    }
    // get the secondaries
    G4TrackingManager* tm = G4EventManager::GetEventManager() 
                                            -> GetTrackingManager();
    G4TrackVector* secondaries = tm->GimmeSecondaries();
    if(not secondaries) {
        return;
    }
    size_t nSeco = secondaries->size();

    bool find_neutron_capture = false;

    for (size_t i = 0; i < nSeco; ++i) {
        G4Track* sectrk = (*secondaries)[i];
        G4ParticleDefinition* secparticle = sectrk->GetDefinition();
        const G4String& sec_name = secparticle->GetParticleName();
        const G4VProcess* creatorProcess = sectrk->GetCreatorProcess();

        if (creatorProcess->GetProcessName()!="nCapture") {
            continue;
        }
        if (not find_neutron_capture) {
            find_neutron_capture = true;
        }
        LogDebug << "*--- Neutron Capture: " << std::endl;
        LogDebug << "*---- particle: " << sec_name << std::endl;
        LogDebug << "*---- kinetic : " << sectrk->GetKineticEnergy() << std::endl;
        // set flag
        NormalTrackInfo* info = (NormalTrackInfo*)(sectrk->GetUserInformation());

        if (not info) {
            LogWarn << "--- get Track Info failed, create new Track Info" << std::endl;
            info = new NormalTrackInfo(aTrack);
            sectrk->SetUserInformation(info);
        }

        info -> setNeutronCapture();
    }

    // save the current neutron
    if (find_neutron_capture) {
        // save the current track
        neutron_only_trkid.push_back(aTrack->GetTrackID());
        neutron_only_kine .push_back(aTrack->GetVertexKineticEnergy());
        neutron_only_t    .push_back(aTrack->GetGlobalTime());
        const G4ThreeVector& start_pos = aTrack->GetVertexPosition();
        const G4ThreeVector& stop_pos  = aTrack->GetPosition();
        neutron_only_start_x.push_back(start_pos.x());
        neutron_only_start_y.push_back(start_pos.y());
        neutron_only_start_z.push_back(start_pos.z());
        neutron_only_stop_x.push_back(stop_pos.x());
        neutron_only_stop_y.push_back(stop_pos.y());
        neutron_only_stop_z.push_back(stop_pos.z());
        neutron_only_track_length.push_back(aTrack->GetTrackLength());

        LogDebug << "***--- Neutron Capture: " << std::endl;
        LogDebug << "***---- Vkinetic: " << aTrack->GetVertexKineticEnergy() << std::endl;
        LogDebug << "***----  kinetic: " << aTrack->GetKineticEnergy() << std::endl;
        LogDebug << "***---- StartPos: " << start_pos.x() << " " 
                                         << start_pos.y() << " "
                                         << start_pos.z() << std::endl;
        LogDebug << "***---- Stop Pos: " << stop_pos.x() << " " 
                                         << stop_pos.y() << " "
                                         << stop_pos.z() << std::endl;
        LogDebug << "***---- Length  : " << aTrack->GetTrackLength() << std::endl;

    }
}
// == save Neutron Capture
void
InteresingProcessAnaMgr::saveNeutronCapture(const G4Track* aTrack) {
    NormalTrackInfo* info = (NormalTrackInfo*)(aTrack->GetUserInformation());
    if ((not info) or (not info->isNeutronCapture())) {
        return;
    }

    // Tag the current track ID will be traced,
    // so that the edep will be accumulated in DepositEnergyAnaMgr

    LogDebug << "!-- Neutron Capture: " << std::endl;
    LogDebug << "!--- TrkID: " << aTrack->GetTrackID() << std::endl;
    LogDebug << "!--- Kine : " << aTrack->GetKineticEnergy() << std::endl;

    neutron_trkid.push_back(aTrack->GetParentID());
    neutron_pdgid.push_back(aTrack->GetDefinition()->GetPDGEncoding());
    neutron_kine .push_back(aTrack->GetKineticEnergy());
    const G4ThreeVector& mom = aTrack->GetMomentum();
    neutron_px   .push_back(mom.x());
    neutron_py   .push_back(mom.y());
    neutron_pz   .push_back(mom.z());
    const G4ThreeVector& pos = aTrack->GetPosition();
    neutron_pos_x.push_back(pos.x());
    neutron_pos_y.push_back(pos.y());
    neutron_pos_z.push_back(pos.z());
    neutron_t    .push_back(aTrack->GetGlobalTime());

    saveSecondaryInit(aTrack);
}

void
InteresingProcessAnaMgr::saveNeutronCaptureSecondary(const G4Track* aTrack) {
    NormalTrackInfo* info = (NormalTrackInfo*)(aTrack->GetUserInformation());
    if ((not info) or (not info->isNeutronCapture())) {
        return;
    }

    saveSecondaryExit(aTrack);
}

// ===========================================================================
// Decay related
// ===========================================================================

void InteresingProcessAnaMgr::selectDecay(const G4Track* aTrack) {
    G4ParticleDefinition* particle = aTrack->GetDefinition();
    // skip optical photon
    if (particle == G4OpticalPhoton::Definition()) {
        return;
    }
    // skip positronium
    if (particle == G4Positronium::Definition()) {
        return;
    }

    // get the secondaries
    G4TrackingManager* tm = G4EventManager::GetEventManager() 
                                            -> GetTrackingManager();
    G4TrackVector* secondaries = tm->GimmeSecondaries();
    if(not secondaries) {
        return;
    }

    size_t nSeco = secondaries->size();

    for (size_t i = 0; i < nSeco; ++i) {
        G4Track* sectrk = (*secondaries)[i];
        G4ParticleDefinition* secparticle = sectrk->GetDefinition();
        // const G4String& sec_name = secparticle->GetParticleName();
        const G4VProcess* creatorProcess = sectrk->GetCreatorProcess();

        if (creatorProcess->GetProcessName()!="Decay") {
            continue;
        }

        NormalTrackInfo* info = (NormalTrackInfo*)(sectrk->GetUserInformation());

        if (not info) {
            LogWarn << "--- get Track Info failed, create new Track Info" << std::endl;
            info = new NormalTrackInfo(aTrack);
            sectrk->SetUserInformation(info);
        }

        info -> setDecayMode();


    }

}

void InteresingProcessAnaMgr::saveDecayInit(const G4Track* aTrack) {
    saveSecondaryInit(aTrack);
}

void InteresingProcessAnaMgr::saveDecayExit(const G4Track* aTrack) {
    saveSecondaryExit(aTrack);
}

void
InteresingProcessAnaMgr::saveSecondaryInit(const G4Track* aTrack) {
    NormalTrackInfo* info = (NormalTrackInfo*)(aTrack->GetUserInformation());
    if ((not info)) {
        return;
    }
    // Create SimTrack for the secondaries
    if (simtracksvc) {
        // NOTE: because pre tracking action could be invoked
        //       multiple times for the same track, we need to
        //       check whether the sim track already exists or not.
        if (simtracksvc->get(aTrack->GetTrackID())) {
            return;
        }


        info->markTracedAncestors(aTrack->GetTrackID());
        JM::SimTrack* jm_trk = new JM::SimTrack();
        jm_trk->setPDGID   (aTrack->GetDefinition()->GetPDGEncoding());
        jm_trk->setTrackID (aTrack->GetTrackID());
        jm_trk->setInitMass(aTrack->GetDefinition()->GetPDGMass());

        jm_trk->setInitPx  (aTrack->GetMomentum().x());
        jm_trk->setInitPy  (aTrack->GetMomentum().y());
        jm_trk->setInitPz  (aTrack->GetMomentum().z());
        jm_trk->setInitX   (aTrack->GetPosition().x());
        jm_trk->setInitY   (aTrack->GetPosition().y());
        jm_trk->setInitZ   (aTrack->GetPosition().z());
        jm_trk->setInitT   (aTrack->GetGlobalTime());

        // Exit info is collected in InteresingProcessAnaMgr::saveNeutronCaptureSecondary
        // Edep/Qedep is collected in DepositEnergyAnaMgr

        simtracksvc->put(jm_trk);
    }


}

void
InteresingProcessAnaMgr::saveSecondaryExit(const G4Track* aTrack) {
    NormalTrackInfo* info = (NormalTrackInfo*)(aTrack->GetUserInformation());
    if ((not info)) {
        return;
    }

    // Create SimTrack for the secondaries
    if (simtracksvc) {
        int current_track_id = aTrack->GetTrackID();

        JM::SimTrack* jm_trk = simtracksvc->get(current_track_id);
        if (!jm_trk) {
            LogWarn << "Failed to find track with trackID "
                    << current_track_id
                    << std::endl;
            return;
        }

        jm_trk->setExitPx  (aTrack->GetMomentum().x());
        jm_trk->setExitPy  (aTrack->GetMomentum().y());
        jm_trk->setExitPz  (aTrack->GetMomentum().z());
        jm_trk->setExitX   (aTrack->GetPosition().x());
        jm_trk->setExitY   (aTrack->GetPosition().y());
        jm_trk->setExitZ   (aTrack->GetPosition().z());
        jm_trk->setExitT   (aTrack->GetGlobalTime());

        jm_trk->setTrackLength(aTrack->GetTrackLength());

    }

}
