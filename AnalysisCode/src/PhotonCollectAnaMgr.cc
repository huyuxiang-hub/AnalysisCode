#include "PhotonCollectAnaMgr.hh"
#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/SniperDataPtr.h"
#include "SniperKernel/SniperLog.h"
#include "SniperKernel/ToolFactory.h"
#include "RootWriter/RootWriter.h"

#include "G4Event.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"

DECLARE_TOOL(PhotonCollectAnaMgr);

PhotonCollectAnaMgr::PhotonCollectAnaMgr(const std::string& name)
    : ToolBase(name)
{

}

PhotonCollectAnaMgr::~PhotonCollectAnaMgr() 
{

}

void
PhotonCollectAnaMgr::BeginOfRunAction(const G4Run* /*run*/)
{
    SniperPtr<RootWriter> svc(*getParent(), "RootWriter");
    if (svc.invalid()) {
        LogError << "Can't Locate RootWriter. If you want to use it, please "
                 << "enalbe it in your job option file."
                 << std::endl;
        return;
    }

    op_cols = svc->bookTree("SIMEVT/opticalphoton", "collection of OPs");
    op_cols->Branch("evtID", &evtid, "evtID/I");
    op_cols->Branch("t", &t, "t/F");
    op_cols->Branch("x", &x, "x/F");
    op_cols->Branch("y", &y, "y/F");
    op_cols->Branch("z", &z, "z/F");

    op_cols->Branch("px", &px, "px/F");
    op_cols->Branch("py", &py, "py/F");
    op_cols->Branch("pz", &pz, "pz/F");

    op_cols->Branch("polx", &polx, "polx/F");
    op_cols->Branch("poly", &poly, "poly/F");
    op_cols->Branch("polz", &polz, "polz/F");
}

void
PhotonCollectAnaMgr::EndOfRunAction(const G4Run* /*run*/)
{

}

void
PhotonCollectAnaMgr::BeginOfEventAction(const G4Event* event)
{
    evtid = event->GetEventID();
    m_proc2cnt.clear();
}

void
PhotonCollectAnaMgr::EndOfEventAction(const G4Event* /*event*/)
{
    for (std::map<G4String, int>::iterator it = m_proc2cnt.begin();
            it != m_proc2cnt.end(); ++it) {
        LogDebug << it->first << " " << it->second << std::endl;
    }
}

void
PhotonCollectAnaMgr::PreUserTrackingAction(const G4Track* aTrack)
{
    if (aTrack->GetDefinition() != G4OpticalPhoton::Definition()) {
        return;
    }
    G4Track* track = const_cast<G4Track*>(aTrack);
    // get the process
    const G4VProcess* proc = track->GetCreatorProcess();
    const G4String& procname = proc->GetProcessName();
    m_proc2cnt[procname] += 1;
    // save the info first
    t = track->GetGlobalTime();
    const G4ThreeVector& pos = track->GetPosition();
    x = pos.x();
    y = pos.y();
    z = pos.z();
    G4ThreeVector mom = track->GetMomentum();
    px = mom.x();
    py = mom.y();
    pz = mom.z();
    const G4ThreeVector& pol = track->GetPolarization();
    polx = pol.x();
    poly = pol.y();
    polz = pol.z();
    op_cols->Fill();

    // kill the track
    track->SetTrackStatus(fStopAndKill);
}

void
PhotonCollectAnaMgr::PostUserTrackingAction(const G4Track* /*aTrack*/)
{

}

void
PhotonCollectAnaMgr::UserSteppingAction(const G4Step* /*step*/)
{

}
