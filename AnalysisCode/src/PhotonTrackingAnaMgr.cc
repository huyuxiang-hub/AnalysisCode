#include "PhotonTrackingAnaMgr.hh"
#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/SniperDataPtr.h"
#include "SniperKernel/SniperLog.h"
#include "SniperKernel/ToolFactory.h"
#include "RootWriter/RootWriter.h"

// force load the dict for vector<float>
#include "TROOT.h"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"
#include "G4SteppingManager.hh"
#include "G4TrackingManager.hh"
#include "G4EventManager.hh"
#include "G4OpBoundaryProcess.hh"


DECLARE_TOOL(PhotonTrackingAnaMgr);

PhotonTrackingAnaMgr::PhotonTrackingAnaMgr(const std::string& name)
    : ToolBase(name), m_steps(0)
{

}

PhotonTrackingAnaMgr::~PhotonTrackingAnaMgr()
{

}

void
PhotonTrackingAnaMgr::BeginOfRunAction(const G4Run* /*run*/)
{
    SniperPtr<RootWriter> svc(*getParent(), "RootWriter");
    if (svc.invalid()) {
        LogError << "Can't Locate RootWriter. If you want to use it, please "
                 << "enalbe it in your job option file."
                 << std::endl;
        return;
    }

    gROOT->ProcessLine("#include <vector>");

    m_steps = svc->bookTree("SIMEVT/opsteps", "collection of OP steps");
    m_steps->Branch("evtID", &evtid, "evtID/I");
    m_steps->Branch("TrackID",  &trackid);
    m_steps->Branch("ParentID", &parentid);

    m_steps->Branch("PreT", &pre_t);
    m_steps->Branch("PreX", &pre_x);
    m_steps->Branch("PreY", &pre_y);
    m_steps->Branch("PreZ", &pre_z);

    m_steps->Branch("PostT", &post_t);
    m_steps->Branch("PostX", &post_x);
    m_steps->Branch("PostY", &post_y);
    m_steps->Branch("PostZ", &post_z);

    m_steps->Branch("OPType", &post_op_subtype);
    m_steps->Branch("BoundaryType", &post_op_boundarytype);
}

void
PhotonTrackingAnaMgr::EndOfRunAction(const G4Run* /*run*/)
{

}

void
PhotonTrackingAnaMgr::BeginOfEventAction(const G4Event* event)
{
    evtid = event->GetEventID(); 

    trackid.clear();
    parentid.clear();

    pre_t.clear();
    pre_x.clear();
    pre_y.clear();
    pre_z.clear();

    post_t.clear();
    post_x.clear();
    post_y.clear();
    post_z.clear();

    post_op_subtype.clear();

}

void
PhotonTrackingAnaMgr::EndOfEventAction(const G4Event* /*event*/)
{
    if (not m_steps) {
        return;
    }

    LogDebug << "summary:" << std::endl;
    LogDebug << " pre_t.size(): " << pre_t.size() << std::endl;
    LogDebug << " pre_x.size(): " << pre_x.size() << std::endl;
    LogDebug << " pre_y.size(): " << pre_y.size() << std::endl;
    LogDebug << " pre_z.size(): " << pre_z.size() << std::endl;

    m_steps->Fill();
}

void
PhotonTrackingAnaMgr::PreUserTrackingAction(const G4Track* /*aTrack*/)
{

}

void
PhotonTrackingAnaMgr::PostUserTrackingAction(const G4Track* /*aTrack*/)
{

}

void
PhotonTrackingAnaMgr::UserSteppingAction(const G4Step* step)
{
    G4Track* track = step->GetTrack();
    // only save the OP
    if (track->GetDefinition() != G4OpticalPhoton::Definition()) {
        return;
    }
    G4StepPoint* prepoint = step->GetPreStepPoint();
    G4StepPoint* postpoint = step->GetPostStepPoint();

    trackid.push_back( track->GetTrackID() );
    parentid.push_back( track->GetParentID() );

    // == pre ==
    const G4ThreeVector& tmp_pre_pos = prepoint->GetPosition();
    float tmp_pre_t = prepoint->GetGlobalTime();

    pre_t.push_back(tmp_pre_t);
    pre_x.push_back(tmp_pre_pos.x());
    pre_y.push_back(tmp_pre_pos.y());
    pre_z.push_back(tmp_pre_pos.z());

    // == post ==
    const G4ThreeVector& tmp_post_pos = postpoint->GetPosition();
    float tmp_post_t = postpoint->GetGlobalTime();

    post_t.push_back(tmp_post_t);
    post_x.push_back(tmp_post_pos.x());
    post_y.push_back(tmp_post_pos.y());
    post_z.push_back(tmp_post_pos.z());

    int tmp_post_op_subtype;
    const G4VProcess* tmp_post_proc = postpoint->GetProcessDefinedStep();
    if (tmp_post_proc) {
        tmp_post_op_subtype = tmp_post_proc->GetProcessSubType();
    } else {
        tmp_post_op_subtype = -1;
    }
    // boundary process
    G4OpBoundaryProcessStatus boundaryStatus=Undefined;
    static G4OpBoundaryProcess* boundary=0;

    //find the boundary process only once
    if(!boundary){
        G4ProcessManager* pm 
            = track->GetDefinition()->GetProcessManager();
        G4int nprocesses = pm->GetProcessListLength();
        G4ProcessVector* pv = pm->GetProcessList();
        G4int i;
        for( i=0;i<nprocesses;i++){
            if((*pv)[i]->GetProcessName()=="OpBoundary"){
                boundary = (G4OpBoundaryProcess*)(*pv)[i];
                boundary->SetVerboseLevel(2);
                break;
            }
        }
    }
    if(postpoint->GetStepStatus()==fGeomBoundary) {
        boundaryStatus = boundary->GetStatus();
    }

    LogDebug << "@@@ op_subtype/boundary_status: " 
             << tmp_post_op_subtype << " / "
             << boundaryStatus
             << " @@@ " << std::endl;

    post_op_subtype.push_back(tmp_post_op_subtype);
    post_op_boundarytype.push_back(boundaryStatus);
}
