#ifndef PhotonCollectAnaMgr_hh
#define PhotonCollectAnaMgr_hh

/*
 * This analysis manager will collect the optical photons generated from 
 * scintillation.
 * * save the info of the OP.
 * * kill the tracking of this OP in geant4.
 */

#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IAnalysisElement.h"
#include "TTree.h"
#include "G4String.hh"
#include <map>

class PhotonCollectAnaMgr: public IAnalysisElement,
                           public ToolBase {
public:
    PhotonCollectAnaMgr(const std::string& name);
    ~PhotonCollectAnaMgr();

    // Run Action
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
    // Event Action
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

    virtual void PreUserTrackingAction(const G4Track* aTrack);
    virtual void PostUserTrackingAction(const G4Track* aTrack);

    virtual void UserSteppingAction(const G4Step* step);

private:
    TTree* op_cols;

    int evtid;
    // time and position
    float t;
    float x;
    float y;
    float z;
    // momentum
    float px;
    float py;
    float pz;
    // polarization
    float polx;
    float poly;
    float polz;

    std::map<G4String, int> m_proc2cnt;
};

#endif
