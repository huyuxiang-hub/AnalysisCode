#ifndef PhotonTrackingAnaMgr_hh
#define PhotonTrackingAnaMgr_hh

/*
 * This analysis manager will save every step of the OPs.
 * The following info will be saved:
 *      * pre/post-step info
 *      * track id
 *      * related physics process
 */

#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IAnalysisElement.h"
#include "TTree.h"
#include "G4String.hh"
#include "G4OpProcessSubType.hh"
#include <vector>

class PhotonTrackingAnaMgr: public IAnalysisElement,
                           public ToolBase {
public:
    PhotonTrackingAnaMgr(const std::string& name);
    ~PhotonTrackingAnaMgr();

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

    TTree* m_steps;
    int evtid;
    std::vector<int> trackid;
    std::vector<int> parentid;
    // == pre ==
    // === time and position ===
    std::vector<float> pre_t;
    std::vector<float> pre_x;
    std::vector<float> pre_y;
    std::vector<float> pre_z;
    // == post ==
    // === time and position ===
    std::vector<float> post_t;
    std::vector<float> post_x;
    std::vector<float> post_y;
    std::vector<float> post_z;

    std::vector<int> post_op_subtype; // Refer to G4OpProcessSubType.hh
    std::vector<int> post_op_boundarytype; // Refer to G4OpBoundaryProcess.hh

};

#endif
