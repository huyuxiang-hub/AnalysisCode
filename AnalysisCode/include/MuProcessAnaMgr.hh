#ifndef MuProcessAnaMgr_hh
#define MuProcessAnaMgr_hh

#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IAnalysisElement.h"
#include "TTree.h"

class MuProcessAnaMgr: public IAnalysisElement,
                               public ToolBase{
public:

    MuProcessAnaMgr(const std::string& name);
    ~MuProcessAnaMgr();
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
    // Evt Data
    TTree* m_evt_tree;
    int m_eventID;

// for muon tree
    TTree*   m_muon_tree;
    bool   mMuonFlag;
    bool   idflag;
    static   const int maxMuN=25;
    int    mMuMult;
    int    maxMuMult;
    int    mPDG[maxMuN]; 
    int    mTrackId[maxMuN];
    int    mParId[maxMuN];  
    double mMuInitPosx[maxMuN];
    double mMuInitPosy[maxMuN];
    double mMuInitPosz[maxMuN] ;
    double mMuInitPx[maxMuN];
    double mMuInitPy[maxMuN];
    double mMuInitPz[maxMuN];
    double mMuInitKine[maxMuN];
    double mTrackLengthInRock[maxMuN];
    double mTrackLengthInTT[maxMuN];
    double mTrackLengthInVetoWater[maxMuN];
    double mTrackLengthInCDWater[maxMuN];
    double mTrackLengthInAcrylic[maxMuN];
    double mTrackLengthInSteel[maxMuN];
    double mTrackLengthInOil[maxMuN];
    double mTrackLengthInScint[maxMuN];
    double mELossInRock[maxMuN];
    double mELossInTT[maxMuN];
    double mELossInVetoWater[maxMuN];
    double mELossInCDWater[maxMuN];
    double mELossInAcrylic[maxMuN];
    double mELossInSteel[maxMuN];
    double mELossInOil[maxMuN];
    double mELossInScint[maxMuN];
    double mMuExitPosx[maxMuN];
    double mMuExitPosy[maxMuN];
    double mMuExitPosz[maxMuN];
    double ExitFlag[maxMuN];
    int    mustpMaterial[maxMuN];
    double mustpx[maxMuN];
    double mustpy[maxMuN];
    double mustpz[maxMuN];
};

#endif
