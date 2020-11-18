#ifndef DepositEnergyAnaMgr_hh
#define DepositEnergyAnaMgr_hh
/*
 * This is used for the tracking of energy deposit of specific primary track
 */
#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IAnalysisElement.h"
#include "TTree.h"
class ISimTrackSvc;

class DepositEnergyAnaMgr : public IAnalysisElement,
                            public ToolBase{
public:

    DepositEnergyAnaMgr(const std::string& name);
    ~DepositEnergyAnaMgr();
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
    double calculateQuenched(const G4Step* step);
    bool save_into_data_model();
private:
    static const int MAX_PARTICLES = 10000;

    double m_BirksConstant1;
    double m_BirksConstant2;

    // Evt Data
    TTree* m_evt_tree;
    int m_eventID;
    int m_init_nparticles; // The total number of primary tracks

    Int_t m_pdgid[MAX_PARTICLES];
    Int_t m_trkid[MAX_PARTICLES];

    Float_t m_edep[MAX_PARTICLES];
    Float_t m_edep_x[MAX_PARTICLES];
    Float_t m_edep_y[MAX_PARTICLES];
    Float_t m_edep_z[MAX_PARTICLES];

    Float_t m_q_edep[MAX_PARTICLES];
    Float_t m_q_edep_x[MAX_PARTICLES];
    Float_t m_q_edep_y[MAX_PARTICLES];
    Float_t m_q_edep_z[MAX_PARTICLES];

    // deposit data not in LS
    Float_t m_edep_notinLS[MAX_PARTICLES];

    bool m_flag_ntuple;

private:
    ISimTrackSvc* simtracksvc;
};

#endif
