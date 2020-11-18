#ifndef DepositEnergyCalibAnaMgr_hh
#define DepositEnergyCalibAnaMgr_hh
/*
 * This is used for the tracking of energy deposit of specific primary track
 */
#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IAnalysisElement.h"
#include "TTree.h"

class DepositEnergyCalibAnaMgr : public IAnalysisElement,
                            public ToolBase{
public:

    DepositEnergyCalibAnaMgr(const std::string& name);
    ~DepositEnergyCalibAnaMgr();
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

    double m_BirksConstant1;
    double m_BirksConstant2;

    // Evt Data
    TTree* m_evt_tree;
    int m_eventID;
    int m_init_nparticles; // The total number of primary tracks

    Int_t m_pdgid[100];
    Int_t m_trkid[100];

    Float_t m_edep[100];
    Float_t m_edep_x[100];
    Float_t m_edep_y[100];
    Float_t m_edep_z[100];

    Float_t m_q_edep[100];
    Float_t m_q_edep_x[100];
    Float_t m_q_edep_y[100];
    Float_t m_q_edep_z[100];

    // deposit data not in LS
    Float_t m_edep_notinLS[100];

    // deposit energy in Different Material
    // - total deposit energy in the LS, window (Mylar) and the pipe (acrylic)
    Float_t m_edep_LS[100];
    Float_t m_edep_Mylar[100];
    Float_t m_edep_Acrylic[100];

    // - total travel length of the particle in the LS, window (Mylar) and the pipe (acrylic)
    Float_t m_travel_len_LS[100];
    Float_t m_travel_len_Mylar[100];
    Float_t m_travel_len_Acrylic[100];

    // - total deposit energy of the gammas (from the e+ annihilation) in the LS, window, and pipe.
    Float_t m_edep_gamma_LS[100];
    Float_t m_edep_gamma_Mylar[100];
    Float_t m_edep_gamma_Acrylic[100];

    // - the start position where the particle hits the window
    Float_t m_hit_wind_pos_x[100];
    Float_t m_hit_wind_pos_y[100];
    Float_t m_hit_wind_pos_z[100];

    // - then end position where the particle stops (for positron, it is where it annihilates)

    bool m_flag_ntuple;
};

#endif
