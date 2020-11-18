#ifndef DepositEnergyTTAnaMgr_hh
#define DepositEnergyTTAnaMgr_hh
/*
 * This is used for the tracking of energy deposit of specific primary track
 */
#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IAnalysisElement.h"
#include "TTree.h"
#include <vector>
#include "TRandom3.h"
#include "DetSimAlg/IDetElement.h"
#include <string>
class DepositEnergyTTAnaMgr : public IAnalysisElement,
                            public ToolBase{
public:

    DepositEnergyTTAnaMgr(const std::string& name);
    ~DepositEnergyTTAnaMgr();
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
    bool save_into_data_model();
private:


    void ComputeDigitTT();

    // Evt Data
    TTree* m_evt_tree;
    TTree* m_evt_treeTT;
    TTree* m_evt_treeTTDigit;

    int m_eventID;
    int m_NDeposits;
    std::vector<int>   m_dep_trackID;
    std::vector<int>   m_dep_pdg;
    std::vector<float> m_dep_E;
    std::vector<float> m_dep_E0;
    std::vector<int>   m_dep_nbar;
    std::vector<int>   m_dep_npanel;
    std::vector<int>   m_dep_nplane;
    std::vector<int>   m_dep_nwall;
    std::vector<float> m_dep_x;
    std::vector<float> m_dep_y;
    std::vector<float> m_dep_z;
    std::vector<float> m_dep_t;
    std::vector<float> m_dep_dL;
    std::vector<float> m_dep_dR;
    std::vector<float> m_dep_peL;
    std::vector<float> m_dep_peR;
    std::vector<float> m_dep_tL;
    std::vector<float> m_dep_tR;

    std::vector<int> m_isMuonDeposits;
    
    float coeff[4][16];
  
    TRandom3 *m_rnd;
    IDetElement* de;
         

    int m_TB_channel[31744];
    int m_TB_DMchannel[63488];

 
    int m_NTouchedBar;
    float m_TB_peL[31744];
    float m_TB_peR[31744];
    float m_TB_timeL[31744];
    float m_TB_timeR[31744];
    int m_TB_ADCL[31744];
    int m_TB_ADCR[31744];
    float m_TB_xc[31744];
    float m_TB_yc[31744];
    float m_TB_zc[31744];
    int m_TB_is_ct[31744];
    int m_TB_isMuonDeposits[31744];
 
  int m_NTouchedChannel;
  float m_TB_pe[63488];
  float m_TB_time[63488];
  int m_TB_ADC[63488];
  int m_TB_channelC[63488];
  
  float m_TB_xcC[63488];
  float m_TB_ycC[63488];
  float m_TB_zcC[63488];
  int m_TB_is_ctC[63488];
  int m_TB_isMuonDepositsC[63488];

  int CTchannels[64][21];
                                
                                
    bool m_flag_ntuple;
};

#endif

/* vim: set noexpandtab shiftwidth=2 tabstop=8 softtabstop=2: */
