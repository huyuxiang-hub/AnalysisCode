#ifndef MuIsoProcessAnaMgr_hh
#define MuIsoProcessAnaMgr_hh

#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IAnalysisElement.h"
#include "TTree.h"

class MuIsoProcessAnaMgr: public IAnalysisElement,
                               public ToolBase{
public:

    MuIsoProcessAnaMgr(const std::string& name);
    ~MuIsoProcessAnaMgr();
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
    void selectMichael(const G4Track* aTrack);
    void selectNeutronCapture(const G4Track* aTrack);

private:
    void saveMichael(const G4Track* aTrack);
    void saveNeutronCapture(const G4Track* aTrack);

private:
    // Evt Data
    TTree* m_evt_tree;
    int m_eventID;
    // - michael electron
    int n_michael;
    int   michael_pdgid[100];
    float michael_kine [100];
    float michael_px   [100];
    float michael_py   [100];
    float michael_pz   [100];
    float michael_pos_x[100];
    float michael_pos_y[100];
    float michael_pos_z[100];
    float michael_t    [100];
    
    // - neutron capture in LS? TODO: Only select the specific material
    // -- detron + gamma (2.2MeV)
    // -- C13    + gamma (4.95MeV)
    // -- C13    + gamma (3.68MeV)
    //           + gamma (1.26MeV)
    TTree* m_neutron_tree;
    int n_neutron_only;    // This only contains the neutron's info
    int neutron_only_trkid[100];
    float neutron_only_kine[100];
    float neutron_only_t[100];
    float neutron_only_start_x[100];
    float neutron_only_start_y[100];
    float neutron_only_start_z[100];
    float neutron_only_stop_x[100];
    float neutron_only_stop_y[100];
    float neutron_only_stop_z[100];
    float neutron_only_track_length[100];

    int n_neutron_capture; // It will contain all secondaries
                           // please use the track id to seperate different
                           // captures (assume it will have several captures)
    int   neutron_trkid[100]; // parent track id (i.e. neutron)
    int   neutron_pdgid[100];
    float neutron_kine [100];
    float neutron_px   [100];
    float neutron_py   [100];
    float neutron_pz   [100];
    float neutron_pos_x[100];
    float neutron_pos_y[100];
    float neutron_pos_z[100];
    float neutron_t    [100];

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
    double mTrackLengthInWater[maxMuN];
    double mTrackLengthInSteel[maxMuN];
    double mTrackLengthInOil[maxMuN];
    double mTrackLengthInScint[maxMuN];
    double mELossInRock[maxMuN];
    double mELossInWater[maxMuN];
    double mELossInSteel[maxMuN];
    double mELossInOil[maxMuN];
    double mELossInScint[maxMuN];
    double mMuExitPosx[maxMuN];
    double mMuExitPosy[maxMuN];
    double mMuExitPosz[maxMuN];
    double ExitFlag[maxMuN];

    static const int maxBinN =20000;
   int      m_tbin; 
   double   m_trigT[maxBinN];
   double   m_posx[maxBinN];
   double   m_posy[maxBinN];
   double   m_posz[maxBinN];
   double   m_dE[maxBinN];
   double   m_QdE[maxBinN];
   double   m_totalQdE;
   // for stop muon information
   int      m_DecayMuCharge;
   double   m_DecayMuE;
   double   m_DecayMuP;
   double   m_DecayMuposx;
   double   m_DecayMuposy;
   double   m_DecayMuposz;
   double   m_DecayMutime;
   TTree*   m_muonIso_tree;

   static const int maxNeuN =100;
   double   m_NeuE[maxNeuN]; 
   double   m_NeuP[maxNeuN];
   double   m_Neuposx[maxNeuN];
   double   m_Neuposy[maxNeuN];
   double   m_Neuposz[maxNeuN];
   double   m_Neutime[maxNeuN];
   int      m_neuN;
  //////////////////////////////////
   int      m_maxSecondary; 
   int      newTrackId;
   int      secSize;
   int      startNum;
   std::vector<TString> m_PartName; //
   std::vector<TString> m_PartProc; //
   std::vector<TString> m_VolMaterial; //
   std::vector<int> m_ParentID; //
   std::vector<double> m_ParentKinE; //
   std::vector<TString> m_ParentPartName; //
   std::vector<int> m_ParentStepNo; //
   std::vector<double>   m_PartPosx;
   std::vector<double>   m_PartPosy;
   std::vector<double>   m_PartPosz;
   std::vector<double>   m_PartTime;
   std::vector<double>   m_PartKinE;
   std::vector<double>   m_ELoss;

   static const int maxIsoN =80000;
   int      m_IsoN;
   int      m_IsoName[maxIsoN];
   int      m_IsoProc[maxIsoN];
   double   m_IsoPosx[maxIsoN];
   double   m_IsoPosy[maxIsoN];
   double   m_IsoPosz[maxIsoN];
   double   m_IsoKinE[maxIsoN];
   double   m_IsoTime[maxIsoN] ;
 ///////////// large enegy loss point of muon
  static const int maxLosN =5000;
   int      m_LosN;
   double   m_LosPosx[maxLosN];
   double   m_LosPosy[maxLosN];
   double   m_LosPosz[maxLosN];
   double   m_LosE[maxLosN];
   double   m_LosStpL[maxLosN];
   double   m_Lostime[maxLosN];

};

#endif
