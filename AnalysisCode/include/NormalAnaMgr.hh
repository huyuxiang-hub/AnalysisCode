#ifndef NormalAnaMgr_hh
#define NormalAnaMgr_hh

#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IAnalysisElement.h"
#include "TTree.h"
#include "TH1I.h"
#include <map>

class NormalAnaMgr: public IAnalysisElement,
                    public ToolBase{
public:

    NormalAnaMgr(const std::string& name);
    ~NormalAnaMgr();
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

    // Evt Data
    TTree* m_evt_tree;
    Int_t m_eventID;
    Int_t m_nPhotons;
    Int_t m_totalPE;
   
    std::vector<int>	 m_nPE ;             
    std::vector<float>   m_energy;
    std::vector<double>  m_hitTime;
    std::vector<int>     m_pmtID;
    std::vector<int>  	 m_peTrackID;
    std::vector<int>  	 m_isCerenkov;
    std::vector<int>     m_isReemission;
    std::vector<int>     m_isOriginalOP;
    std::vector<double>  m_OriginalOPTime;
    
    // PMT Info
    //
    Int_t m_npmts_byPMT;
   
    std::vector<int>  m_nPE_byPMT;
    std::vector<int>  m_PMTID_byPMT;
    
    
    std::map<int, int> m_cache_bypmt;
   
    // - 2015.10.10 Tao Lin
    //   The hit's local position in PMT will be saved.
    //   However, to save the space, Float is enough.
    std::vector<float>  m_localpos_x;  
    std::vector<float>  m_localpos_y;
    std::vector<float>  m_localpos_z;
  

    
    // - 2016.04.17 Tao Lin
    //   hit's local direction in PMT
   
    std::vector<float>  m_localdir_x;   
    std::vector<float>  m_localdir_y;
    std::vector<float>  m_localdir_z;



    // - 2017.03.01 Tao Lin
    //   hit's global position in PMT
    std::vector<float> m_globalpos_x  ;  
    std::vector<float> m_globalpos_y  ;
    std::vector<float> m_globalpos_z  ;
                      
    std::vector<float> m_boundarypos_x;
    std::vector<float> m_boundarypos_y;
    std::vector<float> m_boundarypos_z;
    

    Float_t m_edep;
    Float_t m_edep_x;
    Float_t m_edep_y;
    Float_t m_edep_z;
 
    bool m_flag_ntuple;
    bool m_flag_hitinfo;
    TH1I* m_step_no;
};

#endif
