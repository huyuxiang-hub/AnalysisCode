#ifndef GenEvtInfoAnaMgr_hh
#define GenEvtInfoAnaMgr_hh

#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IAnalysisElement.h"
#include "TTree.h"

/*
 * This Analysis Element is used for saving the initial particle's info
 *
 * At EndOfEventAction, we can get the particle and its track ID in geant4.
 */

class GenEvtInfoAnaMgr: public IAnalysisElement,
                        public ToolBase{
public:

    GenEvtInfoAnaMgr(const std::string& name);
    ~GenEvtInfoAnaMgr();
    // Run Action
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
    // Event Action
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);
    // Tracking Action
    virtual void PreUserTrackingAction(const G4Track* aTrack);
    virtual void PostUserTrackingAction(const G4Track* aTrack);

private:
    bool check_init_same();
    bool save_into_data_model();
private:

    static const int MAX_PARTICLES = 10000;

    TTree* m_evt_tree;
    int m_eventID;
    int m_init_nparticles;
    int m_init_pdgid[MAX_PARTICLES];
    int m_init_trkid[MAX_PARTICLES];
  float m_init_x[MAX_PARTICLES];
  float m_init_y[MAX_PARTICLES];
  float m_init_z[MAX_PARTICLES];
  float m_init_px[MAX_PARTICLES];
  float m_init_py[MAX_PARTICLES];
  float m_init_pz[MAX_PARTICLES];
  float m_init_mass[MAX_PARTICLES];
 double m_init_time[MAX_PARTICLES];
  // exit
  float m_exit_x[MAX_PARTICLES];
  float m_exit_y[MAX_PARTICLES];
  float m_exit_z[MAX_PARTICLES];
 double m_exit_t[MAX_PARTICLES];
  float m_exit_px[MAX_PARTICLES];
  float m_exit_py[MAX_PARTICLES];
  float m_exit_pz[MAX_PARTICLES];
  // track length
  float m_track_length[MAX_PARTICLES]; // please note, the length is not only
                             // in the LS.

private:
  // helper: select info by track id
  bool m_by_track_enable[MAX_PARTICLES];
  float m_by_track_init_x[MAX_PARTICLES];
  float m_by_track_init_y[MAX_PARTICLES];
  float m_by_track_init_z[MAX_PARTICLES];
 double m_by_track_init_t[MAX_PARTICLES];

  float m_by_track_exit_x[MAX_PARTICLES];
  float m_by_track_exit_y[MAX_PARTICLES];
  float m_by_track_exit_z[MAX_PARTICLES];
 double m_by_track_exit_t[MAX_PARTICLES];

  float m_by_track_length[MAX_PARTICLES];

  float m_by_track_exit_px[MAX_PARTICLES];
  float m_by_track_exit_py[MAX_PARTICLES];
  float m_by_track_exit_pz[MAX_PARTICLES];

    bool m_flag_ntuple;
};

#endif

