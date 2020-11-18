#ifndef RadioAnaMgr_hh
#define RadioAnaMgr_hh

/*
 * RadioAnaMgr: handle Geant4 Radioactivity Decay Module related.
 *
 * -- Tao Lin <lintao@ihep.ac.cn>, 20 May 2020
 */


#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IAnalysisElement.h"
#include "TTree.h"
#include "TH1I.h"
#include <map>
#include <vector>
#include <string>

class RadioAnaMgr: public IAnalysisElement,
                   public ToolBase{
public:

    RadioAnaMgr(const std::string& name);
    ~RadioAnaMgr();
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
    // StopAtPa234m: Stop simulation from Pa234m to Pa234.
    //               As Geant4 radioactivity decay will still simulate it,
    //               we just kill the secondaries of decay (Pa234m->Pa234)
    bool m_StopAtPa234m;

private:
    // User data
    TTree* m_radio_tree;
    // Data Structure
    //   Instead of saving all tracks, we only store the radioactivity decay 
    //   including the parent and its decay products.
    //
    //   radioidx is used to identify the decay.
    //
    //   radioidx trkid parentid pdgcode name(str) t ...
    //          0   1      0             U238
    //          0   2      1             Th234
    //          0   3      1             alpha
    //
    // As all the excited states will be simulated, the decay of excited states
    // will be also saved into the data.
    //
    // The time is the generation time of daughter.
    //
    int m_evtid;
    std::vector<int> m_radioidx;
    std::vector<int> m_trkid;
    std::vector<int> m_parentid;
    std::vector<int> m_pdgid;
    std::vector<std::string> m_name;
    std::vector<double> m_time;
    std::vector<double> m_Ek; // Kinetic Energy
private:
    // helper
    int current_radioidx;
};

#endif

