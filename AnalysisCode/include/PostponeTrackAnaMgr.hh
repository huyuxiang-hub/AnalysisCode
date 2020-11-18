#ifndef PostponeTrackAnaMgr_hh
#define PostponeTrackAnaMgr_hh
/*
 * This analysis manager will postpone track(s) to next event.
 * It need to work together with GenEventBuffer in GenTools.
 */

#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IAnalysisElement.h"
#include <vector>

namespace HepMC {
    class GenEvent;
    class GenVertex;
};

enum TrackSplitMode {
    kPrimaryTrack = 0,
    kEveryTrack = 1,
    kTime = 2,
    kUnknown = 0xFFFF
};

class PostponeTrackAnaMgr: public IAnalysisElement,
                           public ToolBase {
public:
    PostponeTrackAnaMgr(const std::string& name);
    ~PostponeTrackAnaMgr();

    // Run Action
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
    // Event Action
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);
    // Tracking Action
    virtual void PreUserTrackingAction(const G4Track* aTrack);
    virtual void PostUserTrackingAction(const G4Track* aTrack);
    // Stepping Action
    virtual void UserSteppingAction(const G4Step* step);

private:
    void postpone(G4Track*);

private:
    int m_eventid;
    std::string m_split_mode_str;
    TrackSplitMode m_split_mode;

    double m_mode_by_time_primarytrack;
    double m_mode_by_time_cut;

    std::vector<HepMC::GenVertex*> m_cache_vertex;
    HepMC::GenEvent* m_cache_event;
};

#endif
