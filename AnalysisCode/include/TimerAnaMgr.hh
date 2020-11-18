#ifndef TimerAnaMgr_hh
#define TimerAnaMgr_hh

#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IAnalysisElement.h"
#include "JunoTimer/IJunoTimerSvc.h"
#include "JunoTimer/JunoTimer.h"

class TimerAnaMgr: public IAnalysisElement,
                   public ToolBase {
public:
    TimerAnaMgr(const std::string& name);
    ~TimerAnaMgr();

    // Run Action
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
    // Event Action
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

private:
    IJunoTimerSvc* m_junotimersvc;

    JunoTimerPtr m_timer_event;
    JunoTimerPtr m_timer_allevent;
};
#endif
