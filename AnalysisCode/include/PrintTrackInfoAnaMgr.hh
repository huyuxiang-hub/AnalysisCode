#ifndef PrintTrackInfoAnaMgr_hh
#define PrintTrackInfoAnaMgr_hh

#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IAnalysisElement.h"
#include "JunoTimer/IJunoTimerSvc.h"
#include "JunoTimer/JunoTimer.h"

class PrintTrackInfoAnaMgr: public IAnalysisElement,
                   public ToolBase {
public:
    PrintTrackInfoAnaMgr(const std::string& name);
    ~PrintTrackInfoAnaMgr();

    // Run Action
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
    // Event Action
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

private:
    
};
#endif
