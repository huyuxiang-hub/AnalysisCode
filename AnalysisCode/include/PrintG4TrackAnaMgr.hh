#ifndef PrintG4TrackAnaMgr_hh
#define PrintG4TrackAnaMgr_hh

#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IAnalysisElement.h"
#include "JunoTimer/IJunoTimerSvc.h"
#include "JunoTimer/JunoTimer.h"

class PrintG4TrackAnaMgr: public IAnalysisElement,
                   public ToolBase {
public:
    PrintG4TrackAnaMgr(const std::string& name);
    ~PrintG4TrackAnaMgr();

    // Run Action
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
    // Event Action
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

private:
    int m_verbose;
    std::vector<int> m_eventID;
    
   
};
#endif
