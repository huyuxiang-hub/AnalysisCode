#include "PrintTrackInfoAnaMgr.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/SniperDataPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperLog.h"
#include "SniperKernel/SniperException.h"

DECLARE_TOOL(PrintTrackInfoAnaMgr);

PrintTrackInfoAnaMgr::PrintTrackInfoAnaMgr(const std::string& name)
    : ToolBase(name)
{

}

PrintTrackInfoAnaMgr::~PrintTrackInfoAnaMgr()
{

}

// ==========================================================================
// Run Action
// ==========================================================================
void
PrintTrackInfoAnaMgr::BeginOfRunAction(const G4Run* /*aRun*/) {
}

void
PrintTrackInfoAnaMgr::EndOfRunAction(const G4Run* /*aRun*/) {
}

// ==========================================================================
// Event Action
// ==========================================================================
void
PrintTrackInfoAnaMgr::BeginOfEventAction(const G4Event* /*aEvent*/) {
}

void
PrintTrackInfoAnaMgr::EndOfEventAction(const G4Event* /*aEvent*/) {
}
