#include "PrintTrackInfoAnaMgr.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/SniperDataPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperLog.h"
#include "SniperKernel/SniperException.h"
#include "G4UImanager.hh"

#include "G4Event.hh"


DECLARE_TOOL(PrintTrackInfoAnaMgr);

PrintTrackInfoAnaMgr::PrintTrackInfoAnaMgr(const std::string& name)
    : ToolBase(name)
{
  declProp("VerBose", m_verbose = 0);
  declProp("EventID", m_eventID);

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
PrintTrackInfoAnaMgr::BeginOfEventAction(const G4Event* evt /*aEvent*/) {

  std::cout<<"PrintTrackInfoAnaMgr!!!!!!!!!!"<<std::endl;
  if( m_verbose > 2 or m_verbose < 0 ){
   
    LogInfo<<"verbose has a invalid value !!"<<std::endl;
    return;
  }

  std::vector<int>::iterator ret;
  ret = std::find(m_eventID.begin(), m_eventID.end(), evt->GetEventID());
     
  if (ret!=m_eventID.end())
        {
          G4UImanager * UImanager = G4UImanager::GetUIpointer();
          if(m_verbose == 0 ){
             UImanager->ApplyCommand("/tracking/verbose 0");
           } 
          else if ( m_verbose == 1 ){
              UImanager->ApplyCommand("/tracking/verbose 1");
            }
          else {
              UImanager->ApplyCommand("/tracking/verbose 2");
             }
          
        }
  else
        {
          G4UImanager * UImanager = G4UImanager::GetUIpointer();
          UImanager->ApplyCommand("/tracking/verbose 0");
        }


}

void
PrintTrackInfoAnaMgr::EndOfEventAction(const G4Event* /*aEvent*/) {
}
