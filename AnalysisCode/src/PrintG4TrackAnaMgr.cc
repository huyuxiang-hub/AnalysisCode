#include "PrintG4TrackAnaMgr.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/SniperDataPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperLog.h"
#include "SniperKernel/SniperException.h"
#include "G4UImanager.hh"

#include "G4Event.hh"


DECLARE_TOOL(PrintG4TrackAnaMgr);

PrintG4TrackAnaMgr::PrintG4TrackAnaMgr(const std::string& name)
    : ToolBase(name)
{
  declProp("VerBose", m_verbose = 0);
  declProp("EventID", m_eventID);

}

PrintG4TrackAnaMgr::~PrintG4TrackAnaMgr()
{

}

// ==========================================================================
// Run Action
// ==========================================================================
void
PrintG4TrackAnaMgr::BeginOfRunAction(const G4Run* /*aRun*/) {
}

void
PrintG4TrackAnaMgr::EndOfRunAction(const G4Run* /*aRun*/) {
}

// ==========================================================================
// Event Action
// ==========================================================================
void
PrintG4TrackAnaMgr::BeginOfEventAction(const G4Event* evt /*aEvent*/) {

  std::cout<<"PrintG4TrackAnaMgr!!!!!!!!!!"<<std::endl;
  if( m_verbose > 5 or m_verbose < 0 ){
   
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
          else if ( m_verbose == 2 ){
              UImanager->ApplyCommand("/tracking/verbose 2");
            }
          else if ( m_verbose == 3 ){
              UImanager->ApplyCommand("/tracking/verbose 3");
            }
          else if ( m_verbose == 4 ){
              UImanager->ApplyCommand("/tracking/verbose 4");
            }
          else if ( m_verbose == 5 ){
              UImanager->ApplyCommand("/tracking/verbose 5");
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
PrintG4TrackAnaMgr::EndOfEventAction(const G4Event* /*aEvent*/) {
}
