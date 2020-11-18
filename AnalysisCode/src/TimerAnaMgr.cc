#include "TimerAnaMgr.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/SniperDataPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperLog.h"
#include "SniperKernel/SniperException.h"

DECLARE_TOOL(TimerAnaMgr);

TimerAnaMgr::TimerAnaMgr(const std::string& name)
    : ToolBase(name), m_junotimersvc(0)
{

}

TimerAnaMgr::~TimerAnaMgr()
{

}

// ==========================================================================
// Run Action
// ==========================================================================
void
TimerAnaMgr::BeginOfRunAction(const G4Run* /*aRun*/) {
    SniperPtr<IJunoTimerSvc> timersvc(*getParent(), "JunoTimerSvc");
    if (timersvc.invalid()) {
        LogError << "Failed to get JunoTimerSvc!" << std::endl;
        throw SniperException("Make sure you have load the JunoTimerSvc.");
    }
    m_junotimersvc = timersvc.data();

    // create timers
    m_timer_event = m_junotimersvc->get("SIMEVT/timerevt");
    m_timer_allevent = m_junotimersvc->get("SIMEVT/timerallevt");
}

void
TimerAnaMgr::EndOfRunAction(const G4Run* /*aRun*/) {
    LogInfo << "summaries: " << std::endl;
    LogInfo << "number of measurements: " 
            << m_timer_event->number_of_measurements() << std::endl;
    LogInfo << "mean time: " << m_timer_event->mean() << " ms" << std::endl;

    m_timer_allevent->stop();
    LogInfo << "total time elapsed: " 
            << m_timer_allevent->elapsed() << " ms" << std::endl;
}

// ==========================================================================
// Event Action
// ==========================================================================
void
TimerAnaMgr::BeginOfEventAction(const G4Event* /*aEvent*/) {
    m_timer_event->start();
    m_timer_allevent->resume();
}

void
TimerAnaMgr::EndOfEventAction(const G4Event* /*aEvent*/) {
    m_timer_event->stop();
    m_timer_allevent->pause();
}
