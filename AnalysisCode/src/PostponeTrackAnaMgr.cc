#include "PostponeTrackAnaMgr.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperLog.h"
#include "SniperKernel/SniperException.h"

#include "GenTools/GenEventBuffer.h"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4SystemOfUnits.hh"
#include <cassert>
#include <algorithm>

DECLARE_TOOL(PostponeTrackAnaMgr);

PostponeTrackAnaMgr::PostponeTrackAnaMgr(const std::string& name) 
    : ToolBase(name)
{
    m_split_mode = kPrimaryTrack;
    declProp("SplitMode", m_split_mode_str="PrimaryTrack");

    declProp("TimeCut", m_mode_by_time_cut=3000.*ns);
    m_mode_by_time_primarytrack = 0.0;
}

PostponeTrackAnaMgr::~PostponeTrackAnaMgr() 
{

}
// Run Action
void
PostponeTrackAnaMgr::BeginOfRunAction(const G4Run*) {
    // convert the string to enum
    if (m_split_mode_str == "PrimaryTrack") {
        m_split_mode = kPrimaryTrack;
    } else if (m_split_mode_str == "EveryTrack") {
        m_split_mode = kEveryTrack;
    } else if (m_split_mode_str == "Time") {
        m_split_mode = kTime;
    } else {
        LogWarn << "Unknown mode: " << m_split_mode_str << std::endl;
        LogWarn << "Use default mode PrimaryTrack" << std::endl;
        m_split_mode = kPrimaryTrack;
    }

    LogInfo << "Track Split Mode: " << m_split_mode_str << std::endl;
}

void
PostponeTrackAnaMgr::EndOfRunAction(const G4Run*) {

}

// Event Action
void
PostponeTrackAnaMgr::BeginOfEventAction(const G4Event* anEvent) {
    m_eventid = anEvent->GetEventID();

    m_mode_by_time_primarytrack = 0.0;

    // select the minimal start time
    G4int nVertex = anEvent -> GetNumberOfPrimaryVertex();
    for (G4int index=0; index < nVertex; ++index) {
        G4PrimaryVertex* vtx = anEvent->GetPrimaryVertex( index );
        if (index==0) {
            m_mode_by_time_primarytrack = vtx->GetT0();
        } else {
            if (vtx->GetT0()<m_mode_by_time_primarytrack) {
                m_mode_by_time_primarytrack = vtx->GetT0();
            }
        }

    }
    LogDebug << "Primary Track start time: " 
             << m_mode_by_time_primarytrack
             << std::endl;
    m_cache_event = 0;
    m_cache_vertex.clear();
}


static bool cmp_vertex (HepMC::GenVertex* i, HepMC::GenVertex* j) { 
    return (i->position().t()<j->position().t()); 
}

void
PostponeTrackAnaMgr::EndOfEventAction(const G4Event*) {
    if (m_cache_vertex.size() == 0) {
        return;
    }

    // first, sort the vertex by time
    std::sort(m_cache_vertex.begin(), m_cache_vertex.end(), cmp_vertex);

    // split them
    double reftime = m_cache_vertex[0]->position().t();
    double subtime = 0;
    std::vector<HepMC::GenEvent*> evts;

    m_cache_event = new HepMC::GenEvent();
    evts.push_back(m_cache_event);


    for (std::vector<HepMC::GenVertex*>::iterator it = m_cache_vertex.begin();
            it != m_cache_vertex.end(); ++it) {

        if ( ((*it)->position().t() - reftime) >= m_mode_by_time_cut ) {
            // == flush the old one into buffer ==

            // == update reftime ==
            subtime = reftime;
            reftime = (*it)->position().t();
            // == create a new event ==
            m_cache_event = new HepMC::GenEvent();
            evts.push_back(m_cache_event);
        }

        // substract the time ref to the last chunk
        HepMC::FourVector vtx = (*it)->position();
        double t = vtx.t() - subtime;
        vtx.setT(t);
        (*it)->set_position(vtx);

        // append the vertex into event
        m_cache_event->set_signal_process_vertex((*it));
    }


    // retrieve the GenEventBuffer
    GenEventBuffer* geb = GenEventBuffer::instance();
    // increase the event 
    Task* curTask = getParent();

    LogInfo << "In the current simevent, split into " << evts.size() << " events." << std::endl;

    for (std::vector<HepMC::GenEvent*>::iterator it = evts.begin();
            it != evts.end(); ++it) {
        m_cache_event = (*it);
        m_cache_event->set_event_number(m_eventid);
        // print the vertex size
        LogDebug << "current cache vertices size: "
                << m_cache_event->vertices_size()
                << std::endl;
        if (logLevel()<=2) {
            m_cache_event->print();
        }
        GenEventPtr hepmc_ptr(m_cache_event);
        geb->push_back(hepmc_ptr);

        int evtMax = curTask->evtMax();
        if (evtMax > 0) {
            curTask->setEvtMax(evtMax+1);
        }

    }
}

// Tracking Action
void
PostponeTrackAnaMgr::PreUserTrackingAction(const G4Track* track) {

    G4Track* aTrack = const_cast<G4Track*>(track);
    // if the track already stop and killed, don't need to postpone
    if (aTrack->GetTrackStatus()==fStopAndKill) {
        return;
    }
    // if the track is optical photon, just skip.
    // that means, we only split non-OP track.
    if (aTrack->GetDefinition()==G4OpticalPhoton::Definition()) {
        return;
    }

    // TODO: Different split mode
    // * primary_track
    //   Only postpone the primary track. But sometimes the muon shower can 
    //   produce a lot of secondaries (non OP).
    // * every_track (non OP)
    // * time (non OP)
    //   
    if (m_split_mode == kPrimaryTrack) {
        return;
    } else if (m_split_mode == kEveryTrack) {
        return;
    } else if (m_split_mode == kTime) {
        double deltatime = aTrack->GetGlobalTime()-m_mode_by_time_primarytrack;
        if (deltatime < 0) {
            LogError << "The deltatime < 0.0; "
                     << "Track Time/Primary Time: "
                     << aTrack->GetGlobalTime() 
                     << " / "
                     << m_mode_by_time_primarytrack
                     << std::endl;
        }
        assert(deltatime>=0);
        if (deltatime<m_mode_by_time_cut) {
            return;
        }
    } else {
        LogWarn << "Unknown split mode or unsupport mode." << std::endl;
        return;
    }

    postpone(aTrack);

    aTrack->SetTrackStatus(fStopAndKill);
}

void
PostponeTrackAnaMgr::PostUserTrackingAction(const G4Track* aTrack) {
    if (aTrack->GetTrackStatus() == fStopAndKill) {
        // try to show the secondaries
    }
}

// Stepping Action
void
PostponeTrackAnaMgr::UserSteppingAction(const G4Step* step) {

    G4Track* aTrack = step->GetTrack();
    // if the track already stop and killed, don't need to postpone
    if (aTrack->GetTrackStatus()==fStopAndKill) {
        return;
    }
    // if the track is optical photon, just skip.
    // that means, we only split non-OP track.
    if (aTrack->GetDefinition()==G4OpticalPhoton::Definition()) {
        return;
    }

    // TODO: Different split mode
    // * primary_track
    //   Only postpone the primary track. But sometimes the muon shower can 
    //   produce a lot of secondaries (non OP).
    // * every_track (non OP)
    // * time (non OP)
    //   
    if (m_split_mode == kPrimaryTrack) {
        // now only postpone the primary track
        if (aTrack->GetParentID() != 0) {
            return;
        }
    } else if (m_split_mode == kEveryTrack) {
        // skip Optical Photon
        if (aTrack->GetDefinition()==G4OpticalPhoton::Definition()) {
            return;
        }
    } else if (m_split_mode == kTime) {
        return;
    } else {
        LogWarn << "Unknown split mode or unsupport mode." << std::endl;
        return;
    }

    postpone(aTrack);

    aTrack->SetTrackStatus(fStopAndKill);
}

// 
void
PostponeTrackAnaMgr::postpone(G4Track* aTrack) {
    // Info needed for hepmc:
    // * position
    // * time
    // * pdgid
    // * momentum
    // * mass
    const G4ThreeVector& pos = aTrack->GetPosition();
    double time = aTrack->GetGlobalTime();
    int pdg_id = aTrack->GetDefinition()->GetPDGEncoding();
    G4ThreeVector mom = aTrack->GetMomentum();
    double energy = aTrack->GetTotalEnergy();

    // == vertex
    HepMC::GenVertex* hepmc_vertex = new HepMC::GenVertex(
                           HepMC::FourVector(pos.x(), pos.y(), pos.z(), time));
    // == particle 
    HepMC::FourVector fourmom(mom.x(),mom.y(),mom.z(),energy);
    HepMC::GenParticle* hepmc_particle = new HepMC::GenParticle(fourmom, pdg_id, 1 /*  status */); 
    hepmc_vertex->add_particle_out(hepmc_particle);

    // == add vertex into cache
    m_cache_vertex.push_back(hepmc_vertex);

}
