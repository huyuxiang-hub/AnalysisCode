#include <boost/python.hpp>
#include <iostream>

#include "GenEvtInfoAnaMgr.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/SniperDataPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperLog.h"
#include "RootWriter/RootWriter.h"

#include "G4Event.hh"
#include "G4Track.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"

#include "EvtNavigator/NavBuffer.h"
#include "BufferMemMgr/IDataMemMgr.h"
#include "Event/SimHeader.h"

DECLARE_TOOL(GenEvtInfoAnaMgr);

GenEvtInfoAnaMgr::GenEvtInfoAnaMgr(const std::string& name)
    : ToolBase(name) {
    declProp("EnableNtuple", m_flag_ntuple=true);
    m_evt_tree = 0;

}

GenEvtInfoAnaMgr::~GenEvtInfoAnaMgr() {

}

void 
GenEvtInfoAnaMgr::BeginOfRunAction(const G4Run* /*aRun*/) {
    if (not m_flag_ntuple) {
        return;
    }
    //---------------//
    SniperPtr<ISimTrackSvc> simtracksvc_ptr(getParent(), "SimTrackSvc");
    if (simtracksvc_ptr.invalid()) {
        simtracksvc = dynamic_cast<ISimTrackSvc*>(getParent()->createSvc("SimTrackSvc"));
    } else {
        simtracksvc = simtracksvc_ptr.data();
    }

    // check the RootWriter is Valid.

    SniperPtr<RootWriter> svc(*getParent(), "RootWriter");
    if (svc.invalid()) {
        LogError << "Can't Locate RootWriter. If you want to use it, please "
                 << "enalbe it in your job option file."
                 << std::endl;
        return;
    }
    m_evt_tree = svc->bookTree("SIMEVT/geninfo", "Generator Info");

    m_evt_tree->Branch("evtID", &m_eventID, "evtID/I");
    m_evt_tree->Branch("nInitParticles", &m_init_nparticles, "nInitParticles/I");
    m_evt_tree->Branch("InitPDGID", m_init_pdgid, "InitPDGID[nInitParticles]/I");
    m_evt_tree->Branch("InitTRKID", m_init_trkid, "InitTRKID[nInitParticles]/I");
    m_evt_tree->Branch("InitX", m_init_x, "InitX[nInitParticles]/F");
    m_evt_tree->Branch("InitY", m_init_y, "InitY[nInitParticles]/F");
    m_evt_tree->Branch("InitZ", m_init_z, "InitZ[nInitParticles]/F");
    m_evt_tree->Branch("InitPX", m_init_px, "InitPX[nInitParticles]/F");
    m_evt_tree->Branch("InitPY", m_init_py, "InitPY[nInitParticles]/F");
    m_evt_tree->Branch("InitPZ", m_init_pz, "InitPZ[nInitParticles]/F");
    m_evt_tree->Branch("InitMass", m_init_mass, "InitMass[nInitParticles]/F");
    m_evt_tree->Branch("InitTime", m_init_time, "InitTime[nInitParticles]/D");

    m_evt_tree->Branch("ExitX",  m_exit_x,  "ExitX[nInitParticles]/F");
    m_evt_tree->Branch("ExitY",  m_exit_y,  "ExitY[nInitParticles]/F");
    m_evt_tree->Branch("ExitZ",  m_exit_z,  "ExitZ[nInitParticles]/F");
    m_evt_tree->Branch("ExitT",  m_exit_t,  "ExitT[nInitParticles]/D");
    m_evt_tree->Branch("ExitPX", m_exit_px, "ExitPX[nInitParticles]/F");
    m_evt_tree->Branch("ExitPY", m_exit_py, "ExitPY[nInitParticles]/F");
    m_evt_tree->Branch("ExitPZ", m_exit_pz, "ExitPZ[nInitParticles]/F");

    m_evt_tree->Branch("TrackLength", m_track_length, "TrackLength[nInitParticles]/F");
}

void 
GenEvtInfoAnaMgr::EndOfRunAction(const G4Run* /*aRun*/) {

}

void
GenEvtInfoAnaMgr::BeginOfEventAction(const G4Event* evt) {
    m_eventID = evt->GetEventID();
    m_init_nparticles = 0;
    for (int i = 0; i < MAX_PARTICLES; ++i) {
        m_init_pdgid[i] = 0;
        m_init_trkid[i] = 0;
        m_init_x[i] = 0.0;
        m_init_y[i] = 0.0;
        m_init_z[i] = 0.0;
        m_init_px[i] = 0.0;
        m_init_py[i] = 0.0;
        m_init_pz[i] = 0.0;
        m_init_mass[i] = 0.0;
        m_init_time[i] = 0.0;

        m_exit_x[i] = 0.0;
        m_exit_y[i] = 0.0;
        m_exit_z[i] = 0.0;
        m_exit_t[i] = 0.0;
        m_exit_px[i] = 0.0;
        m_exit_py[i] = 0.0;
        m_exit_pz[i] = 0.0;

        m_track_length[i] = 0.0;
    }

    for (int i = 0; i < MAX_PARTICLES; ++i) {
        m_by_track_enable[i] = false;
        m_by_track_init_x[i] = 0;
        m_by_track_init_y[i] = 0;
        m_by_track_init_z[i] = 0;
        m_by_track_init_t[i] = 0;

        m_by_track_exit_x[i] = 0;
        m_by_track_exit_y[i] = 0;
        m_by_track_exit_z[i] = 0;
        m_by_track_exit_t[i] = 0;

        m_by_track_length[i] = 0;
    }

}

void
GenEvtInfoAnaMgr::EndOfEventAction(const G4Event* evt) {
    // save the primary vertex info
    LogDebug << " + Begin Event Action" << std::endl;
    G4int nVertex = evt -> GetNumberOfPrimaryVertex();
    for (G4int index=0; index < nVertex; ++index) {
        // Loop Over Particle
        G4PrimaryVertex* vtx = evt->GetPrimaryVertex( index );

        G4PrimaryParticle* pp = vtx -> GetPrimary();

        LogDebug << " ++ Vertex Info: " << std::endl;
        LogDebug << " ++ Position: " << vtx->GetX0() << ", "
                                      << vtx->GetY0() << ", "
                                      << vtx->GetZ0() << std::endl;

        while (pp) {

            int trkid = pp -> GetTrackID();
            int pdgid = pp -> GetPDGcode();
            double px = pp -> GetPx();
            double py = pp -> GetPy();
            double pz = pp -> GetPz();

            LogDebug << " +++ " << std::endl;
            LogDebug << " ++++ TrackID : " << trkid << std::endl;
            LogDebug << " ++++ PDGID   : " << pdgid << std::endl;
            LogDebug << " ++++ Momentum: " << px << ", " << py << ", " << pz << std::endl;
            // Fill the Tree
            //
            if (m_init_nparticles < MAX_PARTICLES) {

                m_init_pdgid[m_init_nparticles] = pdgid;
                m_init_trkid[m_init_nparticles] = trkid;

                m_init_x[m_init_nparticles] = vtx->GetX0();
                m_init_y[m_init_nparticles] = vtx->GetY0();
                m_init_z[m_init_nparticles] = vtx->GetZ0();
                m_init_time[m_init_nparticles] = vtx -> GetT0();

                m_init_px[m_init_nparticles] = px;
                m_init_py[m_init_nparticles] = py;
                m_init_pz[m_init_nparticles] = pz;

                m_init_mass[m_init_nparticles] = pp -> GetMass();
                // check the init value is same
                if (not check_init_same()) {
                    LogWarn << "Check Init Info Consitent failed." << std::endl;
                }

                // if the particle can't be tracking, skip this particle.
                // TODO: do we need set the exit values as the same as the init
                // ones?
                if (trkid < 0) {
                    LogWarn << "Skip the un-tracking particle " << pdgid << std::endl;
                } else {
                    m_exit_x[m_init_nparticles] = m_by_track_exit_x[trkid];
                    m_exit_y[m_init_nparticles] = m_by_track_exit_y[trkid];
                    m_exit_z[m_init_nparticles] = m_by_track_exit_z[trkid];
                    m_exit_t[m_init_nparticles] = m_by_track_exit_t[trkid];

                    m_exit_px[m_init_nparticles] = m_by_track_exit_px[trkid];
                    m_exit_py[m_init_nparticles] = m_by_track_exit_py[trkid];
                    m_exit_pz[m_init_nparticles] = m_by_track_exit_pz[trkid];

                    m_track_length[m_init_nparticles] = m_by_track_length[trkid];
                }

                ++m_init_nparticles;
            }


            pp = pp->GetNext();
        }

    }

    if (m_flag_ntuple and m_evt_tree) {
        m_evt_tree -> Fill();
    }
    save_into_data_model();
}

void
GenEvtInfoAnaMgr::PreUserTrackingAction(const G4Track* /*aTrack*/) {

}

void
GenEvtInfoAnaMgr::PostUserTrackingAction(const G4Track* aTrack) {
    // only select the primary track
    if (aTrack->GetParentID() != 0) {
        return;
    }
    G4int current_trkid = aTrack->GetTrackID();

    // make sure the trkid is valid.
    assert(0 < current_trkid);
    assert(current_trkid < MAX_PARTICLES);
    if (current_trkid <= 0 or current_trkid >= MAX_PARTICLES) {
        return;
    }

    m_by_track_enable[current_trkid] = true;
    const G4ThreeVector& start_pos = aTrack->GetVertexPosition();
    const G4ThreeVector& stop_pos  = aTrack->GetPosition();
    m_by_track_init_x[current_trkid] = start_pos.x();
    m_by_track_init_y[current_trkid] = start_pos.y();
    m_by_track_init_z[current_trkid] = start_pos.z();
    m_by_track_init_t[current_trkid] = aTrack->GetGlobalTime()
                                       - aTrack->GetLocalTime();
    m_by_track_exit_x[current_trkid] = stop_pos.x();
    m_by_track_exit_y[current_trkid] = stop_pos.y();
    m_by_track_exit_z[current_trkid] = stop_pos.z();
    m_by_track_exit_t[current_trkid] = aTrack->GetGlobalTime();

    m_by_track_length[current_trkid] = aTrack->GetTrackLength();
    
    const G4ThreeVector& stop_mom = aTrack->GetMomentum();

    m_by_track_exit_px[current_trkid] = stop_mom.x();
    m_by_track_exit_py[current_trkid] = stop_mom.y();
    m_by_track_exit_pz[current_trkid] = stop_mom.z();

    LogDebug << " + Primary Track [" << current_trkid  << "]"<< std::endl;
    LogDebug << " +-- Init Pos: " << start_pos.x() << " "
                                  << start_pos.y() << " "
                                  << start_pos.z() << std::endl;
    LogDebug << " +-- Init Tim: " << m_by_track_init_t[current_trkid] << std::endl;
    LogDebug << " +-- Exit Pos: " << stop_pos.x() << " "
                                  << stop_pos.y() << " "
                                  << stop_pos.z() << std::endl;
    LogDebug << " +-- Exit Tim: " << m_by_track_exit_t[current_trkid] << std::endl;
    LogDebug << " +-- Exit Mom: " << stop_mom.x() << " "
                                  << stop_mom.y() << " "
                                  << stop_mom.z() << std::endl;
        
    JM::SimTrack* trk=simtracksvc->get(current_trkid);
        if(!trk){
                 trk = new JM::SimTrack();
                 trk->setTrackID(current_trkid);
                 simtracksvc->put(trk);
               }
        trk->setExitPx( m_by_track_exit_px[current_trkid]);
        trk->setExitPy( m_by_track_exit_py[current_trkid]);
        trk->setExitPz( m_by_track_exit_pz[current_trkid]);
        trk->setExitX(m_by_track_exit_x[current_trkid]);
        trk->setExitY(m_by_track_exit_y[current_trkid]);
        trk->setExitZ(m_by_track_exit_z[current_trkid]);
        trk->setExitT(m_by_track_exit_t[current_trkid]);

        trk->setTrackLength(m_by_track_length[current_trkid]);

}

bool
GenEvtInfoAnaMgr::check_init_same() {
    bool flag = true;
    G4int current_trkid = m_init_trkid[m_init_nparticles];
    // skip the un-tracking particles
    if (current_trkid < 0) {
        return false;
    }
    // same init x,y,z,t
    if ( m_init_x[m_init_nparticles] != m_by_track_init_x[current_trkid]) {
        LogWarn << "Init X are not consistent: "
                << m_init_x[m_init_nparticles] << " / "
                << m_by_track_init_x[current_trkid] 
                << std::endl;
        flag = false;
    }
    if ( m_init_y[m_init_nparticles] != m_by_track_init_y[current_trkid]) {
        LogWarn << "Init Y are not consistent: "
                << m_init_y[m_init_nparticles] << " / "
                << m_by_track_init_y[current_trkid] 
                << std::endl;
        flag = false;
    }
    if ( m_init_z[m_init_nparticles] != m_by_track_init_z[current_trkid]) {
        LogWarn << "Init Z are not consistent: "
                << m_init_z[m_init_nparticles] << " / "
                << m_by_track_init_z[current_trkid] 
                << std::endl;
        flag = false;
    }
    if ( m_init_time[m_init_nparticles] != m_by_track_init_t[current_trkid]) {
        LogWarn << "Init Z are not consistent: "
                << m_init_z[m_init_nparticles] << " / "
                << m_by_track_init_z[current_trkid] 
                << std::endl;
        flag = false;
    }

    return flag;
}

bool
GenEvtInfoAnaMgr::save_into_data_model() {
    SniperDataPtr<JM::NavBuffer>  navBuf(*getParent(), "/Event");
    if (navBuf.invalid()) {
        return false;
    }
    LogDebug << "navBuf: " << navBuf.data() << std::endl;
    JM::EvtNavigator* evt_nav = navBuf->curEvt();
    LogDebug << "evt_nav: " << evt_nav << std::endl;
    if (not evt_nav) {
        return false;
    }
   
    const std::vector<std::string>& m_paths= evt_nav->getPath();
    std::vector<std::string>::const_iterator pos = std::find(m_paths.begin(), m_paths.end(), "/Event/Sim");
    if ( m_paths.end() == pos ) {
        return false;
    }   


    JM::SimHeader* m_simheader = dynamic_cast<JM::SimHeader*>(evt_nav->getHeader("/Event/Sim"));
    LogDebug << "simhdr: " << m_simheader << std::endl;
    if (not m_simheader) {
        return false;
    }
    JM::SimEvent* m_simevent = dynamic_cast<JM::SimEvent*>(m_simheader->event());
    LogDebug << "simevt: " << m_simevent << std::endl;
    if (not m_simevent) {
        return false;
    }

    // Fill the Exit Track Info
    for (int i = 0; i < m_init_nparticles; ++i) {
        int trkid = m_init_trkid[i];
        if (trkid <= 0) {continue;}
        JM::SimTrack* jm_trk = m_simevent->findTrackByTrkID(trkid);

        // make sure the init is same ???
        jm_trk->setExitPx( m_exit_px[i] );
        jm_trk->setExitPy( m_exit_py[i] );
        jm_trk->setExitPz( m_exit_pz[i] );
        jm_trk->setExitX( m_exit_x[i] );
        jm_trk->setExitY( m_exit_y[i] );
        jm_trk->setExitZ( m_exit_z[i] );
        jm_trk->setExitT( m_exit_t[i] );

        jm_trk->setTrackLength( m_track_length[i] );

    }
    return true;
}
