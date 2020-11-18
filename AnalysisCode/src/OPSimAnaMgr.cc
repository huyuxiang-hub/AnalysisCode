#include <boost/python.hpp>

#include "OPSimAnaMgr.hh"
//  for event
#include <sstream>
#include <cassert>
#include "junoHit_PMT.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4EventManager.hh"
#include "G4TrackingManager.hh"
#include "G4OpticalPhoton.hh"
#include "G4LossTableManager.hh"
#include "G4Electron.hh"
#include "G4SystemOfUnits.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/SniperDataPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperLog.h"
#include "RootWriter/RootWriter.h"

#include "NormalTrackInfo.hh"
#include "EvtNavigator/NavBuffer.h"
#include "BufferMemMgr/IDataMemMgr.h"
#include "Event/SimHeader.h"

#include "OPSimulator/IOPSimSvc.h"

DECLARE_TOOL(OPSimAnaMgr);

OPSimAnaMgr::OPSimAnaMgr(const std::string& name) 
    : ToolBase(name), m_opsimsvc(0)
{
    declProp("BirksConstant1", m_BirksConstant1 = 6.5e-3*g/cm2/MeV);
    declProp("BirksConstant2", m_BirksConstant2 = 1.5e-6*(g/cm2/MeV)*(g/cm2/MeV));
}

OPSimAnaMgr::~OPSimAnaMgr()
{

}

void
OPSimAnaMgr::BeginOfRunAction(const G4Run* /*aRun*/) {
    // check the RootWriter is Valid.

    SniperPtr<RootWriter> svc(*getParent(), "RootWriter");
    if (svc.invalid()) {
        LogError << "Can't Locate RootWriter. If you want to use it, please "
                 << "enalbe it in your job option file."
                 << std::endl;
        return;
    }

    SniperPtr<IOPSimSvc> opsimsvc(*getParent(), "OPSimSvc");
    if (opsimsvc.valid()) {
        LogInfo << "Get the OPSimSvc" << std::endl;
        m_opsimsvc = opsimsvc.data();
    }

}

void
OPSimAnaMgr::EndOfRunAction(const G4Run* /*aRun*/) {

}

void
OPSimAnaMgr::BeginOfEventAction(const G4Event* evt) {
    m_eventID = evt->GetEventID();
}

void
OPSimAnaMgr::EndOfEventAction(const G4Event* evt) {

    // invoke OPSimSvc
    if (m_opsimsvc) {

        opsimulator::IOPSimulator* opsim = m_opsimsvc->get_opsimulator();

        opsim->simulate();

    }
}


void
OPSimAnaMgr::PreUserTrackingAction(const G4Track* aTrack) {
}

void
OPSimAnaMgr::PostUserTrackingAction(const G4Track* aTrack) {
}

void
OPSimAnaMgr::UserSteppingAction(const G4Step* step) {
    G4Track* track = step->GetTrack();
    G4double edep = step->GetTotalEnergyDeposit();

    // invoke OPSimSvc
    if (m_opsimsvc && track->GetDefinition() != G4OpticalPhoton::Definition()) {

        opsimulator::IOPSimulator* opsim = m_opsimsvc->get_opsimulator();

        opsimulator::IGenStep* gs = opsim->create_genstep();
        gs->eventid(m_eventID);

        gs->trackid(track->GetTrackID());
        gs->pdgcode(track->GetParticleDefinition()->GetPDGEncoding());
        gs->trackst(track->GetTrackStatus());
        gs->matname(track->GetMaterial()->GetName());

        gs->stepno(track->GetCurrentStepNumber());
        gs->edep(edep);

        // calculate quenched energy
        double QuenchedTotalEnergyDeposit = 0.; 
        double dE = step->GetTotalEnergyDeposit();
        double dx = step->GetStepLength();

        if (track->GetParticleDefinition()->GetPDGEncoding()==22) { // gamma
            G4LossTableManager* manager = G4LossTableManager::Instance();
            dx = manager->GetRange(G4Electron::Electron(), dE, track->GetMaterialCutsCouple());
        }
        G4Material* aMaterial = step->GetPreStepPoint()->GetMaterial();
        G4MaterialPropertiesTable* aMaterialPropertiesTable =
            aMaterial->GetMaterialPropertiesTable();
        if (aMaterialPropertiesTable) {

            // There are some properties. Is there a scintillator property?
            const G4MaterialPropertyVector* Fast_Intensity =
                aMaterialPropertiesTable->GetProperty("FASTCOMPONENT");
            const G4MaterialPropertyVector* Slow_Intensity =
                aMaterialPropertiesTable->GetProperty("SLOWCOMPONENT");

            if (Fast_Intensity || Slow_Intensity ) {
                // It's a scintillator.
                double delta = dE/dx/aMaterial->GetDensity();
                //double birk1 = 0.0125*g/cm2/MeV;
                double birk1 = m_BirksConstant1;
                if(track->GetDefinition()->GetPDGCharge()>1.1)//for particle charge greater than 1.
                    birk1 = 0.57*birk1;
                //double birk2 = (0.0031*g/MeV/cm2)*(0.0031*g/MeV/cm2);
                double birk2 = m_BirksConstant2;
                QuenchedTotalEnergyDeposit = dE /(1+birk1*delta+birk2*delta*delta);
            }
        }
        gs->evis(QuenchedTotalEnergyDeposit);
        // end 

        gs->steplen(step->GetStepLength());

        gs->pre_x(step->GetPreStepPoint()->GetPosition().x());
        gs->pre_y(step->GetPreStepPoint()->GetPosition().y());
        gs->pre_z(step->GetPreStepPoint()->GetPosition().z());
        gs->pre_t(step->GetPreStepPoint()->GetGlobalTime());
        gs->pre_velocity(step->GetPreStepPoint()->GetVelocity());

        gs->post_x(step->GetPostStepPoint()->GetPosition().x());
        gs->post_y(step->GetPostStepPoint()->GetPosition().y());
        gs->post_z(step->GetPostStepPoint()->GetPosition().z());
        gs->post_t(step->GetPostStepPoint()->GetGlobalTime());
        gs->post_velocity(step->GetPostStepPoint()->GetVelocity());

        gs->delta_x(step->GetDeltaPosition().x());
        gs->delta_y(step->GetDeltaPosition().y());
        gs->delta_z(step->GetDeltaPosition().z());

    }


}

