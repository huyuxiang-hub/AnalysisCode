#include "MuonToySim.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperLog.h"
#include "SniperKernel/SniperException.h"

#include "DataRegistritionSvc/DataRegistritionSvc.h"

#include "BufferMemMgr/IDataMemMgr.h"
#include "EvtNavigator/NavBuffer.h"
#include "Event/SimHeader.h"
#include "Event/SimEvent.h"
#include "G4SystemOfUnits.hh"
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandExponential.h"
#include "CLHEP/Random/RandFlat.h"

DECLARE_TOOL(MuonToySim);

    MuonToySim::MuonToySim(const std::string& name)
: ToolBase(name)
{
    photonyield      = 10596/MeV;   //photon/MeV
    absorptlength    = 30*m;   //m
    pmt_area         = 3.14*250*250*mm*mm; //mm^2 
    effeciency       = 0.35; 
    clight           = 299.792458*mm/ns ; 
    nLS              = 1.485; 
    yieldratio      = 0.799; 
    fasttimeconstant = 4.93*ns; 
    slowtimeconstant = 20.6*ns; 
    m_pmt_total      = 17746; 

    m_pmtpos.reserve(m_pmt_total); 

}

MuonToySim::~MuonToySim()
{

}

#include <cstdlib>
void
MuonToySim::BeginOfRunAction(const G4Run* /*aRun*/) {
    // SniperPtr<DataRegistritionSvc> drsSvc(getScope(), "DataRegistritionSvc");
    // if ( drsSvc.invalid() ) {
    //     LogError << "Failed to get DataRegistritionSvc!" << std::endl;
    //     throw SniperException("Make sure you have load the DataRegistritionSvc.");
    // }
    // drsSvc->registerData("JM::SimEvent", "/Event/Sim");
    
    std::string base = getenv("RECCDMUONALGROOT"); 
    std::string pmtfile = base+"/share/pmtinfo.txt" ; 
    std::ifstream inf(pmtfile.c_str()); 
    double x, y, z, pmtid; 
    for(int ipmt=0; ipmt<m_pmt_total; ipmt++){
        inf >>pmtid >>x >>y >>z; 
        m_pmtpos.push_back(G4ThreeVector(x, y, z)); 
        /*
        G4cout << m_pmtpos[ipmt].x() << " " 
              << m_pmtpos[ipmt].y() << " " 
              << m_pmtpos[ipmt].z() << G4endl;  
              */
    }
}

void
MuonToySim::EndOfRunAction(const G4Run* /*aRun*/) {

}

void
MuonToySim::BeginOfEventAction(const G4Event* /*evt*/) {

    // create header
    sim_header = new JM::SimHeader;
    // create event
    sim_event = new JM::SimEvent(0);
}

void
MuonToySim::EndOfEventAction(const G4Event* evt) {
    // FIXME: shall we get the navigator first?
    JM::EvtNavigator* nav = new JM::EvtNavigator();
    TTimeStamp ts;
    nav->setTimeStamp(ts);

    SniperPtr<IDataMemMgr> mMgr(*getParent(), "BufferMemMgr");
    mMgr->adopt(nav, "/Event");
    // create header
    //JM::SimHeader* sim_header = new JM::SimHeader;
    // create event
    //JM::SimEvent* sim_event = new JM::SimEvent(0);
    // == fill hits
    //fill_hits(sim_event, evt);
    // == fill tracks
    fill_tracks(sim_event, evt);

    // set the relation
    sim_header->setEvent(sim_event);
    nav->addHeader("/Event/Sim", sim_header);

}

#include "G4Event.hh"
#include "G4Track.hh"
#include "G4SDManager.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"

#include "junoHit_PMT.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"

G4double 
MuonToySim::pmtResponse(G4double cos_theta){
    G4double response; 
    response  = 0.378+0.5093*cos_theta+0.1135*cos_theta*cos_theta; 
    return response; 
}
void 
MuonToySim::UserSteppingAction(const G4Step* astep){

    G4Track* atrack = astep->GetTrack(); 
    G4double edep = astep->GetTotalEnergyDeposit();     
    //G4cout << "edep:" << edep << G4endl; 
    //G4cout << "particle:" << atrack->GetDefinition()->GetParticleName() << G4endl; 
    if( (edep>0) && (atrack->GetDefinition()->GetParticleName() == "mu-")  &&  (atrack->GetMaterial()->GetName() == "LS")){

        G4StepPoint* prePoint = astep->GetPreStepPoint(); 
        G4StepPoint* postPoint = astep->GetPostStepPoint(); 
        G4ThreeVector prePos = prePoint->GetPosition(); 
        G4ThreeVector postPos = postPoint->GetPosition(); 
        G4cout << "STEPSTUDY: (" << postPos.x()/mm  << "mm\t" << postPos.y()/mm << "mm\t" << postPos.z()/mm  << "mm)\t" << edep/MeV << "MeV"<< G4endl; 
        G4double preTime = prePoint->GetGlobalTime(); 
        G4double postTime = postPoint->GetGlobalTime(); 
        G4ThreeVector astepCen = (postPoint->GetPosition()+prePoint->GetPosition())/2; 
        G4double normargu = edep*photonyield*pmt_area/(4*3.1415926); 

        for(int ipmt=0; ipmt<m_pmt_total; ipmt++){
            G4double hittime; 
            G4ThreeVector Rd = m_pmtpos[ipmt]-astepCen; 
            G4double costhetad = Rd.cosTheta(m_pmtpos[ipmt]); 
            G4double exp_npe = normargu*pmtResponse(costhetad)/Rd.mag2()*std::exp(-Rd.mag()/absorptlength)*effeciency; 

            G4int  npe = CLHEP::RandPoisson::shoot(exp_npe); 
            for(G4int ipe=0; ipe<npe; ipe++)
            {

                G4double rand = CLHEP::RandFlat::shoot(); 
                G4ThreeVector emitPos = prePos + (postPos-prePos)*rand; 
                hittime = preTime + (postTime-preTime)*rand + (m_pmtpos[ipmt]-emitPos).mag()*nLS/clight; 
                JM::SimPMTHit* jm_hit = sim_event->addCDHit();
                jm_hit->setPMTID( ipmt );
                jm_hit->setHitTime( hittime/ns );
                jm_hit->setNPE( 1 );
                //G4double timeconst = CLHEP::RandFlat::shoot()>yieldratio?slowtimeconstant:fasttimeconstant;  
                //G4double hittime_real = hittime + CLHEP::RandExponential::shoot(timeconst); 
                //jm_hit->setHitTime( hittime_real/ns );
                /*
                   G4cout << "pmt:" << ipmt << G4endl
                   << "Pos:" << Rd/cm  << "cm"<< G4endl
                   << "hittime" << hittime/ns << "ns" << G4endl
                   << "npe" << npe << G4endl; 
                   */
            }
        }
    }
}
    void 
MuonToySim::fill_hits(JM::SimEvent* dst, const G4Event* evt)
{

    LogDebug << "Begin Fill Hits" << std::endl;
    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    G4int CollID = SDman->GetCollectionID("hitCollection");

    G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
    junoHit_PMT_Collection* col = (junoHit_PMT_Collection*)(HCE->GetHC(CollID));

    // fill evt data
    // int totPE = 0;
    if (col) {
        int n_hit = col->entries();
        // m_nPhotons = n_hit;

        for (int i = 0; i < n_hit; ++i) {
            // create new hit
            JM::SimPMTHit* jm_hit = dst->addCDHit();
            jm_hit->setPMTID( (*col)[i]->GetPMTID() );
            jm_hit->setNPE( (*col)[i]->GetCount() );
            jm_hit->setHitTime( (*col)[i]->GetTime() );
        }

    }
    LogDebug << "End Fill Hits" << std::endl;

}

    void 
MuonToySim::fill_tracks(JM::SimEvent* dst, const G4Event* evt)
{
    LogDebug << "Begin Fill Tracks" << std::endl;

    G4int nVertex = evt -> GetNumberOfPrimaryVertex();
    for (G4int index=0; index < nVertex; ++index) {
        G4PrimaryVertex* vtx = evt->GetPrimaryVertex( index );

        // Vertex info
        double x = vtx->GetX0();
        double y = vtx->GetY0();
        double z = vtx->GetZ0();
        double t = vtx->GetT0();

        // Loop Over Particle
        G4PrimaryParticle* pp = vtx -> GetPrimary();

        while (pp) {

            int trkid = pp -> GetTrackID();
            int pdgid = pp -> GetPDGcode();
            double px = pp -> GetPx();
            double py = pp -> GetPy();
            double pz = pp -> GetPz();
            double mass = pp -> GetMass();

            // new track
            JM::SimTrack* jm_trk = dst->addTrack();
            jm_trk->setPDGID(pdgid);
            jm_trk->setTrackID(trkid);
            jm_trk->setInitPx(px);
            jm_trk->setInitPy(py);
            jm_trk->setInitPz(pz);
            jm_trk->setInitMass(mass);
            jm_trk->setInitX(x);
            jm_trk->setInitY(y);
            jm_trk->setInitZ(z);
            jm_trk->setInitT(t);

            pp = pp->GetNext();
        }
    }
    LogDebug << "End Fill Tracks" << std::endl;
}
