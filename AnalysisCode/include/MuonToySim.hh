#ifndef MuonToySim_hh
#define MuonToySim_hh

#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IAnalysisElement.h"
#include <G4ThreeVector.hh>

class G4Event;
namespace JM {
    class SimEvent;
    class SimHeader; 
}

class MuonToySim: public IAnalysisElement,
                       public ToolBase{

public:
    MuonToySim(const std::string& name);
    ~MuonToySim();

    // Run Action
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
    // Event Action
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);
    // Stepping Action
    virtual void UserSteppingAction(const G4Step*); 
private:

    // fill hits
    void fill_hits(JM::SimEvent* dst, const G4Event* evt);
    // fill tracks
    void fill_tracks(JM::SimEvent* dst, const G4Event* evt);

    G4double pmtResponse(G4double ); 

    int m_pmt_total;
    std::vector<G4ThreeVector> m_pmtpos; 

    G4double photonyield    ;   //photon/MeV
    G4double absorptlength  ;   //m
    G4double pmt_area       ;   //mm^2 
    G4double effeciency     ; 
    G4double clight         ;   //speed of light
    G4double nLS            ;   // refrection effeciency of LS
    G4double yieldratio     ;   //ratio of fast component 
    G4double fasttimeconstant;  //the time constant of fast component
    G4double slowtimeconstant;  //the time constant of slow component

    JM::SimHeader* sim_header; 
    JM::SimEvent* sim_event; 
                       };

#endif
