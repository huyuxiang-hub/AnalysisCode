
#include "MuProcessAnaMgr.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperLog.h"
#include "RootWriter/RootWriter.h"

#include "G4SDManager.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4EventManager.hh"
#include "G4TrackingManager.hh"
#include "G4SteppingManager.hh"
#include "G4VProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "NormalTrackInfo.hh"

DECLARE_TOOL(MuProcessAnaMgr);

MuProcessAnaMgr::MuProcessAnaMgr(const std::string& name)
    : ToolBase(name)
{

}

MuProcessAnaMgr::~MuProcessAnaMgr()
{

}

// == Run Action
void
MuProcessAnaMgr::BeginOfRunAction(const G4Run* /*aRun*/) {
    SniperPtr<RootWriter> svc(*getParent(), "RootWriter");
    m_muon_tree = svc->bookTree("SIMEVT/mu", "muon events");
    m_muon_tree -> Branch("evtID", &m_eventID, "evtID/I");
    m_muon_tree->Branch("MuMult", &maxMuMult,"MuMult/I");
    m_muon_tree->Branch("PDG", mPDG,"PDG[MuMult]/I");
    m_muon_tree->Branch("MuInitPosx", mMuInitPosx,"MuInitPosx[MuMult]/D");
    m_muon_tree->Branch("MuInitPosy", mMuInitPosy,"MuInitPosy[MuMult]/D");
    m_muon_tree->Branch("MuInitPosz", mMuInitPosz,"MuInitPosz[MuMult]/D");
    m_muon_tree->Branch("MuInitPx", mMuInitPx,"MuInitPx[MuMult]/D");
    m_muon_tree->Branch("MuInitPy", mMuInitPy,"MuInitPy[MuMult]/D");
    m_muon_tree->Branch("MuInitPz", mMuInitPz,"MuInitPz[MuMult]/D");
    m_muon_tree->Branch("MuInitKine", mMuInitKine,"MuInitKine[MuMult]/D");
    m_muon_tree->Branch("TrackLengthInRock", mTrackLengthInRock,"TrackLengthInRock[MuMult]/D");
    m_muon_tree->Branch("TrackLengthInTT", mTrackLengthInTT,"TrackLengthInTT[MuMult]/D");
    m_muon_tree->Branch("TrackLengthInVetoWater", mTrackLengthInVetoWater,"TrackLengthInVetoWater[MuMult]/D");
    m_muon_tree->Branch("TrackLengthInCDWater", mTrackLengthInCDWater,"TrackLengthInCDWater[MuMult]/D");
    m_muon_tree->Branch("TrackLengthInAcrylic", mTrackLengthInAcrylic,"TrackLengthInAcrylic[MuMult]/D");
    m_muon_tree->Branch("TrackLengthInSteel", mTrackLengthInSteel,"TrackLengthInSteel[MuMult]/D");
    m_muon_tree->Branch("TrackLengthInScint", mTrackLengthInScint,"TrackLengthInScint[MuMult]/D");
    m_muon_tree->Branch("ELossInRock", mELossInRock,"ELossInRock[MuMult]/D");
    m_muon_tree->Branch("ELossInTT", mELossInTT,"ELossInTT[MuMult]/D");
    m_muon_tree->Branch("ELossInVetoWater", mELossInVetoWater,"ELossInVetoWater[MuMult]/D");
    m_muon_tree->Branch("ELossInCDWater", mELossInCDWater,"ELossInCDWater[MuMult]/D");
    m_muon_tree->Branch("ELossInAcrylic", mELossInAcrylic,"ELossInAcrylic[MuMult]/D");
    m_muon_tree->Branch("ELossInSteel", mELossInSteel,"LossInSteel[MuMult]/D");
    m_muon_tree->Branch("ELossInScint", mELossInScint,"ELossInScint[MuMult]/D");
    m_muon_tree->Branch("MuExitPosx", mMuExitPosx,"MuExitPosx[MuMult]/D");
    m_muon_tree->Branch("MuExitPosy", mMuExitPosy,"MuExitPosy[MuMult]/D");
    m_muon_tree->Branch("MuExitPosz", mMuExitPosz,"MuExitPosz[MuMult]/D");
    m_muon_tree->Branch("mustpMaterial", mustpMaterial,"mustpMaterial[MuMult]/I"); 
    m_muon_tree->Branch("mustpx", mustpx,"mustpx[MuMult]/D");
    m_muon_tree->Branch("mustpy", mustpy,"mustpy[MuMult]/D");
    m_muon_tree->Branch("mustpz", mustpz,"mustpz[MuMult]/D");
}

void
MuProcessAnaMgr::EndOfRunAction(const G4Run* /*aRun*/) {

}

// == Event Action
void
MuProcessAnaMgr::BeginOfEventAction(const G4Event* evt) {
    m_eventID = evt->GetEventID();

  mMuonFlag = false;
//  idflag = false;
  mMuMult = 0;
  maxMuMult = 0;
  for(G4int i=0; i<maxMuN; i++){
    mPDG[i] = 0;
    mTrackId[i] = 0;
    mParId[i] = 0;
    mMuInitPosx[i] =0;
    mMuInitPosy[i] =0;
    mMuInitPosz[i] =0;
    mMuInitPx[i] =0;
    mMuInitPy[i] =0;
    mMuInitPz[i] =0;
    mMuInitKine[i]=0;
    mTrackLengthInRock[i] = 0.;
    mTrackLengthInTT[i] = 0.;
    mTrackLengthInVetoWater[i] = 0.;
    mTrackLengthInCDWater[i] = 0.;
    mTrackLengthInAcrylic[i] = 0.;
    mTrackLengthInSteel[i] = 0.;
    mTrackLengthInScint[i] = 0.;
    mELossInRock[i] = 0.;
    mELossInTT[i] = 0.;
    mELossInVetoWater[i] = 0.;
    mELossInCDWater[i] = 0.;
    mELossInAcrylic[i] = 0.;
    mELossInSteel[i] = 0.;
    mELossInScint[i] = 0.;
    mMuExitPosx[i] =0;
    mMuExitPosy[i] =0;
    mMuExitPosz[i] =0;
    ExitFlag[i]=0; 
    mustpMaterial[i]=0;
    mustpx[i]=0;
    mustpy[i]=0;
    mustpz[i]=0;
 }

}

void
MuProcessAnaMgr::EndOfEventAction(const G4Event* /*evt*/) {

    m_muon_tree->Fill();
}

// == Tracking Action
void
MuProcessAnaMgr::PreUserTrackingAction(const G4Track* /*aTrack*/) {
    // NormalTrackInfo* info = (NormalTrackInfo*)(aTrack->GetUserInformation());
}

void
MuProcessAnaMgr::PostUserTrackingAction(const G4Track* /*aTrack*/) {
}

// == Stepping Action
void
MuProcessAnaMgr::UserSteppingAction(const G4Step* fStep) {
    // cache the material pointer, avoid comparison using string.
    static G4Material* s_mat_LS = nullptr;
    static G4Material* s_mat_Rock = nullptr;
    static G4Material* s_mat_Scintillator = nullptr;
    static G4Material* s_mat_vetoWater = nullptr;
    static G4Material* s_mat_Water = nullptr;
    static G4Material* s_mat_Acrylic = nullptr;
    static G4Material* s_mat_Steel = nullptr;
    if (!s_mat_LS) {
        s_mat_LS = G4Material::GetMaterial("LS");
        s_mat_Rock = G4Material::GetMaterial("Rock");
        s_mat_Scintillator = G4Material::GetMaterial("Scintillator");
        s_mat_vetoWater = G4Material::GetMaterial("vetoWater");
        s_mat_Water = G4Material::GetMaterial("Water");
        s_mat_Acrylic = G4Material::GetMaterial("Acrylic");
        s_mat_Steel = G4Material::GetMaterial("Steel");
    }


  G4Track* fTrack = fStep->GetTrack();
  G4int StepNo = fTrack->GetCurrentStepNumber();

  G4int trackId = fTrack ->GetTrackID();
  G4StepPoint* thePrePoint;
  //G4StepPoint* thePostPoint;
  G4Material* material =  fTrack->GetMaterial();
  //G4Material* material_post;
  G4ParticleDefinition* particle = fTrack->GetDefinition();
  // G4String pname = particle->GetParticleName();

  thePrePoint = fStep->GetPreStepPoint();
  //thePostPoint = fStep->GetPostStepPoint();
  material = thePrePoint->GetMaterial();
  //material_post = thePostPoint->GetMaterial();//WARNING, for the end step of particle, there is no postpoint only prepoint.
  //G4double TotalEnergyDeposit = fStep->GetTotalEnergyDeposit();
  //const G4DynamicParticle* aParticle = fTrack->GetDynamicParticle();
  //G4int parentId = fTrack->GetParentID();
  const CLHEP::Hep3Vector point = fTrack->GetPosition();
  //const G4TrackVector* fSecondary = fpSteppingManager->GetSecondary();
  /// for muon information ...
  if(trackId==1 && (particle == G4MuonMinus::Definition()||particle == G4MuonPlus::Definition())){
    mMuonFlag = true; 
  }

  G4int parId = fTrack->GetParentID();
  if(parId ==0 && (particle == G4MuonMinus::Definition() || particle == G4MuonPlus::Definition())){
   // std::cout<<"parId ="<<parId<<" ; pname ="<<pname<<" ; trackId ="<<trackId<<"; mMuMult="<<mMuMult<<std::endl;
    mMuMult = trackId - 1;
    if(mMuMult>=maxMuN){
      mMuMult = maxMuN - 1;
      G4cout<<" WARNING: maxMuN = "<<maxMuN<<" is not enough!! Sub Muon number = "<<mMuMult<<G4endl;
    }
    mTrackId[mMuMult]= trackId;
    mParId[mMuMult]  = parId;
    mPDG[mMuMult]    = particle->GetPDGEncoding();
    if(StepNo==1)   { 
         mMuInitPosx[mMuMult] =thePrePoint->GetPosition().x()/mm;  
         mMuInitPosy[mMuMult] =thePrePoint->GetPosition().y()/mm;
         mMuInitPosz[mMuMult] =thePrePoint->GetPosition().z()/mm;
         mMuInitPx[mMuMult] =thePrePoint->GetMomentum().x();
         mMuInitPy[mMuMult] =thePrePoint->GetMomentum().y();
         mMuInitPz[mMuMult] =thePrePoint->GetMomentum().z();
         mMuInitKine[mMuMult] =thePrePoint->GetKineticEnergy()/MeV;
     }
    if(material==s_mat_Rock){
      mTrackLengthInRock[mMuMult] += fStep->GetStepLength();
      mELossInRock[mMuMult] += fStep->GetDeltaEnergy();
      if(mTrackLengthInVetoWater[mMuMult]>0&&ExitFlag[mMuMult]==0)
        {
          mMuExitPosx[mMuMult] = point.x()/mm;
          mMuExitPosy[mMuMult] = point.y()/mm;
          mMuExitPosz[mMuMult] = point.z()/mm;
          ExitFlag[mMuMult] =1;
        }
           
     } else if(material==s_mat_Scintillator){
      mTrackLengthInTT[mMuMult] += fStep->GetStepLength();
      mELossInTT[mMuMult] += fStep->GetDeltaEnergy();
    } else if (material==s_mat_vetoWater){
      mTrackLengthInVetoWater[mMuMult] += fStep->GetStepLength();
      mELossInVetoWater[mMuMult] += fStep->GetDeltaEnergy();
    } else if (material==s_mat_Water){
      mTrackLengthInCDWater[mMuMult] += fStep->GetStepLength();
      mELossInCDWater[mMuMult] += fStep->GetDeltaEnergy();
    } else if (material==s_mat_Acrylic){
      mTrackLengthInAcrylic[mMuMult] += fStep->GetStepLength();
      mELossInAcrylic[mMuMult] += fStep->GetDeltaEnergy();
    } else if (material==s_mat_Steel){
      mTrackLengthInSteel[mMuMult] += fStep->GetStepLength();
      mELossInSteel[mMuMult] += fStep->GetDeltaEnergy();
    } else if (material==s_mat_LS){
      mTrackLengthInScint[mMuMult] += fStep->GetStepLength();
      mELossInScint[mMuMult] += fStep->GetDeltaEnergy();
     }
  // for stop muon
   G4TrackingManager* tm = G4EventManager::GetEventManager()     -> GetTrackingManager();
   G4SteppingManager* sm = tm->GetSteppingManager();
   //G4Step* stp = sm->GetStep();
   const G4TrackVector* secs = sm->GetSecondary();
   bool isMuStop = false;
   for (size_t i = 0; i < secs->size(); ++i) {
      G4Track* sectrk = (*secs)[i];
      const G4VProcess* creatorProcess = sectrk->GetCreatorProcess();
      if(creatorProcess->GetProcessName() =="muMinusCaptureAtRest" || creatorProcess->GetProcessName() == "Decay"){
           isMuStop = true;
           break;
     }
   }
   if(isMuStop == true){
     if(material==s_mat_LS)            mustpMaterial[mMuMult] = 1;
     else if(material==s_mat_Acrylic)  mustpMaterial[mMuMult] = 2;
     else if(material==s_mat_Water)    mustpMaterial[mMuMult] = 3; // inner water
     else if(material==s_mat_Steel)    mustpMaterial[mMuMult] = 4;//support structure
     else if(material==s_mat_vetoWater)mustpMaterial[mMuMult] = 5; // outer veto water
     else if(material==s_mat_Rock)     mustpMaterial[mMuMult] = 6;
     else mustpMaterial[mMuMult] = -1;
     mustpx[mMuMult] = point.x()/mm;
     mustpy[mMuMult] = point.y()/mm;
     mustpz[mMuMult] = point.z()/mm;
   }
    mMuMult += 1;
    if(maxMuMult < mMuMult)
      maxMuMult = mMuMult;
  }

}

