
#include "MuFastnProcessAnaMgr.hh"

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
#include "G4OpticalPhoton.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"

#include "G4SystemOfUnits.hh"
#include "NormalTrackInfo.hh"

DECLARE_TOOL(MuFastnProcessAnaMgr);

MuFastnProcessAnaMgr::MuFastnProcessAnaMgr(const std::string& name)
    : ToolBase(name)
{

}

MuFastnProcessAnaMgr::~MuFastnProcessAnaMgr()
{

}

// == Run Action
void
MuFastnProcessAnaMgr::BeginOfRunAction(const G4Run* /*aRun*/) {
    SniperPtr<RootWriter> svc(*getParent(), "RootWriter");
    m_fastn_tree = svc->bookTree("SIMEVT/mufastn", "muon induced fast neutron");
    m_fastn_tree -> Branch("evtID", &m_eventID, "evtID/I");
    m_fastn_tree -> Branch("stepNumber", &mStepNumber, "stepNumber/I");
    m_fastn_tree -> Branch("fnPDG", mfnPDG, "fnPDG[stepNumber]/I");
    m_fastn_tree -> Branch("fnPositionx", mfnPositionx, "fnPositionx[stepNumber]/D");
    m_fastn_tree -> Branch("fnPositiony", mfnPositiony, "fnPositiony[stepNumber]/D");
    m_fastn_tree -> Branch("fnPositionz", mfnPositionz, "fnPositionz[stepNumber]/D");
    m_fastn_tree -> Branch("fnGlobalTime", mfnGlobalTime, "fnGlobalTime[stepNumber]/D");
    m_fastn_tree -> Branch("fnLocalTime", mfnLocalTime, "fnLocalTime[stepNumber]/D");
    m_fastn_tree -> Branch("fnTotalEnergy", mfnTotalEnergy, "fnTotalEnergy[stepNumber]/D");
    m_fastn_tree -> Branch("fnKineticEnergy", mfnKineticEnergy, "fnKineticEnergy[stepNumber]/D");
    m_fastn_tree -> Branch("fnEnergyDeposit", mfnEnergyDeposit, "fnEnergyDeposit[stepNumber]/D");
    m_fastn_tree -> Branch("fnQEnergyDeposit", mfnQEnergyDeposit, "fnQEnergyDeposit[stepNumber]/D");
    m_fastn_tree -> Branch("fnCaptime", mfnCaptime, "fnCaptime[stepNumber]/D");

}

void
MuFastnProcessAnaMgr::EndOfRunAction(const G4Run* /*aRun*/) {

}

// == Event Action
void
MuFastnProcessAnaMgr::BeginOfEventAction(const G4Event* evt) {
  m_eventID = evt->GetEventID();
  mMuonFlag = false;
  mStepNumber = 0;
  for(int i = 0; i < maxStepNumber; i++){
    mfnPDG[i] = 0;
    mfnParId[i] = 0;
    mfnPositionx[i] = 0.;
    mfnPositiony[i] = 0.;
    mfnPositionz[i] = 0.;
    mfnGlobalTime[i] = 0.;
    mfnLocalTime[i] = 0.;
    mfnTotalEnergy[i] = 0.;
    mfnKineticEnergy[i] = 0.;
    mfnEnergyDeposit[i] = 0.;
    mfnQEnergyDeposit[i] = 0.;
    mfnCaptime[i] = 0.;
  }

}

void
MuFastnProcessAnaMgr::EndOfEventAction(const G4Event* /*evt*/) {
  //std::cout<<" ---------------------- Before sorting -------------------"<<std::endl;
  //  for(int i=0; i < mStepNumber; i++)
  //    std::cout<<" "<<i<<", "<<mfnGlobalTime[i]<<", "<<mfnLocalTime[i]<<", "<<mfnCaptime[i]<<", E "<<mfnEnergyDeposit[i]<<", "<<mfnQEnergyDeposit[i]<<", "<<mfnKineticEnergy[i]<<", "<<mfnTotalEnergy[i]<<", "<<mfnPDG[i]<<", "<<mfnParId[i]<<", "<<mfnPositionx[i]<<", "<<mfnPositiony[i]<<", "<<mfnPositionz[i]<<std::endl;
  //std::cout<<" ---------------------- Before sorting END -------------------"<<std::endl;
  //////  heap sort/////////
  /// sort globaltime,other parameters follow globaltime:
  heap(mfnGlobalTime, mfnLocalTime, mfnCaptime,
       mfnEnergyDeposit, mfnQEnergyDeposit, mfnKineticEnergy, mfnTotalEnergy,
       mfnPDG, mfnParId, 
       mfnPositionx, mfnPositiony, mfnPositionz,
       mStepNumber);
  //std::cout<<" ---------------------- After sorting -------------------"<<std::endl;
  //  for(int i=0; i < mStepNumber; i++)
  //    std::cout<<" "<<i<<", "<<mfnGlobalTime[i]<<", "<<mfnLocalTime[i]<<", "<<mfnCaptime[i]<<", E "<<mfnEnergyDeposit[i]<<", "<<mfnQEnergyDeposit[i]<<", "<<mfnKineticEnergy[i]<<", "<<mfnTotalEnergy[i]<<", "<<mfnPDG[i]<<", "<<mfnParId[i]<<", "<<mfnPositionx[i]<<", "<<mfnPositiony[i]<<", "<<mfnPositionz[i]<<std::endl;
  //std::cout<<" ---------------------- After sorting END -------------------"<<std::endl;

  //---------- begin merge steps each time window (50ns) ------------
  const int maxsub = 50000;
  Double_t tbeg[maxsub] = {0};
  //  Double_t tend[maxsub] = {0};
  Double_t tlocal[maxsub];
  Double_t esub[maxsub];
  Double_t eqde[maxsub];
  Double_t eke[maxsub];
  Double_t etotal[maxsub];
  Double_t posx[maxsub], posy[maxsub], posz[maxsub];

  Int_t khit[maxsub] = {0}; 
  Int_t ksub;
  ksub = -1;
  Double_t twin = 25;
  //initialize
  if(mStepNumber > 0)
  {
    ksub = 0;
    tbeg[ksub] = mfnGlobalTime[0];
    //tend[ksub] = mfnGlobalTime[0];
    tlocal[ksub] = mfnLocalTime[0];
    esub[ksub]   = mfnEnergyDeposit[0];
    eqde[ksub]   = mfnQEnergyDeposit[0];
    eke[ksub]    = mfnKineticEnergy[0];
    etotal[ksub] = mfnTotalEnergy[0];
    khit[ksub]   = 1;
    posx[ksub]   = mfnPositionx[0]*mfnEnergyDeposit[0];
    posy[ksub]   = mfnPositiony[0]*mfnEnergyDeposit[0];
    posz[ksub]   = mfnPositionz[0]*mfnEnergyDeposit[0];
    if(mfnCaptime[ksub] > 0.01) 
    {mfnPDG[ksub] = 2112;}
    else
    {mfnPDG[ksub] = 0; mfnCaptime[ksub] = 0;}
  }
  //loop merge
  for(int istep = 1; istep < mStepNumber; istep++)
  {
    if(twin < 25) twin = 25;
    if(istep > 0)
    {
      if(mfnGlobalTime[istep] < mfnGlobalTime[istep-1])
        std::cout<<"sort error!!!!!!!!!!!!"<<std::endl;
    }
    if(mfnCaptime[istep] > 0.01)
    {
      if((mfnCaptime[istep]-tbeg[ksub]) < 0) std::cout<<"captime error!"<<std::endl;
      else
      { if((mfnCaptime[istep]-tbeg[ksub]) < twin)
        twin = mfnCaptime[istep]-tbeg[ksub];
      }
    }
    if((mfnGlobalTime[istep]-tbeg[ksub]) < twin )
    {
      esub[ksub] = esub[ksub] +  mfnEnergyDeposit[istep];
      eqde[ksub] = eqde[ksub] +  mfnQEnergyDeposit[istep];
      eke[ksub] = eke[ksub] +  mfnKineticEnergy[istep];
      etotal[ksub] = etotal[ksub] +  mfnTotalEnergy[istep];
      posx[ksub] = posx[ksub]+mfnPositionx[istep]*mfnEnergyDeposit[istep];
      posy[ksub] = posy[ksub]+mfnPositiony[istep]*mfnEnergyDeposit[istep];
      posz[ksub] = posz[ksub]+mfnPositionz[istep]*mfnEnergyDeposit[istep];
      khit[ksub] = khit[ksub] + 1;
    }
    else
      if((ksub < maxsub))
      {
        ksub = ksub + 1;
        khit[ksub] = 1;

        tbeg[ksub] = mfnGlobalTime[istep];
        //tend[ksub-1] = mfnGlobalTime[istep-1];
        tlocal[ksub] = mfnLocalTime[istep];
        esub[ksub] = mfnEnergyDeposit[istep];
        eqde[ksub] = mfnQEnergyDeposit[istep];
        eke[ksub]  = mfnKineticEnergy[istep];
        etotal[ksub] =  mfnTotalEnergy[istep];
        posx[ksub] = mfnPositionx[istep]*mfnEnergyDeposit[istep];
        posy[ksub] = mfnPositiony[istep]*mfnEnergyDeposit[istep];
        posz[ksub] = mfnPositionz[istep]*mfnEnergyDeposit[istep];

        if(khit[ksub-1] != 0){
          eke[ksub-1] = eke[ksub-1]/khit[ksub-1];
          etotal[ksub-1] = etotal[ksub-1]/khit[ksub-1];
        }else { eke[ksub-1] = 0; etotal[ksub-1] = 0;}
        if(esub[ksub-1] != 0){
          posx[ksub-1] = posx[ksub-1]/esub[ksub-1];
          posy[ksub-1] = posy[ksub-1]/esub[ksub-1];
          posz[ksub-1] = posz[ksub-1]/esub[ksub-1];
        }
        else{posx[ksub-1] = 0; posy[ksub-1] = 0; posz[ksub-1] = 0;}
        if(mfnCaptime[istep] > 0.01)
        {mfnPDG[ksub] = 2112; mfnCaptime[ksub] = mfnCaptime[istep];}
        else
        {mfnPDG[ksub] = 0; mfnCaptime[ksub] = 0;}
        //          std::cout<<"ksub="<<ksub<<std::endl;
      }
      else std::cout<<" error: maxsub step is too small, is "<<maxsub<<std::endl;
    if(istep == (mStepNumber-1)) 
    {   //tend[ksub] = mfnGlobalTime[istep];
      eke[ksub]    = eke[ksub]/khit[ksub];
      etotal[ksub] = etotal[ksub]/khit[ksub];
      if(esub[ksub]!=0){
        posx[ksub] = posx[ksub]/esub[ksub];
        posy[ksub] = posy[ksub]/esub[ksub];
        posz[ksub] = posz[ksub]/esub[ksub];
      }
      else {posx[ksub] = 0; posy[ksub] = 0; posz[ksub] = 0;}

    }
  }//for
  mStepNumber = ksub + 1;

  for(int isub=0;isub<mStepNumber;isub++)
  {
    mfnGlobalTime[isub]    = tbeg[isub];
    mfnLocalTime[isub]     = tlocal[isub];
    mfnEnergyDeposit[isub] = esub[isub];
    mfnQEnergyDeposit[isub]= eqde[isub];
    mfnKineticEnergy[isub] = eke[isub];
    mfnTotalEnergy[isub]   = etotal[isub];
    mfnPositionx[isub]     = posx[isub]; 
    mfnPositiony[isub]     = posy[isub]; 
    mfnPositionz[isub]     = posz[isub]; 
  }
  ////------------- merge END -----------------
    //std::cout<<"------------- After  merging ------------------- "<<std::endl;
    //for(int i=0; i<mStepNumber; i++)
    //  std::cout<<" "<<i<<", "<<mfnGlobalTime[i]<<", "<<mfnLocalTime[i]<<", "<<mfnCaptime[i]<<", E "<<mfnEnergyDeposit[i]<<", "<<mfnQEnergyDeposit[i]<<", "<<mfnKineticEnergy[i]<<", "<<mfnTotalEnergy[i]<<", "<<mfnPDG[i]<<", "<<mfnParId[i]<<", "<<mfnPositionx[i]<<", "<<mfnPositiony[i]<<", "<<mfnPositionz[i]<<std::endl;
    //std::cout<<"------------- After merging END ------------------- "<<std::endl;
  if(mStepNumber>0){
    m_fastn_tree->Fill();}
}

// == Tracking Action
void
MuFastnProcessAnaMgr::PreUserTrackingAction(const G4Track* /*aTrack*/) {
    //NormalTrackInfo* info = (NormalTrackInfo*)(aTrack->GetUserInformation());
}

void
MuFastnProcessAnaMgr::PostUserTrackingAction(const G4Track* /*aTrack*/) {
}

// == Stepping Action
void
MuFastnProcessAnaMgr::UserSteppingAction(const G4Step* fStep) {
  G4Track* fTrack = fStep->GetTrack();
  //G4int StepNo = fTrack->GetCurrentStepNumber();
    G4double edep = fStep->GetTotalEnergyDeposit();

    // cache the material pointer, avoid comparison using string.
    static G4Material* s_mat_LS = nullptr;
    if (!s_mat_LS) {
        s_mat_LS = G4Material::GetMaterial("LS");
    }

    bool needToSim = false;
    if (edep > 0 and fTrack->GetDefinition() != G4OpticalPhoton::Definition()
                 and fTrack->GetMaterial() == s_mat_LS) {
        needToSim = true;
    }
    if (not needToSim) {
        return;
    }

  G4int trackId = fTrack ->GetTrackID();
  //G4StepPoint* thePrePoint;
  //G4StepPoint* thePostPoint;
  //G4Material* material =  fTrack->GetMaterial();
  //G4Material* material_post;
  G4ParticleDefinition* particle = fTrack->GetDefinition();
  // G4String pname = particle->GetParticleName();

  //thePrePoint = fStep->GetPreStepPoint();
  //thePostPoint = fStep->GetPostStepPoint();
  //material = thePrePoint->GetMaterial();
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
  // for neutron capture info 
  if(mStepNumber>=maxStepNumber)
  {
    mStepNumber=maxStepNumber-1;
  }
  if(fTrack->GetMaterial() == s_mat_LS){//=prePoint
    // Birks' law
    G4double dE_dx = edep/fStep->GetStepLength();
    G4double delta = dE_dx/fStep->GetTrack()->GetMaterial()->GetDensity();
    G4double birk1 = 0.0125*g/cm2/MeV;
    if(particle->GetPDGCharge()>1.1)//for particle charge greater than 1.
      birk1 = 0.57*birk1;
    G4double birk2 = 0;
    //birk2 = (0.0031*g/MeV/cm2)*(0.0031*g/MeV/cm2);
    G4double QuenchedTotalEnergyDeposit
      = edep/(1+birk1*delta+birk2*delta*delta);
    if(QuenchedTotalEnergyDeposit>0.001 ||
        (particle->GetPDGEncoding()==2112&&
         (fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="nCapture")))
    {
      mfnPDG[mStepNumber]    = particle->GetPDGEncoding();
      mfnParId[mStepNumber] = parId;
      mfnPositionx[mStepNumber]  = fTrack->GetPosition().x()/mm;//=post point
      mfnPositiony[mStepNumber]  = fTrack->GetPosition().y()/mm;
      mfnPositionz[mStepNumber]  = fTrack->GetPosition().z()/mm;
      mfnGlobalTime[mStepNumber] = fTrack->GetGlobalTime()/ns;
      mfnLocalTime[mStepNumber]  = fTrack->GetLocalTime()/ns;
      mfnTotalEnergy[mStepNumber]   = fTrack->GetTotalEnergy()/MeV;
      mfnKineticEnergy[mStepNumber] = fTrack->GetKineticEnergy()/MeV;
      mfnEnergyDeposit[mStepNumber] = fStep->GetTotalEnergyDeposit()/MeV;
      mfnQEnergyDeposit[mStepNumber] = QuenchedTotalEnergyDeposit;
      if(particle->GetPDGEncoding()==2112 && 
          fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "nCapture")
        mfnCaptime[mStepNumber] = fTrack->GetGlobalTime()/ns;
      else
        mfnCaptime[mStepNumber] = 0;
      mStepNumber++;
      if(mStepNumber==maxStepNumber-1)
        G4cout<<"WARNING: maxStepNumber not enough!!!!!!!!!!"<<G4endl;
    }
  }

}

