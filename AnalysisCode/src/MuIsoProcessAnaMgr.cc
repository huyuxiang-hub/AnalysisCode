
#include "MuIsoProcessAnaMgr.hh"
#include <sstream>
#include <cassert>
#include "junoHit_PMT.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SystemOfUnits.hh"
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
#include "G4Neutron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"

#include "NormalTrackInfo.hh"
#include "G4OpticalPhoton.hh"
#include <vector>
DECLARE_TOOL(MuIsoProcessAnaMgr);

MuIsoProcessAnaMgr::MuIsoProcessAnaMgr(const std::string& name)
    : ToolBase(name)
{

}

MuIsoProcessAnaMgr::~MuIsoProcessAnaMgr()
{

}

// == Run Action
void
MuIsoProcessAnaMgr::BeginOfRunAction(const G4Run* /*aRun*/) {
    SniperPtr<RootWriter> svc(*getParent(),"RootWriter");
   m_muonIso_tree = svc->bookTree("SIMEVT/muIso", "isotope and neutron from muon spallation.");
   m_muonIso_tree->Branch("evtID", &m_eventID, "evtID/I");
   m_muonIso_tree->Branch("IsoNum", &m_IsoN, "IsoNum/I");
   m_muonIso_tree->Branch("IsoName", m_IsoName, "IsoName[IsoNum]/I");
   m_muonIso_tree->Branch("IsoProc", m_IsoProc, "IsoProc[IsoNum]/I");
   m_muonIso_tree->Branch("IsoPosx", m_IsoPosx, "IsoPosx[IsoNum]/D");
   m_muonIso_tree->Branch("IsoPosy", m_IsoPosy, "IsoPosy[IsoNum]/D");
   m_muonIso_tree->Branch("IsoPosz", m_IsoPosz, "IsoPosz[IsoNum]/D");
   m_muonIso_tree->Branch("IsoKinE", m_IsoKinE, "IsoKinE[IsoNum]/D");
   m_muonIso_tree->Branch("Isotime", m_IsoTime, "IsoTime[IsoNum]/D");

   m_muonIso_tree->Branch("maxSecondary", &m_maxSecondary, "maxSecondary/I");
   m_muonIso_tree->Branch("PartName",&m_PartName);
   m_muonIso_tree->Branch("PartProc",&m_PartProc);
   m_muonIso_tree->Branch("VolMaterial",&m_VolMaterial);
   m_muonIso_tree->Branch("PartPosx",&m_PartPosx);
   m_muonIso_tree->Branch("PartPosy",&m_PartPosy);
   m_muonIso_tree->Branch("PartPosz",&m_PartPosz);
   m_muonIso_tree->Branch("PartTime",&m_PartTime);
   m_muonIso_tree->Branch("PartKine",&m_PartKinE);
   m_muonIso_tree->Branch("ParentKine",&m_ParentKinE);
   m_muonIso_tree->Branch("ParentID",&m_ParentID);
   m_muonIso_tree->Branch("ParentPartName", &m_ParentPartName);
    m_muonIso_tree->Branch("ParentStepNo", &m_ParentStepNo);
   m_muonIso_tree->Branch("ParentELoss", &m_ELoss);
   
   m_muonIso_tree->Branch("LosN", &m_LosN, "LosN/I");
   m_muonIso_tree->Branch("LosPosx", m_LosPosx, "LosPosx[LosN]/D");
   m_muonIso_tree->Branch("LosPosy", m_LosPosy, "LosPosy[LosN]/D");
   m_muonIso_tree->Branch("LosPosz", m_LosPosz, "LosPosz[LosN]/D");
   m_muonIso_tree->Branch("LosE", m_LosE, "LosE[LosN]/D");
   m_muonIso_tree->Branch("LosStpL", m_LosStpL, "LosStpL[LosN]/D");
   m_muonIso_tree->Branch("Lostime", m_Lostime, "Lostime[LosN]/D");
}

void
MuIsoProcessAnaMgr::EndOfRunAction(const G4Run* /*aRun*/) {

}

// == Event Action
void
MuIsoProcessAnaMgr::BeginOfEventAction(const G4Event* evt) {
    m_eventID = evt->GetEventID();
   ///////////////isotope //////////////////////////////
  
     m_IsoN =0;
  //   m_isoName->clear();
  //   m_isoProc->clear();
    for(int j=0;j<maxIsoN;j++) {
      m_IsoName[j] =0;
      m_IsoProc[j] =0;
      m_IsoPosx[j]  = 0;
      m_IsoPosy[j] =  0;
      m_IsoPosz[j] =  0;
      m_IsoKinE[j] = 0;
      m_IsoTime[j] = 0;
     }
   // for vector storage
    newTrackId =0;
    secSize =0;
    startNum =0;
    m_maxSecondary =0;  
    m_PartName.clear();
    m_PartProc.clear();
    m_VolMaterial.clear();
    m_PartPosx.clear();
    m_PartPosy.clear();
    m_PartPosz.clear();
    m_PartTime.clear();
    m_PartKinE.clear(); 
    m_ParentKinE.clear();
    m_ParentID.clear();
    m_ParentPartName.clear(); 
    m_ParentStepNo.clear();
    m_ELoss.clear();  
 
 //large energy loss
   m_LosN=0;
   for(int j=0;j<maxLosN;j++) {
      m_LosPosx[j]  = 0;
      m_LosPosy[j] =  0;
      m_LosPosz[j] =  0;
      m_LosE[j] = 0;
      m_LosStpL[j] = 0;
      m_Lostime[j] = 0;
     }

}

void
MuIsoProcessAnaMgr::EndOfEventAction(const G4Event* /*evt*/) {
    m_muonIso_tree->Fill();
}

// == Tracking Action
void
MuIsoProcessAnaMgr::PreUserTrackingAction(const G4Track* aTrack) {
    if(aTrack->GetParentID()==0 && aTrack->GetUserInformation()==0)
    {
        NormalTrackInfo* anInfo = new NormalTrackInfo(aTrack);
        G4Track* theTrack = (G4Track*)aTrack;
        theTrack->SetUserInformation(anInfo);
    }
    //NormalTrackInfo* info = (NormalTrackInfo*)(aTrack->GetUserInformation());

//save the parent track name
    //if(info && aTrack->GetDefinition()!=G4OpticalPhoton::Definition()) {
    //     const G4String& parent_name = info->getParentName(); 
//         G4cout << "PreUserAction The parent of " << aTrack->GetDefinition()->GetParticleName()<<",trackId="<<aTrack->GetTrackID() 
  //                << " is " << parent_name<<"; parentId="<<aTrack->GetParentID() 
    //             <<";kine= "<<aTrack->GetKineticEnergy()/MeV
      //           << G4endl; 
         
      //}
}
void 
MuIsoProcessAnaMgr::PostUserTrackingAction(const G4Track* aTrack) {
////////////////////////////////////////////////////////
  G4TrackingManager* tm = G4EventManager::GetEventManager()
                                            -> GetTrackingManager();
    G4TrackVector* secondaries = tm->GimmeSecondaries();
    if(secondaries)
    {
        NormalTrackInfo* info = (NormalTrackInfo*)(aTrack->GetUserInformation());
        size_t nSeco = secondaries->size();
        if(nSeco>0)
        {
            for(size_t i=0;i<nSeco;i++)
            {
               //  make sure the secondaries' track info is empty
                // if already attached track info, skip it.
                if ((*secondaries)[i]->GetUserInformation()) {
                    LogDebug << "The secondary already has user track info. skip creating new one" << std::endl;
                    continue;
                }
                
                NormalTrackInfo* infoNew = new NormalTrackInfo(info);
                // save parent track info 
                infoNew->setParentName(aTrack->GetDefinition()->GetParticleName());
                (*secondaries)[i]->SetUserInformation(infoNew);
            }
        }
    }
/*
//////////////////////////////////
  G4int StepNo = aTrack->GetCurrentStepNumber();
  G4int trackId = aTrack ->GetTrackID();
  G4ParticleDefinition* particle = aTrack->GetDefinition();
  G4String pname = particle->GetParticleName();
G4TrackVector* fSecondary = tm->GimmeSecondaries();
  size_t nSeco = fSecondary->size();
if(nSeco>0) { 
G4cout<<"---Mother Id ="<<trackId<<"; name="<<pname<< "; stepNo ="<<StepNo<<"; secondary size="<<(*fSecondary).size()<<" KinE="<<aTrack->GetKineticEnergy()/MeV<<G4endl;
}
// all particle information save
 if((*fSecondary).size()>0) {
   for (size_t lp=0; lp<(*fSecondary).size(); lp++) {
     if(((*fSecondary)[lp]->GetDefinition()->GetParticleName()=="e+"||(*fSecondary)[lp]->GetDefinition()->GetParticleName()=="e-"||(*fSecondary)[lp]->GetDefinition()->GetParticleName()=="gamma")&&(*fSecondary)[lp]->GetKineticEnergy()/MeV<20)
      {
        continue;
      }
      m_PartName.push_back((*fSecondary)[lp]->GetDefinition()->GetParticleName());
    G4cout<<"secondary size ="<<(*fSecondary).size()<<"; j ="<<lp<<"; secondary pname ="<<(*fSecondary)[lp]->GetDefinition()->GetParticleName()<<"; kine="<<(*fSecondary)[lp]->GetKineticEnergy()/MeV<<" ;";
     if((*fSecondary)[lp]->GetCreatorProcess()) {
      m_PartProc.push_back((*fSecondary)[lp]->GetCreatorProcess()->GetProcessName());
      G4cout<<"secondary proc ="<<(*fSecondary)[lp]->GetCreatorProcess()->GetProcessName()<<G4endl;
      }
     else { m_PartProc.push_back(-1);}

      m_maxSecondary++;

   }//lp
  }// Size>0
/////////////////////////////////////////////////
*/
}

// == Stepping Action
void
MuIsoProcessAnaMgr::UserSteppingAction(const G4Step* fStep) {

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
    //G4double edep = fStep->GetTotalEnergyDeposit();

    bool needToSim = false;
    if (fTrack->GetDefinition() != G4OpticalPhoton::Definition()) {
        needToSim = true;
    }
    if (not needToSim) {
        return;
    }
  G4int StepNo = fTrack->GetCurrentStepNumber();
  G4int trackId = fTrack ->GetTrackID();
  G4StepPoint* thePrePoint;
  //G4StepPoint* thePostPoint;
  G4Material* material =  fTrack->GetMaterial();
  //G4Material* material_post;
  G4ParticleDefinition* particle = fTrack->GetDefinition();
  G4String pname = particle->GetParticleName();

  thePrePoint = fStep->GetPreStepPoint();
  //thePostPoint = fStep->GetPostStepPoint();
  material = thePrePoint->GetMaterial();
  //material_post = thePostPoint->GetMaterial();//WARNING, for the end step of particle, there is no postpoint only prepoint.
  //G4double TotalEnergyDeposit = fStep->GetTotalEnergyDeposit();
  //const G4DynamicParticle* aParticle = fTrack->GetDynamicParticle();
  //G4int parentId = fTrack->GetParentID();
  const CLHEP::Hep3Vector point = thePrePoint->GetPosition();
  G4int parId = fTrack->GetParentID();
  G4TrackingManager* tm = G4EventManager::GetEventManager()-> GetTrackingManager();
  G4TrackVector* fSecondary = tm->GimmeSecondaries();
  //size_t nSeco = fSecondary->size();
 int isoflag =0;
 int neuflag =0;
 int stpIso =0;
 int stpNeu =0;
   // if(pname == "mu-"||pname == "mu+")   G4cout<<"matername ="<<material->GetName()<<"; pos x ="<<point.x()/mm<<"; pos y="<<point.y()/mm<<"; pos z ="<<point.z()/mm<<"; px ="<<fTrack->GetMomentumDirection().x()<<"; py ="<<fTrack->GetMomentumDirection().y()<<"; pz ="<<fTrack->GetMomentumDirection().z()<<" ; theta ="<<180/3.14*atan(sqrt(fTrack->GetMomentumDirection().x()*fTrack->GetMomentumDirection().x()+fTrack->GetMomentumDirection().y()*fTrack->GetMomentumDirection().y())/fTrack->GetMomentumDirection().z())<<G4endl;
//  if(material->GetName()!="Rock" )
  if(material==s_mat_LS ) 
   { 
       if((particle==G4MuonPlus::Definition()||particle==G4MuonMinus::Definition())&&(-1*fStep->GetDeltaEnergy()/fStep->GetStepLength()*10>6&&-1*fStep->GetDeltaEnergy()>20)) 
        {
           m_LosPosx[m_LosN] = point.x()/mm;
           m_LosPosy[m_LosN] = point.y()/mm;
           m_LosPosz[m_LosN] = point.z()/mm; 
           m_LosE[m_LosN] =   (-1)*fStep->GetDeltaEnergy();
           m_LosStpL[m_LosN] = fStep->GetStepLength()/mm;
           m_Lostime[m_LosN]   =fTrack->GetGlobalTime()/(1000*ns) ; 
           m_LosN++;
           if(m_LosN>4998) {m_LosN=4999; G4cout<<"Warning: the energy loss reach the max. Los Bin number!"<<G4endl;}
        // cout<<"matername ="<<material->GetName()<<"; pos x ="<<point.x()/mm<<"; pos y="<<point.y()/mm<<"; pos z ="<<point.z()/mm<<"; quenchEd ="<<fStep->GetDeltaEnergy()<<"; stepLength ="<<fStep->GetStepLength()<<"; avg loss per MeV ="<<fStep->GetDeltaEnergy()/fStep->GetStepLength()*10.<<" ; n Secondary ="<<(*fSecondary).size()<<endl;
   //for (size_t lp=0; lp<(*fSecondary).size(); lp++) {
   //  if((*fSecondary)[lp]->GetDefinition()->GetParticleName()=="e-"||(*fSecondary)[lp]->GetDefinition()->GetParticleName()=="e+"||(*fSecondary)[lp]->GetDefinition()->GetParticleName()=="gamma"||(*fSecondary)[lp]->GetDefinition()->GetParticleName()=="mu+"||(*fSecondary)[lp]->GetDefinition()->GetParticleName()=="mu-") continue;
     // cout<<"SS  pname="<<(*fSecondary)[lp]->GetDefinition()->GetParticleName()<<" ; KinE energy ="<< (*fSecondary)[lp]->GetKineticEnergy()/MeV<<endl;
  // }//for 

     }
  for (size_t lp=0; lp<(*fSecondary).size(); lp++) {
    if ((*fSecondary)[lp]->GetDefinition()->GetParticleName()=="O15[0.0]"||
     (*fSecondary)[lp]->GetDefinition()->GetParticleName()=="B13[0.0]"||(*fSecondary)[lp]->GetDefinition()->GetParticleName()=="N13[0.0]"||
     (*fSecondary)[lp]->GetDefinition()->GetParticleName()=="B12[0.0]"||(*fSecondary)[lp]->GetDefinition()->GetParticleName()=="N12[0.0]"||(*fSecondary)[lp]->GetDefinition()->GetParticleName()=="C12[0.0]"||
     (*fSecondary)[lp]->GetDefinition()->GetParticleName()=="C11[0.0]"||(*fSecondary)[lp]->GetDefinition()->GetParticleName()=="Be11[0.0]"||(*fSecondary)[lp]->GetDefinition()->GetParticleName()=="Li11[0.0]"||
     (*fSecondary)[lp]->GetDefinition()->GetParticleName()=="C10[0.0]"||(*fSecondary)[lp]->GetDefinition()->GetParticleName()=="Be10[0.0]"||(*fSecondary)[lp]->GetDefinition()->GetParticleName()=="B10[0.0]"||
     (*fSecondary)[lp]->GetDefinition()->GetParticleName()=="C9[0.0]"||(*fSecondary)[lp]->GetDefinition()->GetParticleName()=="Li9[0.0]"||(*fSecondary)[lp]->GetDefinition()->GetParticleName()=="Be9[0.0]"||(*fSecondary)[lp]->GetDefinition()->GetParticleName()=="B9[0.0]"||
  (*fSecondary)[lp]->GetDefinition()->GetParticleName()=="He8[0.0]"||(*fSecondary)[lp]->GetDefinition()->GetParticleName()=="Li8[0.0]"||(*fSecondary)[lp]->GetDefinition()->GetParticleName()=="Be8[0.0]"||(*fSecondary)[lp]->GetDefinition()->GetParticleName()=="B8[0.0]"||

 (*fSecondary)[lp]->GetDefinition()->GetParticleName()=="Be7[0.0]"||(*fSecondary)[lp]->GetDefinition()->GetParticleName()=="Li7[0.0]"||
 (*fSecondary)[lp]->GetDefinition()->GetParticleName()=="He6[0.0]"||(*fSecondary)[lp]->GetDefinition()->GetParticleName()=="Li6[0.0]"||
 (*fSecondary)[lp]->GetDefinition()->GetParticleName()=="Li5[0.0]"||
 (*fSecondary)[lp]->GetDefinition()->GetParticleName()=="H4[0.0]"||
 (*fSecondary)[lp]->GetDefinition()->GetParticleName()=="alpha"||
 (*fSecondary)[lp]->GetDefinition()->GetParticleName()=="He3"||
 (*fSecondary)[lp]->GetDefinition()->GetParticleName()=="triton"||
 (*fSecondary)[lp]->GetDefinition()->GetParticleName()=="deutron"||
  (*fSecondary)[lp]->GetDefinition()->GetParticleName()=="proton")
  {
         isoflag =1;
         stpIso =StepNo; 
   }// if
      if ((*fSecondary)[lp]->GetDefinition()==G4Neutron::Definition()){
        neuflag =1;          
        stpNeu = StepNo;
     }
 }// for lp
  if((StepNo==1&&pname=="neutron"&&neuflag==0)||(StepNo==1&&(pname=="O15[0.0]"||pname=="B13[0.0]"||pname=="N13[0.0]"||pname=="B12[0.0]"||pname=="N12[0.0]"||pname=="C12[0.0]"||pname=="C11[0.0]"||pname=="Be11[0.0]"||pname=="Li11[0.0]"||pname=="C10[0.0]"||pname=="Be10[0.0]"||pname=="B10[0.0]"||pname=="C9[0.0]"||pname=="Li9[0.0]"||pname=="B9[0.0]"||pname=="Be9[0.0]"||pname=="He8[0.0]"||pname=="Li8[0.0]"||pname=="B8[0.0]"||pname=="Be8[0.0]"||pname=="Be7[0.0]"||pname=="Li7[0.0]"||pname=="He6[0.0]"||pname=="Li6[0.0]"||
pname=="Li5[0.0]"||pname=="H4[0.0]"||pname=="alpha"||pname=="He3"||pname=="triton"||pname=="deutron"||pname=="proton")))
{
       if(pname=="neutron")        {m_IsoName[m_IsoN] = 1;}
      else if(pname=="B12[0.0]")  {m_IsoName[m_IsoN] = 2;}
      else if(pname=="Li9[0.0]")  {m_IsoName[m_IsoN] = 3;}
      else if(pname=="He8[0.0]")  {m_IsoName[m_IsoN] = 4;}
      else if(pname=="N12[0.0]")  {m_IsoName[m_IsoN] = 5;}

      else if(pname=="O15[0.0]")  {m_IsoName[m_IsoN] = 6;}
      else if(pname=="B13[0.0]")   {m_IsoName[m_IsoN] = 7;}
      else if(pname=="N13[0.0]")   {m_IsoName[m_IsoN] = 8;}
      else if(pname=="C12[0.0]")  {m_IsoName[m_IsoN] = 9;}
      else if(pname=="C11[0.0]")  {m_IsoName[m_IsoN] = 10;}
      else if(pname=="Be11[0.0]") {m_IsoName[m_IsoN] = 11;}
      else if(pname=="Li11[0.0]")  {m_IsoName[m_IsoN] = 12;}
      else if(pname=="C10[0.0]")  {m_IsoName[m_IsoN] = 13;}
      else if(pname=="Be10[0.0]")  {m_IsoName[m_IsoN] = 14;}
      else if(pname=="B10[0.0]")  {m_IsoName[m_IsoN] = 15;}

      else if(pname=="C9[0.0]") {m_IsoName[m_IsoN] = 16;}
      else if(pname=="B9[0.0]") {m_IsoName[m_IsoN] = 17;}
      else if(pname=="Be9[0.0]") {m_IsoName[m_IsoN] = 18;}
      else if(pname=="Li8[0.0]") {m_IsoName[m_IsoN] = 19;}
      else if(pname=="B8[0.0]") {m_IsoName[m_IsoN] = 20;}
      else if(pname=="Be8[0.0]") {m_IsoName[m_IsoN] = 21;}
      else if(pname=="Be7[0.0]") {m_IsoName[m_IsoN] = 22;}
      else if(pname=="Li7[0.0]") {m_IsoName[m_IsoN] = 23;}
      else if(pname=="He6[0.0]") {m_IsoName[m_IsoN] = 24;}
      else if(pname=="Li6[0.0]") {m_IsoName[m_IsoN] = 25;}
      else if(pname=="Li5[0.0]") {m_IsoName[m_IsoN] = 26;}
      else if(pname=="H4[0.0]") {m_IsoName[m_IsoN] = 27;}
      else if(pname=="He3") {m_IsoName[m_IsoN] = 28;}
      else if(pname=="alpha") {m_IsoName[m_IsoN] = 29;}
      else if(pname=="triton") {m_IsoName[m_IsoN] = 30;}
      else if(pname=="deutron") {m_IsoName[m_IsoN] = 31;}
      else if(pname=="proton") {m_IsoName[m_IsoN] = 32;}

      else {G4cout<<"un-record particle name ="<<pname<<G4endl;  m_IsoName[m_IsoN] = 100;}
 
      if (fTrack->GetCreatorProcess()) { //Avoid bug happenes if generator is not muon. Jilei, 29-04-2020

 //     G4cout<<"Iso tree inofor. particle name ="<<pname<<"; gen proc ="<<fTrack->GetCreatorProcess()->GetProcessName()<<G4endl;
      if(fTrack->GetCreatorProcess()->GetProcessName()=="PhotonInelastic")          {m_IsoProc[m_IsoN] =1;}
      else if(fTrack->GetCreatorProcess()->GetProcessName()=="PionMinusInelastic")  {m_IsoProc[m_IsoN] =2;}
      else if(fTrack->GetCreatorProcess()->GetProcessName()=="NeutronInelastic")    {m_IsoProc[m_IsoN] =3;}
      else if(fTrack->GetCreatorProcess()->GetProcessName()=="PionInelastic")       {m_IsoProc[m_IsoN] =4;}
      else if(fTrack->GetCreatorProcess()->GetProcessName()=="ProtonInelastic")     {m_IsoProc[m_IsoN] =5;}
      else if(fTrack->GetCreatorProcess()->GetProcessName()=="MuMinusCaptureAtRest"){m_IsoProc[m_IsoN] =6;}
      else if(fTrack->GetCreatorProcess()->GetProcessName()=="ElectroNuclear")      {m_IsoProc[m_IsoN] =7;}
      else if(fTrack->GetCreatorProcess()->GetProcessName()=="KaonZeroLInelastic")  {m_IsoProc[m_IsoN] =8;}
      else if(fTrack->GetCreatorProcess()->GetProcessName()=="KaonPlusInelastic")   {m_IsoProc[m_IsoN] =9;}
      else if(fTrack->GetCreatorProcess()->GetProcessName()=="KaonMinusInelastic")  {m_IsoProc[m_IsoN] =10;}
      else if(fTrack->GetCreatorProcess()->GetProcessName()=="DeuteronInelastic")   {m_IsoProc[m_IsoN] =11;}
      else if(fTrack->GetCreatorProcess()->GetProcessName()=="RadioactiveDecay")    {m_IsoProc[m_IsoN] =12;}
      else if(fTrack->GetCreatorProcess()->GetProcessName()=="PionPlusInelastic")    {m_IsoProc[m_IsoN] =13;}
      else if(fTrack->GetCreatorProcess()->GetProcessName()=="muNuclear")            {m_IsoProc[m_IsoN] =14;}
      else if(fTrack->GetCreatorProcess()->GetProcessName()=="PositionNuclear")      {m_IsoProc[m_IsoN] =15;}
      else if(fTrack->GetCreatorProcess()->GetProcessName()=="AlphaInelastic")       {m_IsoProc[m_IsoN] =16;}
      else if(fTrack->GetCreatorProcess()->GetProcessName()=="muMinusCaptureAtRest") {m_IsoProc[m_IsoN] =17;}
      else if(fTrack->GetCreatorProcess()->GetProcessName()=="Decay")                {m_IsoProc[m_IsoN] =18;}
      else if(fTrack->GetCreatorProcess()->GetProcessName()=="neutronElastic")       {m_IsoProc[m_IsoN] =19;}
      else if(fTrack->GetCreatorProcess()->GetProcessName()=="LElastic")             {m_IsoProc[m_IsoN] =20;}
      else {G4cout<<"unrecord process name ="<<fTrack->GetCreatorProcess()->GetProcessName()<<G4endl; 
            m_IsoProc[m_IsoN] =100;
      }
      }
      else {
       G4cout<<"Unrecord process due to nullptr." << G4endl;
       m_IsoProc[m_IsoN] =100;
      }
      m_IsoPosx[m_IsoN]  = point.x()/mm;
      m_IsoPosy[m_IsoN] =  point.y()/mm;
      m_IsoPosz[m_IsoN] =  point.z()/mm;
      m_IsoKinE[m_IsoN]  =thePrePoint->GetKineticEnergy()/MeV;
      m_IsoTime[m_IsoN] = thePrePoint->GetGlobalTime()/(1000*ns);  // in unit of us
      m_IsoN++;
     if(m_IsoN>79998) {m_IsoN=79999;std::cout<<"Warning: isotope number reach the max.number!"<<std::endl; }
  }

 }// for material LS
//if((*fSecondary).size()>0)
 //{
 //G4cout<<"---Stepping motherId ="<<trackId<<"; name="<<pname<< "; stepNo ="<<StepNo<<"; secondary="<<(*fSecondary).size()<<"; Kine ="<<fTrack->GetKineticEnergy()/MeV<<"; pos x="<<point.x()/mm<<G4endl;
//}
//
//
if(material!=s_mat_Rock ){
  // for all particle information save
 if((*fSecondary).size()>0) {
 if(trackId!=newTrackId) {
   newTrackId = trackId;
   secSize = (*fSecondary).size(); 
   startNum = 0;
//   G4cout<<"Track!=newTrack;new TrackId="<<newTrackId<<"; secSize ="<<secSize<<"; startnum="<<startNum<<G4endl;
  }
 else if(trackId==newTrackId){
   startNum = secSize; 
   secSize = (*fSecondary).size();
 // G4cout<<"track=newTrack, new TrackId="<<newTrackId<<"; secSize ="<<secSize<<"; startnum="<<startNum<<G4endl;
  }
 }//size>0
 if((*fSecondary).size()>0) {
   for (size_t lp=startNum; lp<(*fSecondary).size(); lp++) {
     if(((*fSecondary)[lp]->GetCreatorProcess()) 
        and (
          ((*fSecondary)[lp]->GetCreatorProcess()->GetProcessName()=="muIoni")
          ||((*fSecondary)[lp]->GetCreatorProcess()->GetProcessName()=="muPairProd")
          ||((*fSecondary)[lp]->GetCreatorProcess()->GetProcessName()=="LowEnergyIoni")
          ||((*fSecondary)[lp]->GetCreatorProcess()->GetProcessName()=="LowEnBrem")
          ||((*fSecondary)[lp]->GetCreatorProcess()->GetProcessName()=="LowEnCompton")
          ||((*fSecondary)[lp]->GetCreatorProcess()->GetProcessName()=="eIoni")
          ||((*fSecondary)[lp]->GetCreatorProcess()->GetProcessName()=="eBrem")
          ||((*fSecondary)[lp]->GetCreatorProcess()->GetProcessName()=="annihil")
          ||((*fSecondary)[lp]->GetCreatorProcess()->GetProcessName()=="LowEnConversion")
          ||((*fSecondary)[lp]->GetCreatorProcess()->GetProcessName()=="muBrems")
          ||((*fSecondary)[lp]->GetCreatorProcess()->GetProcessName()=="LowEnPhotoElec")
          ||((*fSecondary)[lp]->GetCreatorProcess()->GetProcessName()=="neutronElastic")
          ||((*fSecondary)[lp]->GetCreatorProcess()->GetProcessName()=="hLowEIoni") 
          ||((*fSecondary)[lp]->GetCreatorProcess()->GetProcessName()=="phot") //skip phot process
          ||((*fSecondary)[lp]->GetCreatorProcess()->GetProcessName()=="compt") //skip compt process
          ||((*fSecondary)[lp]->GetDefinition() == G4OpticalPhoton::Definition()) //skip opticalphoton saving, 29-04-2020
        )
      )
    {
     continue; //These processes will be not saved.
     }
      m_PartName.push_back((*fSecondary)[lp]->GetDefinition()->GetParticleName());
  //  G4cout<<"secondary size ="<<(*fSecondary).size()<<"; j ="<<lp<<";secondary pname ="<<(*fSecondary)[lp]->GetDefinition()->GetParticleName()<<"; kine="<<(*fSecondary)[lp]->GetKineticEnergy()/MeV<<"; parent Pname="<<pname<<"; StepNo="<<StepNo<<" ; parent kinE="<<fTrack->GetKineticEnergy()/MeV<<"; material="<<material->GetName()<<"; E loss ="<<-1*fStep->GetDeltaEnergy()/MeV;
     if((*fSecondary)[lp]->GetCreatorProcess()) {
      m_PartProc.push_back((*fSecondary)[lp]->GetCreatorProcess()->GetProcessName());
   //   G4cout<<"secondary proc ="<<(*fSecondary)[lp]->GetCreatorProcess()->GetProcessName()<<G4endl;
      } 
     else { m_PartProc.push_back(-1);}
       m_PartKinE.push_back((*fSecondary)[lp]->GetKineticEnergy()/MeV);
       m_PartTime.push_back(thePrePoint->GetGlobalTime()/(1000*ns));// in unit of us
       m_ParentID.push_back(parId);
       m_ParentPartName.push_back(pname);
       m_PartPosx.push_back(point.x()/mm);
       m_PartPosy.push_back(point.y()/mm);
       m_PartPosz.push_back(point.z()/mm);
       m_VolMaterial.push_back(material->GetName());
       m_ParentKinE.push_back(fTrack->GetKineticEnergy()/MeV);
       m_ParentStepNo.push_back(StepNo);
       m_ELoss.push_back(-1*fStep->GetDeltaEnergy()/MeV);
       m_maxSecondary++;
  //G4cout<<"pname ="<<pname<<"; proc ="<<fproc<<"; kine ="<<thePrePoint->GetKineticEnergy()/MeV<<G4endl;
   // G4cout<<"m_maxSecondary ="<<m_maxSecondary<<G4endl;
   }//lp
  }// Size>0
 }// Rock
}
