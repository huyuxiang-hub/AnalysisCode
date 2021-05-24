#include <boost/python.hpp>
#include <iostream>

#include "DepositEnergyTTAnaMgr.hh"
#include "NormalTrackInfo.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/SniperDataPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "RootWriter/RootWriter.h"
#include "G4Event.hh"
#include "G4Box.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4LossTableManager.hh"
// using CLHEP's random engine
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4OpticalPhoton.hh"

#include "EvtNavigator/NavBuffer.h"
#include "BufferMemMgr/IDataMemMgr.h"
#include "Event/SimHeader.h"
#include "Event/SimTTHit.h"

#include "DetSimAlg/DetSimAlg.h"
#include "DetSimAlg/IDetElement.h"
#include "Geometry/TTGeomSvc.h"
#include "Identifier/TtID.h"

#include "TROOT.h"
using namespace CLHEP;
DECLARE_TOOL(DepositEnergyTTAnaMgr);

DepositEnergyTTAnaMgr::DepositEnergyTTAnaMgr(const std::string& name)
  : ToolBase(name)
{
  declProp("EnableNtuple", m_flag_ntuple=false);
  m_evt_treeTT = 0;
  m_evt_treeTTDigit = 0;
  
  m_rnd = 0;
  coeff[0][0]=0.80;
  coeff[1][0]=0.013;
  coeff[2][0]=8.9;
  coeff[3][0]=1.92E-3;
  
  coeff[0][1]=2.24;
  coeff[1][1]=9.52E-3;
  coeff[2][1]=8.38;
  coeff[3][1]=1.73E-3;
  
  coeff[0][2]=3.76;
  coeff[1][2]=9.515E-3;
  coeff[2][2]=7.61;
  coeff[3][2]=1.57E-3;
  
  coeff[0][3]=5.65;
  coeff[1][3]=5.12E-3;
  coeff[2][3]=5.08;
  coeff[3][3]=1.015E-3;
  
  coeff[0][4]=4.08;
  coeff[1][4]=8.48E-3;
  coeff[2][4]=7.39;
  coeff[3][4]=1.49E-3;
  
  coeff[0][5]=3.73;
  coeff[1][5]=7.04E-3;
  coeff[2][5]=7.56;
  coeff[3][5]=1.54E-3;
  
  coeff[0][6]=3.2;
  coeff[1][6]=7.37E-3;
  coeff[2][6]=8.09;
  coeff[3][6]=1.63E-3;
  
  coeff[0][7]=4.73;
  coeff[1][7]=6.30E-3;
  coeff[2][7]=6.63;
  coeff[3][7]=1.31E-3;
  
  coeff[0][8]=4.7;
  coeff[1][8]=6.23E-3;
  coeff[2][8]=6.76;
  coeff[3][8]=1.34E-3;
  
  coeff[0][9]=4.19;
  coeff[1][9]=7.27E-3;
  coeff[2][9]=7.67;
  coeff[3][9]=1.53E-3;
  
  coeff[0][10]=4.48;
  coeff[1][10]=8.34E-3;
  coeff[2][10]=7.8;
  coeff[3][10]=1.52E-3;
  
  coeff[0][11]=5.6;
  coeff[1][11]=7.39E-3;
  coeff[2][11]=6.85;
  coeff[3][11]=1.317E-3;
  
  coeff[0][12]=4.54;
  coeff[1][12]=7.4E-3;
  coeff[2][12]=7.7;
  coeff[3][12]=1.52E-3;
  
  coeff[0][13]=4.83;
  coeff[1][13]=0.012;
  coeff[2][13]=8.31;
  coeff[3][13]=1.64E-3;
  
  coeff[0][14]=4.44;
  coeff[1][14]=9.8E-3;
  coeff[2][14]=8.57;
  coeff[3][14]=1.67E-3;
  
  coeff[0][15]=4.9;
  coeff[1][15]=0.0107;
  coeff[2][15]=8.3;
  coeff[3][15]=1.63E-3;
    
}

DepositEnergyTTAnaMgr::~DepositEnergyTTAnaMgr()
{
  
}

void
DepositEnergyTTAnaMgr::BeginOfRunAction(const G4Run* /*aRun*/) {
  m_rnd=new TRandom3();
  
  //read CrossTalk
  
  TTGeomSvc* tt_ct_svc = 0;
  SniperPtr<TTGeomSvc> svcct(getParent(), "TTGeomSvc");
  if (svcct.invalid()) {
    LogError << "Can't get TTGeomSvc. We can't initialize Cross Talk." << std::endl;
    assert(0);
  } else {
    LogInfo << "Retrieve TTGeomSvc successfully." << std::endl;
    tt_ct_svc = svcct.data();
  }
  
  for(int ii=0;ii<64;ii++)
    for(int kk=0;kk<21;kk++)
      {
	CTchannels[ii][kk]=tt_ct_svc->getCT(ii,kk);
      }
  
  SniperPtr<DetSimAlg> detsimalg(getParent(), "DetSimAlg");
  if (detsimalg.invalid()) {
    // std::cout << "Can't Load DetSimAlg" << std::endl;
    assert(0);
  }
  std::string name = "TopTrackerConstruction";
  ToolBase* t = detsimalg->findTool(name);
  assert(t);
  
  de=dynamic_cast<IDetElement*>(t);
  
#if 0
  std::ofstream WallCorrespondenceTXT("wall_correspondence.txt");
  for(int iwall=0;iwall<3*7*3; iwall++){
    int wall   = de->geom_info("WallNum", iwall);
    int layer  = de->geom_info("Wall.Layer" ,wall);
    int column = de->geom_info("Wall.Column",wall);
    int row    = de->geom_info("Wall.Row"   ,wall);
    WallCorrespondenceTXT << iwall << "\t" << layer << "\t" << column << "\t" << row << std::endl;
  }
  WallCorrespondenceTXT.close();

  // Positioning following code in ComputeDigitTT...
  std::ofstream ChannelPositionTXT("ChannelPosition.txt");
  G4LogicalVolume * TT_vol = de->getLV();
  for(int iwall=0; iwall < TT_vol->GetNoDaughters(); ++iwall){
    G4VPhysicalVolume * Wall_vol = TT_vol->GetDaughter(iwall);
    G4AffineTransform Wall_transf(Wall_vol->GetRotation(), Wall_vol->GetTranslation());
    int n_wall;
    sscanf(Wall_vol->GetName().c_str(),"pWall_%d",&n_wall);
    for(int iplane=0; iplane < Wall_vol->GetLogicalVolume()->GetNoDaughters(); ++iplane){
      G4VPhysicalVolume * Plane_vol = Wall_vol->GetLogicalVolume()->GetDaughter(iplane);
      G4AffineTransform Plane_transf(Plane_vol->GetRotation(), Plane_vol->GetTranslation());
      int n_plane;
      sscanf(Plane_vol->GetName().c_str(),"pPlane_%d",&n_plane);
      for(int ipannel=0; ipannel < Plane_vol->GetLogicalVolume()->GetNoDaughters(); ++ipannel){
        G4VPhysicalVolume * Pannel_vol = Plane_vol->GetLogicalVolume()->GetDaughter(ipannel);
        G4AffineTransform Pannel_transf(Pannel_vol->GetRotation(), Pannel_vol->GetTranslation());
        int n_pannel;
        sscanf(Pannel_vol->GetName().c_str(),"pPanel_%d",&n_pannel);
        G4VPhysicalVolume * Adhesive_vol = Pannel_vol->GetLogicalVolume()->GetDaughter(0);
        for(int istrip=0; istrip < Adhesive_vol->GetLogicalVolume()->GetNoDaughters(); ++istrip){
          G4VPhysicalVolume * Strip_vol = Adhesive_vol->GetLogicalVolume()->GetDaughter(istrip);
          G4AffineTransform Strip_transf(Strip_vol->GetRotation(), Strip_vol->GetTranslation());
          int n_strip;
          sscanf(Strip_vol->GetName().c_str(),"pCoating_%d",&n_strip);

          // G4AffineTransform TT_transf = Wall_transf * Plane_transf * Pannel_transf * Strip_transf;
          G4AffineTransform TT_transf = Strip_transf * Pannel_transf * Plane_transf * Wall_transf;
          G4ThreeVector TT_global_pos = TT_transf.TransformPoint(G4ThreeVector(0,0,0));

          G4ThreeVector x_disp = TT_transf.TransformAxis(G4ThreeVector(1,0,0));

          if(std::abs(x_disp.y()) < 0.5){ // Should be 0 or +-1 if everything aligned
            TT_global_pos.setX(1e9);
          }
          else{
            TT_global_pos.setY(1e9);
          }
          for(int side=0; side < 2; side++){
            int n_pmt  = TtID::pmt_from_geant(side, n_plane, n_pannel);
            unsigned int chid = TtID::getIntID(n_wall, n_pmt, n_strip);
            ChannelPositionTXT << chid << " " << TT_global_pos.x() << " " << TT_global_pos.y() << " " << TT_global_pos.z() << std::endl;
          }
        }
      }
    }
  }
  ChannelPositionTXT.close();
#endif

  if (not m_flag_ntuple) {
    return;
  }
  // check the RootWriter is Valid.
  SniperPtr<RootWriter> svc(getParent(), "RootWriter");
  if (svc.invalid()) {
    LogError << "Can't Locate RootWriter. If you want to use it, please "
	     << "enable it in your job option file."
	     << std::endl;
    return;
  }
  
  gROOT->ProcessLine("#include <vector>");
  
  m_evt_treeTT = svc->bookTree("SIMEVT/TT", "Deposit Energy TT");
  m_evt_treeTT->Branch("evtID", &m_eventID, "evtID/I");
  m_evt_treeTT->Branch("NDeposits",&m_NDeposits,"NDeposits/I");
  m_evt_treeTT->Branch("dep_trackID",&m_dep_trackID);
  m_evt_treeTT->Branch("dep_pdg",&m_dep_pdg);
  m_evt_treeTT->Branch("dep_E",&m_dep_E);
  m_evt_treeTT->Branch("dep_E0",&m_dep_E0);
  m_evt_treeTT->Branch("dep_nbar",&m_dep_nbar);
  m_evt_treeTT->Branch("dep_npanel",&m_dep_npanel);
  m_evt_treeTT->Branch("dep_nplane",&m_dep_nplane);
  m_evt_treeTT->Branch("dep_nwall",&m_dep_nwall);
  m_evt_treeTT->Branch("dep_x",&m_dep_x);
  m_evt_treeTT->Branch("dep_y",&m_dep_y);
  m_evt_treeTT->Branch("dep_z",&m_dep_z);
  m_evt_treeTT->Branch("dep_t",&m_dep_t);
  m_evt_treeTT->Branch("dep_dL",&m_dep_dL);
  m_evt_treeTT->Branch("dep_dR",&m_dep_dR);
  m_evt_treeTT->Branch("dep_peL",&m_dep_peL);
  m_evt_treeTT->Branch("dep_peR",&m_dep_peR);
  m_evt_treeTT->Branch("dep_tL",&m_dep_tL);
  m_evt_treeTT->Branch("dep_tR",&m_dep_tR);
  m_evt_treeTT->Branch("isMuonDeposits",&m_isMuonDeposits);
  
  //new with channel information (no Left or right anymore)
  m_evt_treeTTDigit = svc->bookTree("SIMEVT/TTDigit", "PE TT");
  m_evt_treeTTDigit->Branch("evtID", &m_eventID, "evtID/I");
  m_evt_treeTTDigit->Branch("NTouchedChannel",&m_NTouchedChannel,"NTouchedChannel/I");
  m_evt_treeTTDigit->Branch("TB_channelC",m_TB_channelC,"TB_channelC[NTouchedChannel]/I");
  m_evt_treeTTDigit->Branch("TB_DMchannel",m_TB_DMchannel,"TB_DMchannel[NTouchedChannel]/I");
  m_evt_treeTTDigit->Branch("TB_pe",m_TB_pe,"TB_pe[NTouchedChannel]/F");
  m_evt_treeTTDigit->Branch("TB_time",m_TB_time,"TB_time[NTouchedChannel]/F");
  m_evt_treeTTDigit->Branch("TB_ADC",m_TB_ADC,"TB_ADC[NTouchedChannel]/I");
  m_evt_treeTTDigit->Branch("TB_xcC",m_TB_xcC,"TB_xcC[NTouchedChannel]/F");
  m_evt_treeTTDigit->Branch("TB_ycC",m_TB_ycC,"TB_ycC[NTouchedChannel]/F");
  m_evt_treeTTDigit->Branch("TB_zcC",m_TB_zcC,"TB_zcC[NTouchedChannel]/F");
  m_evt_treeTTDigit->Branch("TB_is_ctC",&m_TB_is_ctC,"TB_is_ctC[NTouchedChannel]/I");
  m_evt_treeTTDigit->Branch("TB_isMuonDepositsC",m_TB_isMuonDepositsC,"TB_isMuonDepositsC[NTouchedChannel]/I");
}

void
DepositEnergyTTAnaMgr::EndOfRunAction(const G4Run* /*aRun*/) {
  
}

void
DepositEnergyTTAnaMgr::BeginOfEventAction(const G4Event* evt) {
  m_eventID = evt->GetEventID();
  
  m_NDeposits=0;
  m_dep_trackID.clear();
  m_dep_pdg.clear();
  m_dep_E.clear();
  m_dep_E0.clear();
  m_dep_nbar.clear();
  m_dep_npanel.clear();
  m_dep_nplane.clear();
  m_dep_nwall.clear();
  m_dep_x.clear();
  m_dep_y.clear();
  m_dep_z.clear();
  m_dep_t.clear();
  m_dep_dL.clear();
  m_dep_dR.clear();
  m_dep_peL.clear();
  m_dep_peR.clear();
  m_dep_tL.clear();
  m_dep_tR.clear();
  m_isMuonDeposits.clear();

  for(int i=0;i<31744;i++){
    m_TB_isMuonDeposits[i]=0;
  }
  for(int i=0;i<63488;i++){
    m_TB_isMuonDepositsC[i]=0;
  }
  
  m_NTouchedBar = 0;
  m_NTouchedChannel = 0;
}

void
DepositEnergyTTAnaMgr::EndOfEventAction(const G4Event* /*evt*/) {
  
  std::cout<<"end of event action "<<std::endl;
  ComputeDigitTT();//fill digitization tree
  
  if (m_flag_ntuple and m_evt_treeTT) {
    m_evt_treeTT -> Fill();
  }
  
  save_into_data_model();
}

void
DepositEnergyTTAnaMgr::PreUserTrackingAction(const G4Track* /*aTrack*/) {
  
}

void
DepositEnergyTTAnaMgr::PostUserTrackingAction(const G4Track* /*aTrack*/) {
  
}

void
DepositEnergyTTAnaMgr::UserSteppingAction(const G4Step* step) {
  G4Track* track = step->GetTrack();
  
  //NormalTrackInfo* trackinfo = (NormalTrackInfo*)track->GetUserInformation();
  //if (not trackinfo) {
  //  return;
  //}
  
  
  
  G4double edep = step->GetTotalEnergyDeposit();
  G4ThreeVector pos = (step->GetPreStepPoint()->GetPosition() + step->GetPostStepPoint()->GetPosition())/2;
  if (edep > 0.001 and track->GetDefinition() != G4OpticalPhoton::Definition() and track->GetMaterial()->GetName() == "Scintillator")
    {
      G4VPhysicalVolume* mother1 = track->GetTouchable()->GetVolume(1);
      G4String mother1name="";
      if(mother1)
	mother1name = mother1->GetName();
      
      G4VPhysicalVolume* mother3 = track->GetTouchable()->GetVolume(3);
      G4String mother3name="";
      if(mother3)
	mother3name = mother3->GetName();
      
      G4VPhysicalVolume* mother4 = track->GetTouchable()->GetVolume(4);
      G4String mother4name="";
      if(mother4)
	mother4name = mother4->GetName();
      
      G4VPhysicalVolume* mother5 = track->GetTouchable()->GetVolume(5);
      G4String mother5name="";
      if(mother5)
	mother5name = mother5->GetName();
      
      int barnum;
      sscanf(mother1name.c_str(),"pCoating_%d",&barnum);
      
      int panelnum;
      sscanf(mother3name.c_str(),"pPanel_%d",&panelnum);
      
      int planenum;
      sscanf(mother4name.c_str(),"pPlane_%d",&planenum);
      
      int wallnum;
      sscanf(mother5name.c_str(),"pWall_%d",&wallnum);
      
      m_dep_trackID.push_back(track->GetTrackID());
      m_dep_pdg.push_back(track->GetDefinition()->GetPDGEncoding());
      m_dep_E.push_back(edep);
      m_dep_E0.push_back(track->GetVertexKineticEnergy ());
      m_dep_nbar.push_back(barnum);
      m_dep_npanel.push_back(panelnum);
      m_dep_nplane.push_back(planenum);
      m_dep_nwall.push_back(wallnum);
      m_dep_x.push_back(pos.x());
      m_dep_y.push_back(pos.y());
      m_dep_z.push_back(pos.z());
      m_dep_t.push_back((step->GetPreStepPoint()->GetGlobalTime()+step->GetPostStepPoint()->GetGlobalTime())/2);
      int isMuon(0);
      if(m_dep_pdg[m_NDeposits]==13 || m_dep_pdg[m_NDeposits]==-13)
	isMuon = 1;
      m_isMuonDeposits.push_back(isMuon);

      //std::cout<<"wall number "<<wallnum<<std::endl;
      
      // add the extra length out from the bar to the PMT, they are different depending on the bar num
      // the fitting result is len=-2.45263806631*barnum+ 91.6040322685 (cm) only one side
      int barnum32 = barnum < 32 ? barnum:63-barnum;
      double a_fiber_len = -2.45263806631;
      double b_fiber_len = 91.6040322685;
      double fiber_len =  a_fiber_len*barnum32+ b_fiber_len; // cm

      // Transform coordinate system from JUNO global, to the one of the scintillator strip.
      // With this conversion, any rotation they might have suffered is taken into account
      // and we don't need to treat separately the x-aligned and y-aligned strips.
      G4ThreeVector local_pos = step->GetPreStepPoint()->GetTouchable()->GetHistory()->GetTopTransform().TransformPoint(pos);
      double length_scintillator = ((G4Box*) track->GetTouchable()->GetVolume(0)->GetLogicalVolume()->GetSolid())->GetXHalfLength();

      double distanceEL = ( length_scintillator + local_pos.x() )/cm + fiber_len;
      double distanceER = ( length_scintillator - local_pos.x() )/cm + fiber_len;
      
      double a0(0);
      double a1(0);
      double b0(0);
      double b1(0);
      
      while(a0<=0 || a1<=0 || b0<=0 || b1<=0)
        {
	  if(barnum<32)
            {
	      a0=G4RandGauss::shoot(coeff[0][int(barnum/2)],0.2*coeff[0][int(barnum/2)]);
	      b0=G4RandGauss::shoot(coeff[1][int(barnum/2)],0.1*coeff[1][int(barnum/2)]);
	      a1=G4RandGauss::shoot(coeff[2][int(barnum/2)],0.2*coeff[2][int(barnum/2)]);
	      b1=G4RandGauss::shoot(coeff[3][int(barnum/2)],0.1*coeff[3][int(barnum/2)]);	      
            }
	  else
            {
	      a0=G4RandGauss::shoot(coeff[0][int((63 - barnum)/2)],0.2*coeff[0][int((63 - barnum)/2)]);
	      b0=G4RandGauss::shoot(coeff[1][int((63 - barnum)/2)],0.1*coeff[1][int((63 - barnum)/2)]);
	      a1=G4RandGauss::shoot(coeff[2][int((63 - barnum)/2)],0.2*coeff[2][int((63 - barnum)/2)]);
	      b1=G4RandGauss::shoot(coeff[3][int((63 - barnum)/2)],0.1*coeff[3][int((63 - barnum)/2)]);
	      
            }
        }
      
      m_dep_dL.push_back(distanceEL);
      m_dep_dR.push_back(distanceER);
      m_dep_peL.push_back((edep/2.15) * (a0 * exp (-b0*distanceEL) + a1 * exp (-b1*distanceEL)));
      m_dep_peR.push_back((edep/2.15) * (a0 * exp (-b0*distanceER) + a1 * exp (-b1*distanceER)));

      //      std::cout<<"Nb p.e Left: "<<(edep/2.15) * (a0 * exp (-b0*distanceEL) + a1 * exp (-b1*distanceEL))<<" pe right: "<<(edep/2.15) * (a0 * exp (-b0*distanceER) + a1 * exp (-b1*distanceER))<<std::endl;
      
      // n of WLS fiber: 1.4 - 1.6 [http://www.phenix.bnl.gov/WWW/publish/donlynch/RXNP/Safety%20Review%206_22_06/Kuraray-PSF-Y11.pdf]
      // n of polystyrene: 1.55-1.59 [google]
      // => use n=1.6 to be conservative
      const double c_scintilator=c_light/1.6;
      const double time_spread_slope = -0.2265;
      // The random component below comes from measurements done where muons were traversing a
      // horizontal scintillator strip vertically, at the center of the strip. The measured
      // time difference between the time on each end followed an exponential with slope -0.2265
      // (time_spread_slope). Given that no other measures were made, the positional dependence
      // is assumed to be linear for now.
      m_dep_tL.push_back((distanceEL*cm/c_scintilator)/ns + G4RandExponential::shoot(-1./time_spread_slope) * distanceEL/((distanceEL+distanceER)/2));
      m_dep_tR.push_back((distanceER*cm/c_scintilator)/ns + G4RandExponential::shoot(-1./time_spread_slope) * distanceER/((distanceEL+distanceER)/2));

      m_NDeposits++;
      //SetDeposits(m_deposit);

    }
  //track->GetVolume()->GetMotherLogical()->GetName()
}

bool DepositEnergyTTAnaMgr::save_into_data_model() {
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
  
  unsigned int DMChannel;
  int Nnwall(0);
  int Nnplane(0);
  int Nnpanel(0);
  int Nnbar(0);
  
  // save the hits
  for (int i = 0; i < m_NTouchedBar; ++i) {
    for(int hh=0;hh<2;hh++){ //hh=0 Left, hh=1 right
      
      JM::SimTTHit* hit = m_simevent->addTTHit();
      //hit->setBarID( m_TB_channel[i] );//old channel number
      Nnwall=int(m_TB_channel[i]/10000)-1;
      Nnplane=int(int(m_TB_channel[i]/1000)-int(m_TB_channel[i]/10000)*10)-1;
      Nnpanel=int(int(m_TB_channel[i]/100)-int(m_TB_channel[i]/1000)*10)-1;
      Nnbar=int(m_TB_channel[i]-int(m_TB_channel[i]/100)*100)-1;
            
      //int layer = de->geom_info("Wall.Layer",Nnwall);
      //int column = de->geom_info("Wall.Column",Nnwall);
      //int row = de->geom_info("Wall.Row",Nnwall);
      int pmt = TtID::pmt_from_geant(hh, Nnplane, Nnpanel);

      int strip = Nnbar;
      
      DMChannel = TtID::getIntID(Nnwall, pmt, strip);
      
      hit->setChannelID(DMChannel);
      
      if(hh==0)
	{
	  hit->setPE( m_TB_peL[i] );
	  hit->setTime( m_TB_timeL[i] );
	  hit->setADC( m_TB_ADCL[i] );
	}
      else
	{
	  hit->setPE( m_TB_peR[i] );
	  hit->setTime( m_TB_timeR[i] );
	  hit->setADC( m_TB_ADCR[i] );
	}
      
      hit->setX( m_TB_xc[i] );
      hit->setY( m_TB_yc[i] );
      hit->setZ( m_TB_zc[i] );
    }
  }
  
  return true;
}

void DepositEnergyTTAnaMgr::ComputeDigitTT(){
  
  //read deposits and group them per channels using Poisson on the umber of pe
  int ch(0);
  int ntouch(0);
  int found(0);
  
  for(int i=0;i<m_NDeposits;i++)
    {
      ch=10000*(m_dep_nwall[i]+1) + 1000*(m_dep_nplane[i]+1) + 100*(m_dep_npanel[i]+1) + (m_dep_nbar[i]+1);
      found=0;
      
      for(int kk=0;kk<ntouch;kk++)
        {
	  if(ch==m_TB_channel[kk])
            {
	      // Threshold effect for timing calculation approximated at this point
	      // as this is on number of expected PE, rather than the PE observed
	      // after ADC->PE conversion. To take into account the threshold effect
	      // properly a different flow is needed, and accounting for x-talk becomes
	      // non trivial.
	      if(m_TB_peL[kk] < 1./3) m_TB_timeL[kk] = m_dep_t[i] + m_dep_tL[i];
	      m_TB_peL[kk]+=G4Poisson(m_dep_peL[i]);
	      if(m_TB_peR[kk] < 1./3) m_TB_timeR[kk] = m_dep_t[i] + m_dep_tR[i];
	      m_TB_peR[kk]+=G4Poisson(m_dep_peR[i]);

          if(m_isMuonDeposits[i]==1)
            m_TB_isMuonDeposits[kk]=1;
	      
	      found=1;
	      break;
            }
	  
        }
      
      if(found==0)
        {
	  m_TB_channel[ntouch]=ch;
	  m_TB_timeL[ntouch] = m_dep_t[i] + m_dep_tL[i];
	  m_TB_peL[ntouch]=G4Poisson(m_dep_peL[i]);
	  m_TB_timeR[ntouch] = m_dep_t[i] + m_dep_tR[i];
	  m_TB_peR[ntouch]=G4Poisson(m_dep_peR[i]);
	  m_TB_is_ct[ntouch]=0;
      
      if(m_isMuonDeposits[i]==1)
        m_TB_isMuonDeposits[ntouch]=1;
	  ntouch++;
        }
    }

  //cross talk
  int chPM(0);
  double probCT(0);
  int newchannel(0);
  int newchannelcode(0);
  int ntouchloop=ntouch;
  
  int Nnwall(0);
  int Nnplane(0);
  int Nnpanel(0);
  int Nnbar(0);
  
  int npeloop(0);
  
  for(int i=0;i<ntouchloop;i++)
    {
      npeloop=m_TB_peL[i];
      for(int kkk=0;kkk<npeloop;kkk++)//each pe at a time first Left and then Right
        {
	  chPM=int(m_TB_channel[i]-int(m_TB_channel[i]/100)*100)-1;
	  
            probCT=G4UniformRand();
	    
            if(probCT<=0.92496){newchannel=0;}
            if(probCT>0.92496 && probCT<=0.93806){newchannel=CTchannels[chPM][1];}
            if(probCT>0.93806 &&probCT<=0.95116){newchannel=CTchannels[chPM][2];}
            if(probCT>0.95116 && probCT<=0.96426){newchannel=CTchannels[chPM][3];}
            if(probCT>0.96426 && probCT<=0.97736){newchannel=CTchannels[chPM][4];}
            if(probCT>0.97736 && probCT<=0.98046){newchannel=CTchannels[chPM][5];}
            if(probCT>0.98046 && probCT<=0.98356){newchannel=CTchannels[chPM][6];}
            if(probCT>0.98356 && probCT<=0.98666){newchannel=CTchannels[chPM][7];}
            if(probCT>0.98666 && probCT<=0.98976){newchannel=CTchannels[chPM][8];}
            if(probCT>0.98976 && probCT<=0.99106){newchannel=CTchannels[chPM][9];}
            if(probCT>0.99106 && probCT<=0.99236){newchannel=CTchannels[chPM][10];}
            if(probCT>0.99236 && probCT<=0.99366){newchannel=CTchannels[chPM][11];}
            if(probCT>0.99366 && probCT<=0.99496){newchannel=CTchannels[chPM][12];}
            if(probCT>0.99496 && probCT<=0.99559){newchannel=CTchannels[chPM][13];}
            if(probCT>0.99559 && probCT<=0.99622){newchannel=CTchannels[chPM][14];}
            if(probCT>0.99622 && probCT<=0.99685){newchannel=CTchannels[chPM][15];}
            if(probCT>0.99685 && probCT<=0.99748){newchannel=CTchannels[chPM][16];}
            if(probCT>0.99748 && probCT<=0.99811){newchannel=CTchannels[chPM][17];}
            if(probCT>0.99811 && probCT<=0.99874){newchannel=CTchannels[chPM][18];}
            if(probCT>0.99874 && probCT<=0.99937){newchannel=CTchannels[chPM][19];}
            if(probCT>0.99937 && probCT<=1){newchannel=CTchannels[chPM][20];}
	    
            if(newchannel!=0&&newchannel!=99)
            {
	      m_TB_peL[i]-=1;
	      
	      newchannelcode=int(int(m_TB_channel[i]/100)*100)+newchannel+1;
	      
	      found=0;
	      for(int kk=0;kk<ntouch;kk++)
                {
		  
		  if(newchannelcode==m_TB_channel[kk])
                    {
		      m_TB_peL[kk]+=1;
		      
		      found=1;
		      break;
                    }
                }
	      
	      if(found==0)
                {
		  m_TB_channel[ntouch]=newchannelcode;
		  m_TB_timeL[ntouch] = m_TB_timeL[i]; // not exactly... but better than nothing...
		  m_TB_peL[ntouch]=1;
		  m_TB_timeR[ntouch] = m_TB_timeR[i]; // not exactly... but better than nothing...
		  m_TB_peR[ntouch]=0;  //new
		  m_TB_is_ct[ntouch]=1;
		  ntouch++;
                }
            }
	    
        }
      
        //NOW RIGHT...
              
      npeloop=m_TB_peR[i];
      for(int kkk=0;kkk<npeloop;kkk++)
        {
	  chPM=int(m_TB_channel[i]-int(m_TB_channel[i]/100)*100)-1;

	  probCT=G4UniformRand();
          
	  if(probCT<=0.92496){newchannel=0;}
	  if(probCT>0.92496 && probCT<=0.93806){newchannel=CTchannels[chPM][1];}
	  if(probCT>0.93806 &&probCT<=0.95116){newchannel=CTchannels[chPM][2];}
	  if(probCT>0.95116 && probCT<=0.96426){newchannel=CTchannels[chPM][3];}
	  if(probCT>0.96426 && probCT<=0.97736){newchannel=CTchannels[chPM][4];}
	  if(probCT>0.97736 && probCT<=0.98046){newchannel=CTchannels[chPM][5];}
	  if(probCT>0.98046 && probCT<=0.98356){newchannel=CTchannels[chPM][6];}
	  if(probCT>0.98356 && probCT<=0.98666){newchannel=CTchannels[chPM][7];}
	  if(probCT>0.98666 && probCT<=0.98976){newchannel=CTchannels[chPM][8];}
	  if(probCT>0.98976 && probCT<=0.99106){newchannel=CTchannels[chPM][9];}
	  if(probCT>0.99106 && probCT<=0.99236){newchannel=CTchannels[chPM][10];}
	  if(probCT>0.99236 && probCT<=0.99366){newchannel=CTchannels[chPM][11];}
	  if(probCT>0.99366 && probCT<=0.99496){newchannel=CTchannels[chPM][12];}
	  if(probCT>0.99496 && probCT<=0.99559){newchannel=CTchannels[chPM][13];}
	  if(probCT>0.99559 && probCT<=0.99622){newchannel=CTchannels[chPM][14];}
	  if(probCT>0.99622 && probCT<=0.99685){newchannel=CTchannels[chPM][15];}
	  if(probCT>0.99685 && probCT<=0.99748){newchannel=CTchannels[chPM][16];}
	  if(probCT>0.99748 && probCT<=0.99811){newchannel=CTchannels[chPM][17];}
	  if(probCT>0.99811 && probCT<=0.99874){newchannel=CTchannels[chPM][18];}
	  if(probCT>0.99874 && probCT<=0.99937){newchannel=CTchannels[chPM][19];}
	  if(probCT>0.99937 && probCT<=1){newchannel=CTchannels[chPM][20];}
	  
	  if(newchannel!=0&&newchannel!=99)
            {	      
	      m_TB_peR[i]-=1;
	      
	      newchannelcode=int(int(m_TB_channel[i]/100)*100)+newchannel+1;
	      
	      found=0;
	      for(int kk=0;kk<ntouch;kk++)
                {		  
		  if(newchannelcode==m_TB_channel[kk])
                    {
		      m_TB_peR[kk]+=1;
                      
		      found=1;
		      break;
                    }
                }
	      
	      if(found==0)
                {
		  m_TB_channel[ntouch]=newchannelcode;
		  m_TB_timeR[ntouch] = m_TB_timeR[i]; // not exactly... but better than nothing...
		  m_TB_peR[ntouch]=1;
		  m_TB_timeL[ntouch] = m_TB_timeL[i]; // not exactly... but better than nothing...
		  m_TB_peL[ntouch]=0; //new
		  m_TB_is_ct[ntouch]=1;
		  ntouch++;
                }
            }
	  
        }
      
    }
  
  //write bar position, ADC with mean gain 40 AC units and 50% resolution
  //drop events not triggered (1/3 p.e. in OR between left and right)
  
  int dropped(0);
  
  for(int i=0;i<ntouch;i++)
    {
      
      Nnwall=int(m_TB_channel[i]/10000)-1;
      Nnplane=int(int(m_TB_channel[i]/1000)-int(m_TB_channel[i]/10000)*10)-1;
      Nnpanel=int(int(m_TB_channel[i]/100)-int(m_TB_channel[i]/1000)*10)-1;
      Nnbar=int(m_TB_channel[i]-int(m_TB_channel[i]/100)*100)-1;
      
      //get bar coordinates
      double x0wall = de->geom_info("Wall.X",Nnwall);
      double y0wall = de->geom_info("Wall.Y",Nnwall);
      double z0wall = de->geom_info("Wall.Z",Nnwall);
      double xy0pan = de->geom_info("Panel.XY",Nnpanel);
      double xy0bar = de->geom_info("Bar.XY",Nnbar);
      double z0plane = de->geom_info("Plane.Z",Nnplane);
      
      if(Nnplane==1)
        {
	  m_TB_yc[i-dropped]=xy0bar+xy0pan+y0wall;
	  m_TB_xc[i-dropped]=1E9;
        }
      else
        {
	  m_TB_xc[i-dropped]=xy0bar+xy0pan+x0wall;
	  m_TB_yc[i-dropped]=1E9;
        }
      
      m_TB_zc[i-dropped]=z0plane+z0wall;
      
      
      double meanGain(40); //mean gain in ADC units
      
      m_TB_ADCL[i-dropped]=G4RandGauss::shoot(m_TB_peL[i]*meanGain,m_TB_peL[i]*0.5*meanGain);
      if(m_TB_ADCL[i-dropped]<0)
	m_TB_ADCL[i-dropped]=0;
      if(m_TB_ADCL[i-dropped]>4096)
	m_TB_ADCL[i-dropped]=4096;
      
      
      m_TB_ADCR[i-dropped]=G4RandGauss::shoot(m_TB_peR[i]*meanGain,m_TB_peR[i]*0.5*meanGain);
      if(m_TB_ADCR[i-dropped]<0)
	m_TB_ADCR[i-dropped]=0;
      if(m_TB_ADCR[i-dropped]>4096)
	m_TB_ADCR[i-dropped]=4096;
      
      m_TB_peL[i-dropped]= m_TB_ADCL[i-dropped]/meanGain;
      m_TB_peR[i-dropped]= m_TB_ADCR[i-dropped]/meanGain;
      
      m_TB_channel[i-dropped]=m_TB_channel[i];
      m_TB_is_ct[i-dropped]=m_TB_is_ct[i];
      
      m_TB_isMuonDeposits[i-dropped] = m_TB_isMuonDeposits[i];
      
      if( !( m_TB_peL[i-dropped]>0.34  ||  m_TB_peR[i-dropped]>0.34 ) )
	dropped++;
    }
  
  
  m_NTouchedBar=ntouch-dropped;
  
  
  for(int jj=0;jj<ntouch-dropped;jj++)
    {      
      //convert left and right bar into channel info
      
      Nnwall=int(m_TB_channel[jj]/10000)-1;
      Nnplane=int(int(m_TB_channel[jj]/1000)-int(m_TB_channel[jj]/10000)*10)-1;
      Nnpanel=int(int(m_TB_channel[jj]/100)-int(m_TB_channel[jj]/1000)*10)-1;
      Nnbar=int(m_TB_channel[jj]-int(m_TB_channel[jj]/100)*100)-1;
      
      //int layer = de->geom_info("Wall.Layer",Nnwall);
      //int column = de->geom_info("Wall.Column",Nnwall);
      //int row = de->geom_info("Wall.Row",Nnwall);
      int strip = Nnbar;
      
      if(m_TB_peL[jj]>0)
	{
	  
	  int pmt = TtID::pmt_from_geant(0, Nnplane, Nnpanel);
	  
	  m_TB_channelC[m_NTouchedChannel]=m_TB_channel[jj];
	  m_TB_pe[m_NTouchedChannel]=m_TB_peL[jj];
	  m_TB_time[m_NTouchedChannel]=m_TB_timeL[jj];
	  m_TB_ADC[m_NTouchedChannel]=m_TB_ADCL[jj];
	  m_TB_xcC[m_NTouchedChannel]=m_TB_xc[jj];
	  m_TB_ycC[m_NTouchedChannel]=m_TB_yc[jj];
	  m_TB_zcC[m_NTouchedChannel]=m_TB_zc[jj];
	  m_TB_is_ctC[m_NTouchedChannel]=m_TB_is_ct[jj];
      m_TB_isMuonDepositsC[m_NTouchedChannel]=m_TB_isMuonDeposits[jj];
	  m_TB_DMchannel[m_NTouchedChannel]=TtID::getIntID(Nnwall, pmt, strip);
	  m_NTouchedChannel++;
        }
      if(m_TB_peR[jj]>0)
        {
	  int pmt = TtID::pmt_from_geant(1, Nnplane, Nnpanel);
	  	  
	  m_TB_channelC[m_NTouchedChannel]=m_TB_channel[jj];
	  m_TB_pe[m_NTouchedChannel]=m_TB_peR[jj];
	  m_TB_time[m_NTouchedChannel]=m_TB_timeR[jj];
	  m_TB_ADC[m_NTouchedChannel]=m_TB_ADCR[jj];
	  m_TB_xcC[m_NTouchedChannel]=m_TB_xc[jj];
	  m_TB_ycC[m_NTouchedChannel]=m_TB_yc[jj];
	  m_TB_zcC[m_NTouchedChannel]=m_TB_zc[jj];
	  m_TB_is_ctC[m_NTouchedChannel]=m_TB_is_ct[jj];
      m_TB_isMuonDepositsC[m_NTouchedChannel]=m_TB_isMuonDeposits[jj];
	  m_TB_DMchannel[m_NTouchedChannel]=TtID::getIntID(Nnwall, pmt, strip);
	  m_NTouchedChannel++;
        }      
    }

  if (m_flag_ntuple and m_evt_treeTTDigit) {
    m_evt_treeTTDigit -> Fill();
  }
  
}

/* vim: set noexpandtab shiftwidth=2 tabstop=8 softtabstop=2: */
