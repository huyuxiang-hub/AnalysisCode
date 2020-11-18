#ifndef MuFastnProcessAnaMgr_hh
#define MuFastnProcessAnaMgr_hh

#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IAnalysisElement.h"
#include "TTree.h"

class MuFastnProcessAnaMgr: public IAnalysisElement,
                               public ToolBase{
public:

    MuFastnProcessAnaMgr(const std::string& name);
    ~MuFastnProcessAnaMgr();
    // Run Action
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
    // Event Action
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

    virtual void PreUserTrackingAction(const G4Track* aTrack);
    virtual void PostUserTrackingAction(const G4Track* aTrack);

    virtual void UserSteppingAction(const G4Step* step);

    //heap sort++++++++++++++++++++++++++++
    void shift(double a[], double loct[], double capt[],
               double deposit[], double qdeposit[], double ke[], double totale[], 
               int pdg[],int parid[],
               double posx[],double posy[],double posz[], int i , int m)
    {
      int k;
      double t, tloct, tcapt, tdeposit, tqdeposit, tke, ttotale, tx, ty, tz;
      int tpdg, tparid;
      t = a[i]; 
      k = 2 * i + 1;
      tloct = loct[i]; tcapt = capt[i];
      tdeposit = deposit[i]; tqdeposit = qdeposit[i]; tke = ke[i]; ttotale = totale[i];
      tpdg = pdg[i]; tparid = parid[i];
      tx = posx[i]; ty = posy[i]; tz = posz[i];
      while (k < m)
      {
        if ((k < m - 1) && (a[k] < a[k+1])) k ++;
        if (t < a[k]) {
          a[i] = a[k];
          loct[i] = loct[k]; capt[i] = capt[k];
          deposit[i] = deposit[k]; qdeposit[i] = qdeposit[k]; ke[i]=ke[k]; totale[i]=totale[k];
          pdg[i] = pdg[k]; parid[i] = parid[k];
          posx[i] = posx[k]; posy[i] = posy[k]; posz[i] = posz[k];
          i = k; k = 2 * i + 1;
        }
        else break;
      }
      a[i] = t;
      loct[i] = tloct; capt[i] = tcapt; 
      deposit[i] = tdeposit; qdeposit[i] = tqdeposit; ke[i] = tke; totale[i] = ttotale;
      pdg[i] = tpdg; parid[i] = tparid;
      posx[i] = tx; posy[i] = ty; posz[i] = tz;
    }

    void heap(double a[], double loct[], double capt[], 
              double deposit[], double qdeposit[], double ke[], double totale[], 
              int pdg[], int parid[], 
              double posx[],double posy[],double posz[], int n){ 
      int i ;
      double k, kloct, kcapt, kdeposit, kqdeposit, kke, ktotale, kpos;
      int kpdg, kparid;
      for (i = n/2-1; i >= 0; i --){
        shift(a ,loct, capt, deposit, qdeposit, ke, totale, pdg, parid, posx, posy, posz, i, n);
      }
      for (i = n-1; i >= 1; i --)
      {
        k = a[0]; a[0] = a[i]; a[i] = k;
        kloct = loct[0]; loct[0] = loct[i]; loct[i] = kloct;
        kcapt = capt[0]; capt[0] = capt[i]; capt[i] = kcapt;
        kdeposit = deposit[0]; deposit[0] = deposit[i]; deposit[i] = kdeposit;
        kqdeposit = qdeposit[0]; qdeposit[0] = qdeposit[i]; qdeposit[i] = kqdeposit;
        kke = ke[0]; ke[0] = ke[i]; ke[i] = kke;
        ktotale = totale[0]; totale[0] = totale[i]; totale[i] = ktotale;
        kpdg = pdg[0]; pdg[0] = pdg[i]; pdg[i] = kpdg;
        kparid = parid[0]; parid[0] = parid[i]; parid[i] = kparid;
        kpos = posx[0]; posx[0] = posx[i]; posx[i] = kpos;
        kpos = posy[0]; posy[0] = posy[i]; posy[i] = kpos;
        kpos = posz[0]; posz[0] = posz[i]; posz[i] = kpos;
        shift(a, loct, capt, deposit, qdeposit, ke, totale, pdg ,parid, posx, posy, posz, 0, i);
      }
    }


private:
    // Evt Data
    int m_eventID;

    TTree*   m_fastn_tree;
    bool mMuonFlag;
    //step tree for fast neutron
    static const int maxStepNumber = 700000;
    int    mStepNumber;
    int    mfnPDG[maxStepNumber];
    int    mfnParId[maxStepNumber];
    double mfnPositionx[maxStepNumber];
    double mfnPositiony[maxStepNumber];
    double mfnPositionz[maxStepNumber];
    double mfnGlobalTime[maxStepNumber];
    double mfnLocalTime[maxStepNumber];
    double mfnTotalEnergy[maxStepNumber];
    double mfnKineticEnergy[maxStepNumber];
    double mfnEnergyDeposit[maxStepNumber];
    double mfnQEnergyDeposit[maxStepNumber];
    double mfnCaptime[maxStepNumber];

};

#endif
