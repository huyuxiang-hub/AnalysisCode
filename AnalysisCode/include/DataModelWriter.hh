#ifndef DataModelWriter_hh
#define DataModelWriter_hh

#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IAnalysisElement.h"

#include "junoHit_PMT.hh"
#include "junoHit_PMT_muon.hh"
#include "Event/SimHeader.h"
#include "Event/SimEvent.h"

class ISimTrackSvc;
class G4Event;
namespace JM {
    class SimEvent;
}

template<typename T>
struct Holder {
    const static int is_normal_hit = 1;
    const static int is_muon_hit = 0;
};
template<> struct Holder<junoHit_PMT_muon> {
    const static int is_normal_hit = 0;
    const static int is_muon_hit = 1;
};


class DataModelWriter: public IAnalysisElement,
                       public ToolBase{

public:
    DataModelWriter(const std::string& name);
    ~DataModelWriter();

    // Run Action
    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);
    // Event Action
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);
private:
    int m_nPhotons;
    double m_timewindow;
    int m_nPhotons_muon;
    double m_timewindow_muon;

    ISimTrackSvc* simtracksvc;

    // fill hits
    void fill_hits(JM::SimEvent* dst, const G4Event* evt);
    template<typename T>
    void fill_hits_tmpl(G4THitsCollection<T>* col, JM::SimEvent* dst) {
        if (col) {
   
            int n_hit = col->entries();

            if (Holder<T>::is_normal_hit==1) {
                m_nPhotons = n_hit;
                m_timewindow = 0.0;
            } else if (Holder<T>::is_muon_hit==1) {
                m_nPhotons_muon = 0;
                m_timewindow_muon = 0.0;
            }

            bool minmax_initialized = false;
            double max_CDLPMT_hittime = 0;
            double min_CDLPMT_hittime = 0;
            for (int i = 0; i < n_hit; ++i) {
                // create new hit
                // The PMT Hit can be from WP (Water Pool) or CD (Central
                // Detector). 
                // Please use the copy no to distinguish the PMT.
                int copyno = (*col)[i]->GetPMTID();
                JM::SimPMTHit* jm_hit = 0;
                // FIXME: hard code the copy no
                if ((copyno < 30000) or (copyno >= 300000)) {
                    // TODO because in current Data Model, the 3inch and the 20inch
                    // PMTs are in the same collection.
                    jm_hit = dst->addCDHit();
                } else if (copyno >= 30000) {
                    jm_hit = dst->addWPHit();
                }
                if(copyno < 30000) 
                {
                    if (!minmax_initialized) { // initialize for the first hit
                        max_CDLPMT_hittime = ((*col)[i]->GetTime());
                        min_CDLPMT_hittime = ((*col)[i]->GetTime());
                        minmax_initialized = true;
                    } else {
                        if ((*col)[i]->GetTime() < min_CDLPMT_hittime) {
                            min_CDLPMT_hittime = (*col)[i]->GetTime();
                        }
                        if ((*col)[i]->GetTime() > max_CDLPMT_hittime) {
                            max_CDLPMT_hittime = (*col)[i]->GetTime();
                        }
                    }
                }
                jm_hit->setPMTID( (*col)[i]->GetPMTID() );
                jm_hit->setNPE( (*col)[i]->GetCount() );
                jm_hit->setHitTime( (*col)[i]->GetTime() );
                jm_hit->setTrackID( (*col)[i]->GetProducerID() );
                jm_hit->setLocalTheta( (*col)[i]->GetTheta() );
                jm_hit->setLocalPhi( (*col)[i]->GetPhi() );

                if (Holder<T>::is_muon_hit==1 && (copyno < 30000 || copyno >= 300000)) {
                    m_nPhotons_muon += (*col)[i]->GetCount();
                }
            }

            if (Holder<T>::is_normal_hit==1) {
                m_timewindow = max_CDLPMT_hittime - min_CDLPMT_hittime;
            } else if (Holder<T>::is_muon_hit==1) {
                m_timewindow_muon = max_CDLPMT_hittime - min_CDLPMT_hittime;
            }

        }
    }
    // fill tracks
    void fill_tracks(JM::SimEvent* dst, const G4Event* evt);
    
    // get the additional tracks from SimTrackSvc.
    void fill_additional_tracks(JM::SimEvent* dst);
};

#endif
