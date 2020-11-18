/*
 * MuTrackingAnaMgr was created to make it possible to save the full muon
 * path in the simulation, taking into account deviations from the generated
 * MC direction.
 *
 * The first and last steps of the muon are always saved, so in all cases
 * you are guaranteed to have at least 2 sets of position and direction
 * for each muon. For all saved steps the 3D position and momentum are
 * saved.
 *
 * Since the goal is to track deviations that would be significant in the
 * context of JUNO, a minimal deviation threshold is defined, and if the
 * 'kink' in the muon tracjectory is smaller than that we won't save that
 * information. This is set by default to 1 mrad (ie about 0.06 degrees).
 * As the kinks are always compared to the last recorded information,
 * multiple small kinks eventually show up as one entry adding them up
 * together once the sum goes above the set threshold.
 * The threshold can be changed by setting the DeviationThreshold property.
 *
 * In addition to a threshold to saving kinks, a momentum threshold was
 * also added to avoid saving low-energy secondary muons produced at
 * low-energy. Primary muons are always tracked independently of their
 * energy.
 * The momentum threshold is by default set to 0 MeV and can be changed
 * by changing the MomentumThreshold property.
 *
 * Note: If no information on secondary muons is desired, they can be
 * removed by setting a very high MomentumThreshold.
 *
 * To enable the MuTrackingAnaMgr when using tut_detsim.py please add to
 * the command line the option '--anamgr-list MuTrackingAnaMgr' or add it
 * to an anamgr-config-file.
 *
 * Enabling the MuTrackingAnaMgr will add a 'mu_tracking' TTree to the
 * user ROOT file. It is useful to know that each entry of the TTree will
 * correspond to a unique (evtID, MuTrackID) pair. This means that the
 * entry number in this TTree will not correspond to the entry number on
 * other TTrees produce (like the geninfo TTree), so do to a matching of
 * the information the evtID needs to be used. This was done this way to
 * make it possible to track an a priori unknown number of muons through an
 * a priori unknown number of steps, while keeping things sufficiently
 * simple.
 *
 * Author: J. P. A. M. de Andre <jpandre+junosvn@iphc.cnrs.fr>
 *
 */

#ifndef MuTrackingAnaMgr_hh
#define MuTrackingAnaMgr_hh

#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IAnalysisElement.h"

#include <TTree.h>
#include <vector>
#include <map>
#include <set>

class MuTrackingAnaMgr: public IAnalysisElement, public ToolBase{
    public:
        MuTrackingAnaMgr(const std::string& name);
        ~MuTrackingAnaMgr();

        // Run Action
        virtual void BeginOfRunAction(const G4Run*);
        virtual void EndOfRunAction(const G4Run*);

        // Event Action
        virtual void BeginOfEventAction(const G4Event*);
        virtual void EndOfEventAction(const G4Event*);

        virtual void UserSteppingAction(const G4Step* step);

    private:
        // Parameters
        float deviation_threshold;
        float momentum_threshold;

        // The output muon tracking tree
        TTree*   m_muon_track_tree;
        int m_eventID;
        int m_MuTrackID;
        int m_MuParentID;
        std::vector<float> * m_MuPosx;
        std::vector<float> * m_MuPosy;
        std::vector<float> * m_MuPosz;
        std::vector<float> * m_MuPx  ;
        std::vector<float> * m_MuPy  ;
        std::vector<float> * m_MuPz  ;

        // Storage needed during processing
        std::vector<int> m_v_MuTrackID  ;
        std::vector<int> m_v_MuParentID ;
        std::map<int, std::vector<float> *> m_v_MuPosx;
        std::map<int, std::vector<float> *> m_v_MuPosy;
        std::map<int, std::vector<float> *> m_v_MuPosz;
        std::map<int, std::vector<float> *> m_v_MuPx  ;
        std::map<int, std::vector<float> *> m_v_MuPy  ;
        std::map<int, std::vector<float> *> m_v_MuPz  ;
        std::map<int, float> m_v_MuExitPosx;
        std::map<int, float> m_v_MuExitPosy;
        std::map<int, float> m_v_MuExitPosz;
        std::map<int, float> m_v_MuExitPx  ;
        std::map<int, float> m_v_MuExitPy  ;
        std::map<int, float> m_v_MuExitPz  ;

        // List of low momentum tracks to skip them afterwards
        std::set<int> m_trackID_lowP;
};

#endif // MuTrackingAnaMgr_hh

// vim: et:sts=4:ts=4:sw=4
