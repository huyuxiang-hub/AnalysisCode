#ifndef OpticalParameterAnaMgr_hh
#define OpticalParameterAnaMgr_hh

#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IAnalysisElement.h"
#include <string>
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"

class TTree;

/*  
 * This analysis element is used to dump optical properties
 * defined in 
 * * http://juno.ihep.ac.cn/mediawiki/index.php/Analysis:Basic_Distributions_of_JUNO
 **/

class OpticalParameterAnaMgr: public IAnalysisElement,
                              public ToolBase {
public:
    OpticalParameterAnaMgr(const std::string& name);
    ~OpticalParameterAnaMgr();

    // Run Action
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
    // Event Action
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
private:
    // return N, x[N], y[N]
    void get_matprop(const G4MaterialPropertyVector*, int& N, double*& x, double*& y);
    void get_matprop(G4MaterialPropertiesTable*, const char*, int& N, double*& x, double*& y);
private:
    TTree* m_op_tree;
};

#endif
