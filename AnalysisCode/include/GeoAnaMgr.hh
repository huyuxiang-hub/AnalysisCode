#ifndef GeoAnaMgr_hh
#define GeoAnaMgr_hh

#include "SniperKernel/ToolBase.h"
#include "DetSimAlg/IAnalysisElement.h"
#include <string>

class TGeoManager;

class GeoAnaMgr: public IAnalysisElement,
                 public ToolBase{
public:
    GeoAnaMgr(const std::string& name);
    ~GeoAnaMgr();
    // Run Action
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);
    // Event Action
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

private:
    void saveGeom();
    void deleteGdmlFile();

    void save_use_rootio(TGeoManager* geom);
private:
    std::string m_filename;
    bool m_flag_deletegdml;
    bool m_flag_gdml2root;
private:
    bool m_flag_output;
    bool m_flag_first;
    bool m_store_refs; // Save 0x... in GDML
#ifdef G4DAE
    bool m_flag_dae;
    std::string m_filename_dae;
#endif
};

#endif
