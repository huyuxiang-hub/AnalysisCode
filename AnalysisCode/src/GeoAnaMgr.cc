
#include <boost/python.hpp>
#include "GeoAnaMgr.hh"
#include "G4Run.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperLog.h"
#include "SniperKernel/SniperException.h"
#include "RootWriter/RootWriter.h"

#include "G4GDMLParser.hh"

#include <fstream>
#include "TGeoManager.h"

#include <boost/filesystem.hpp>
#include <exception>

#ifdef G4DAE
#include "G4DAEParser.hh"
#endif

DECLARE_TOOL(GeoAnaMgr);

GeoAnaMgr::GeoAnaMgr(const std::string& name) 
    : ToolBase(name) {
    declProp("GdmlEnable", m_flag_output=false);
    declProp("GdmlOutput", m_filename="test.gdml");
    declProp("GdmlStoreRefs", m_store_refs=true);
    declProp("Gdml2RootEnable", m_flag_gdml2root=false);
    declProp("GdmlDeleteFile", m_flag_deletegdml=false);
#ifdef G4DAE
    declProp("DAEEnable", m_flag_dae=false);
    declProp("DAEOutput", m_filename_dae="test.dae");
#endif

    m_flag_first = true;
}

GeoAnaMgr::~GeoAnaMgr() {

}

void
GeoAnaMgr::BeginOfRunAction(const G4Run* /*run*/) {
}

void
GeoAnaMgr::EndOfRunAction(const G4Run* /*run*/) {

}

void
GeoAnaMgr::BeginOfEventAction(const G4Event*) {
    // get the world
    //
    if (not m_flag_first) {
        return;
    }
    m_flag_first = false;
    G4Navigator* gNavigator =
        G4TransportationManager::GetTransportationManager()
        ->GetNavigatorForTracking();

    G4VPhysicalVolume* pworld = gNavigator->GetWorldVolume();
    G4LogicalVolume* lworld = pworld->GetLogicalVolume();


    if (m_flag_output) {
        G4GDMLReadStructure* fReader = new G4GDMLReadStructure;
        G4GDMLWriteStructure* fWriter = new G4GDMLWriteStructure;
        G4GDMLParser* fParser = new G4GDMLParser(fReader, fWriter);

        // fParser->Write(m_filename, lworld);

        // FIXME: Convert GDML to ROOT

        // TODO: 
        // To make sure the jobs run in parallel don't crash, 
        // we save the result into a temporary file, then use 
        // root to convert this temporary into ROOT object.
        std::string real_filename = m_filename;
        boost::filesystem::path temp = boost::filesystem::unique_path();
        const std::string tempstr    = temp.native();
        m_filename = tempstr+".gdml";
        fParser->Write(m_filename, lworld, m_store_refs);
        saveGeom();
        // * copy the temp file to the user want?
        if (!m_flag_deletegdml) {
            try {
                boost::filesystem::copy_file( m_filename, real_filename );
            } catch(const std::runtime_error& e) {
                LogWarn << e.what() << std::endl;
            }
        }
        // * remove the temp one
        // force remove the old one
        deleteGdmlFile();

        m_filename = real_filename;
    }

#ifdef G4DAE
    if (m_flag_dae) {
        G4DAEParser* fParserDAE = new G4DAEParser;

        fParserDAE->Write(m_filename_dae, pworld, true, false, -1);

    }
#endif

}

void
GeoAnaMgr::EndOfEventAction(const G4Event*) {

}

void
GeoAnaMgr::saveGeom() {
    if (!m_flag_gdml2root) return;
    std::streambuf *backup;
    std::ofstream fout;
    //fout.open("data.out");
    fout.open("/dev/null");
    backup = std::cout.rdbuf();    // back up cout's streambuf
    std::cout.rdbuf(fout.rdbuf()); // assign file's streambuf to cout
    TGeoManager* geom = TGeoManager::Import(m_filename.c_str(), "gdml geom");
    geom->SetName("JunoGeom");
    std::cout.rdbuf(backup);     // restore cin's original streambuf

    // TODO: save the geom to where?
    // or using the RootIO service??
    geom->Write();

    // * using RootIO
    save_use_rootio(geom);
}

void
GeoAnaMgr::deleteGdmlFile() {
    std::ifstream ifile(m_filename.c_str());
    if (ifile && remove(m_filename.c_str()) != 0 ) {
        perror( "Error deleting gdml file" );
    }
    else {
        std::cout << "Gdml file " << m_filename << " deleted" << std::endl;
    }

}

#include "RootIOSvc/RootOutputSvc.h"
void
GeoAnaMgr::save_use_rootio(TGeoManager* geom) {
    SniperPtr<RootOutputSvc> output_svc(*getParent(), "OutputSvc");
    if (output_svc.invalid()) {
        LogWarn << "Can't find the OutputSvc in current task." << std::endl;
        return;
    }
    bool status = output_svc->attachObj("/Event/Sim", geom);
    LogDebug << "Attach Geometry: " << status << std::endl;
}
