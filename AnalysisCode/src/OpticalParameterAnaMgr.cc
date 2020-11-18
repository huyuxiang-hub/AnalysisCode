
#include "OpticalParameterAnaMgr.hh"
#include "G4Run.hh"

#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/ToolFactory.h"
#include "SniperKernel/SniperLog.h"
#include "SniperKernel/SniperException.h"
#include "RootWriter/RootWriter.h"
#include "TTree.h"

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"

DECLARE_TOOL(OpticalParameterAnaMgr);

OpticalParameterAnaMgr::OpticalParameterAnaMgr(const std::string& name)
    : ToolBase(name) {
    m_op_tree = 0;
}

OpticalParameterAnaMgr::~OpticalParameterAnaMgr() {

}

void
OpticalParameterAnaMgr::BeginOfRunAction(const G4Run* /*run*/) {
    SniperPtr<RootWriter> svc(*getParent(), "RootWriter");
    if (svc.invalid()) {
        LogError << "Can't Locate RootWriter. If you want to use it, please "
                 << "enalbe it in your job option file."
                 << std::endl;
        return;
    }

    m_op_tree = svc->bookTree("SIMEVT/opticalparam", "Optical Parameters");

    // NOTE: 
    // we use geant4 9.4 now, so the property vector is using the old style.
    // If we use new version, this code needs to be upgraded.
    // -- 2015.10.14 Tao Lin
    // ======================================================================
    // * LS light yield, refractive index. 
    // * LS emission spectrum, re-emission prob. 
    // * LS absorption length and Rayleigh scattering length. 
    // ======================================================================
    G4Material* mat_LS = G4Material::GetMaterial("LS");
    G4int LS_LY_n = 0;         G4double* LS_LY_energy = 0;         G4double* LS_LY_ly = 0;
    G4int LS_RI_n = 0;         G4double* LS_RI_energy = 0;         G4double* LS_RI_idx = 0;
    G4int LS_Emission_n = 0;   G4double* LS_Emission_energy = 0;   G4double* LS_Emission_spec = 0;
    G4int LS_ReEmission_n = 0; G4double* LS_ReEmission_energy = 0; G4double* LS_ReEmission_prob = 0;
    G4int LS_ABS_n = 0;        G4double* LS_ABS_energy = 0;        G4double* LS_ABS_len = 0;
    G4int LS_Rayleigh_n = 0;   G4double* LS_Rayleigh_energy = 0;   G4double* LS_Rayleigh_len = 0;
    if (mat_LS) {
        G4MaterialPropertiesTable* tbl_LS = mat_LS->GetMaterialPropertiesTable();
        // SCINTILLATIONYIELD
        get_matprop(tbl_LS, "SCINTILLATIONYIELD", LS_LY_n, LS_LY_energy, LS_LY_ly);
        // RINDEX
        get_matprop(tbl_LS, "RINDEX", LS_RI_n, LS_RI_energy, LS_RI_idx);
        // emission including Fast and Slow, here we save Fast (TODO)
        get_matprop(tbl_LS, "FASTCOMPONENT", LS_Emission_n, LS_Emission_energy, LS_Emission_spec);
        // REEMISSIONPROB
        get_matprop(tbl_LS, "REEMISSIONPROB", LS_ReEmission_n, LS_ReEmission_energy, LS_ReEmission_prob);
        // ABSLENGTH
        get_matprop(tbl_LS, "ABSLENGTH", LS_ABS_n, LS_ABS_energy, LS_ABS_len);
        // RAYLEIGH
        get_matprop(tbl_LS, "RAYLEIGH", LS_Rayleigh_n, LS_Rayleigh_energy, LS_Rayleigh_len);
    }
    m_op_tree->Branch("LS_LY_n", &LS_LY_n, "LS_LY_n/I");
    m_op_tree->Branch("LS_LY_energy", LS_LY_energy, "LS_LY_energy[LS_LY_n]/D");
    m_op_tree->Branch("LS_LY_ly", LS_LY_ly, "LS_LY_ly[LS_LY_n]/D");

    m_op_tree->Branch("LS_RI_n", &LS_RI_n, "LS_RI_n/I");
    m_op_tree->Branch("LS_RI_energy", LS_RI_energy, "LS_RI_energy[LS_RI_n]/D");
    m_op_tree->Branch("LS_RI_idx", LS_RI_idx, "LS_RI_idx[LS_RI_n]/D");

    m_op_tree->Branch("LS_Emission_n", &LS_Emission_n, "LS_Emission_n/I");
    m_op_tree->Branch("LS_Emission_energy", LS_Emission_energy, "LS_Emission_energy[LS_Emission_n]/D");
    m_op_tree->Branch("LS_Emission_spec", LS_Emission_spec, "LS_Emission_spec[LS_Emission_n]/D");

    m_op_tree->Branch("LS_ReEmission_n", &LS_ReEmission_n, "LS_ReEmission_n/I");
    m_op_tree->Branch("LS_ReEmission_energy", LS_ReEmission_energy, "LS_ReEmission_energy[LS_ReEmission_n]/D");
    m_op_tree->Branch("LS_ReEmission_prob", LS_ReEmission_prob, "LS_ReEmission_prob[LS_ReEmission_n]/D");

    m_op_tree->Branch("LS_ABS_n", &LS_ABS_n, "LS_ABS_n/I");
    m_op_tree->Branch("LS_ABS_energy", LS_ABS_energy, "LS_ABS_energy[LS_ABS_n]/D");
    m_op_tree->Branch("LS_ABS_len", LS_ABS_len, "LS_ABS_len[LS_ABS_n]/D");

    m_op_tree->Branch("LS_Rayleigh_n", &LS_Rayleigh_n, "LS_Rayleigh_n/I");
    m_op_tree->Branch("LS_Rayleigh_energy", LS_Rayleigh_energy, "LS_Rayleigh_energy[LS_Rayleigh_n]/D");
    m_op_tree->Branch("LS_Rayleigh_len", LS_Rayleigh_len, "LS_Rayleigh_len[LS_Rayleigh_n]/D");

    // ======================================================================
    // * Acrylic absorption length and refractive index.
    // ======================================================================
    G4Material* mat_Acrylic = G4Material::GetMaterial("Acrylic");
    G4int Acrylic_ABS_n = 0; G4double* Acrylic_ABS_energy = 0; G4double* Acrylic_ABS_len = 0;
    G4int Acrylic_RI_n = 0;  G4double* Acrylic_RI_energy = 0;  G4double* Acrylic_RI_idx = 0;
    if (mat_Acrylic) {
        G4MaterialPropertiesTable* tbl_Acrylic = mat_Acrylic->GetMaterialPropertiesTable();
        // ABSLENGTH
        get_matprop(tbl_Acrylic, "ABSLENGTH", Acrylic_ABS_n, Acrylic_ABS_energy, Acrylic_ABS_len);
        // RINDEX
        get_matprop(tbl_Acrylic, "RINDEX", Acrylic_RI_n, Acrylic_RI_energy, Acrylic_RI_idx);
    }

    m_op_tree->Branch("Acrylic_ABS_n", &Acrylic_ABS_n, "Acrylic_ABS_n/I");
    m_op_tree->Branch("Acrylic_ABS_energy", Acrylic_ABS_energy, "Acrylic_ABS_energy[Acrylic_ABS_n]/D");
    m_op_tree->Branch("Acrylic_ABS_len", Acrylic_ABS_len, "Acrylic_ABS_len[Acrylic_ABS_n]/D");

    m_op_tree->Branch("Acrylic_RI_n", &Acrylic_RI_n, "Acrylic_RI_n/I");
    m_op_tree->Branch("Acrylic_RI_energy", Acrylic_RI_energy, "Acrylic_RI_energy[Acrylic_RI_n]/D");
    m_op_tree->Branch("Acrylic_RI_idx", Acrylic_RI_idx, "Acrylic_RI_idx[Acrylic_RI_n]/D");

    // ======================================================================
    // * Water absorption length and refractive index.
    // ======================================================================
    G4Material* mat_Water = G4Material::GetMaterial("Water");
    G4int Water_ABS_n = 0; G4double* Water_ABS_energy = 0; G4double* Water_ABS_len = 0;
    G4int Water_RI_n = 0;  G4double* Water_RI_energy = 0;  G4double* Water_RI_idx = 0;
    if (mat_Water) {
        G4MaterialPropertiesTable* tbl_Water = mat_Water->GetMaterialPropertiesTable();
        // ABSLENGTH
        get_matprop(tbl_Water, "ABSLENGTH", Water_ABS_n, Water_ABS_energy, Water_ABS_len);
        // RINDEX
        get_matprop(tbl_Water, "RINDEX", Water_RI_n, Water_RI_energy, Water_RI_idx);
    }

    m_op_tree->Branch("Water_ABS_n", &Water_ABS_n, "Water_ABS_n/I");
    m_op_tree->Branch("Water_ABS_energy", Water_ABS_energy, "Water_ABS_energy[Water_ABS_n]/D");
    m_op_tree->Branch("Water_ABS_len", Water_ABS_len, "Water_ABS_len[Water_ABS_n]/D");

    m_op_tree->Branch("Water_RI_n", &Water_RI_n, "Water_RI_n/I");
    m_op_tree->Branch("Water_RI_energy", Water_RI_energy, "Water_RI_energy[Water_RI_n]/D");
    m_op_tree->Branch("Water_RI_idx", Water_RI_idx, "Water_RI_idx[Water_RI_n]/D");

    // ======================================================================
    // * PMT QE (photocathode)
    // * refractive index of PMT glass (Pyrex)
    // ======================================================================
    G4int PC_EFF_n = 0;    G4double* PC_EFF_energy = 0;    G4double* PC_EFF_eff = 0;
    G4int Pyrex_RI_n = 0;  G4double* Pyrex_RI_energy = 0;  G4double* Pyrex_RI_idx = 0;
    G4Material* mat_Photocathode = G4Material::GetMaterial("photocathode");
    if (mat_Photocathode) {
        G4MaterialPropertiesTable* tbl_pc = mat_Photocathode->GetMaterialPropertiesTable();
        // EFFICIENCY
        get_matprop(tbl_pc, "EFFICIENCY", PC_EFF_n, PC_EFF_energy, PC_EFF_eff);
    }
    G4Material* mat_Pyrex = G4Material::GetMaterial("Pyrex");
    if (mat_Pyrex) {
        G4MaterialPropertiesTable* tbl_pyrex = mat_Pyrex->GetMaterialPropertiesTable();
        // RINDEX
        get_matprop(tbl_pyrex, "RINDEX", Pyrex_RI_n, Pyrex_RI_energy, Pyrex_RI_idx);
    }

    m_op_tree->Branch("PC_EFF_n", &PC_EFF_n, "PC_EFF_n/I");
    m_op_tree->Branch("PC_EFF_energy", PC_EFF_energy, "PC_EFF_energy[PC_EFF_n]/D");
    m_op_tree->Branch("PC_EFF_eff", PC_EFF_eff, "PC_EFF_eff[PC_EFF_n]/D");

    m_op_tree->Branch("Pyrex_RI_n", &Pyrex_RI_n, "Pyrex_RI_n/I");
    m_op_tree->Branch("Pyrex_RI_energy", Pyrex_RI_energy, "Pyrex_RI_energy[Pyrex_RI_n]/D");
    m_op_tree->Branch("Pyrex_RI_idx", Pyrex_RI_idx, "Pyrex_RI_idx[Pyrex_RI_n]/D");

    // Finale
    m_op_tree->Fill();
}

void
OpticalParameterAnaMgr::EndOfRunAction(const G4Run* /*run*/) {

}

void
OpticalParameterAnaMgr::BeginOfEventAction(const G4Event* /*event*/) {

}

void
OpticalParameterAnaMgr::EndOfEventAction(const G4Event* /*event*/) {

}

void
OpticalParameterAnaMgr::get_matprop(const G4MaterialPropertyVector* propvec,
        int& N, double*& x, double*& y) {
    if (!propvec) {
        N = 0;
        x = 0;
        y = 0;
        return;
    }
    // int entries = propvec->Entries();
    // N = entries;
    // x = new double[entries];
    // y = new double[entries];
    // for (int i = 0; i < entries; ++i) {
    //     G4MPVEntry entry = propvec->GetEntry(i);
    //     x[i] = entry.GetPhotonEnergy();
    //     y[i] = entry.GetProperty();
    // }

}

void
OpticalParameterAnaMgr::get_matprop(G4MaterialPropertiesTable* table, 
        const char* label, int& N, double*& x, double*& y) {
    G4MaterialPropertyVector* propvec = table->GetProperty(label);
    return get_matprop(propvec, N, x, y);
}
