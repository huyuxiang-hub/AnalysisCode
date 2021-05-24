#ifndef MuonFastSimVoxelHelper_hh
#define MuonFastSimVoxelHelper_hh

#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TVector3.h>
#include <TProfile2D.h>
#include <TH1I.h>
#include <TMath.h>
#include <TRandom.h>

#include <iostream>
#include <vector>

namespace VoxelMethodHelper {
struct Geom {

    Geom() {
        PMT_x.clear();
        PMT_y.clear();
        PMT_z.clear();

        N = 0;
    }

    void load() {
        const char* preDir = gDirectory->GetPath();
        TFile* f = TFile::Open(m_geom_file.c_str());
        gDirectory->cd(preDir);
        if (!f) { return; }

        TTree* t = dynamic_cast<TTree*>(f->Get("geom"));
        if (!t) { return; }

        int pmttype;
        float x,y,z;
        t->SetBranchAddress("pmttype", &pmttype);
        t->SetBranchAddress("x", &x);
        t->SetBranchAddress("y", &y);
        t->SetBranchAddress("z", &z);

        for (int i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            // only select 20inch PMT
            if (pmttype != 20) {
                continue;
            }
            PMT_type.push_back(pmttype);
            PMT_x.push_back(x);
            PMT_y.push_back(y);
            PMT_z.push_back(z);

            PMT_pos.push_back(TVector3(x,y,z));
        }

        N = PMT_pos.size();
    }

    void info() {
        std::cout << "summary of geom service: " << std::endl;
        std::cout << "size of x: " << PMT_x.size() << std::endl;
        std::cout << "size of y: " << PMT_y.size() << std::endl;
        std::cout << "size of z: " << PMT_z.size() << std::endl;
    }

    std::string m_geom_file;
    std::vector<int> PMT_type;
    std::vector<double> PMT_x;
    std::vector<double> PMT_y;
    std::vector<double> PMT_z;
    std::vector<TVector3> PMT_pos;

    int N;
};

// * 2015.11.26 update the dim of hit time.
// * 2015.11.30 use flags to switch r3 or r, cos(theta) or theta
// The original HitTimeLoader and NPELoader are replaced by V2 since r4309.
// struct HitTimeLoader {};

// * 2015.11.25 update the dim of ths_npe to 100*180
// struct NPELoader {};

/*
 * Optimize the code by creating a base class DataLoader and two derived class
 */

struct BaseDataLoader {
    enum Mode {
        kUnknown = 0,
        kR3 = 1, kR = 2,
        kCosTheta = 3, kTheta = 4
    };

    BaseDataLoader(Mode _xmode, Int_t _xbin,
                   Mode _ymode, Int_t _ybin,
                   bool _lazy) {
        // initialized with default mode and bin
        xmode = _xmode;
        xbinnum = _xbin;
        ymode = _ymode;
        ybinnum = _ybin;

        m_lazy_loading = _lazy;

        bins_r3 = 0;
        bins_r = 0;
        bins_costheta = 0;
        bins_theta = 0;
        
        f_main = 0;
    }

    void load_from(TString fn) {
        if (f_main) {
            return;
        }

        const char* preDir = gDirectory->GetPath();
        f_main = TFile::Open(fn);
        gDirectory->cd(preDir);
        if (!f_main) { 
            return;
        }

        // ##############################################################################
        // # Load the aux data
        // ##############################################################################

        // X

        if (f_main->Get("bins_r3")) {
            bins_r3 = dynamic_cast<TAxis*>(f_main->Get("bins_r3"));
            xmode = kR3;
            xbinnum = bins_r3->GetNbins();
            std::cout << "===== Load bins_r3 (" << xbinnum << ") " << std::endl;
        } else if (f_main->Get("bins_r")) {
            bins_r = dynamic_cast<TAxis*>(f_main->Get("bins_r"));
            xmode = kR;
            xbinnum = bins_r->GetNbins();
            std::cout << "===== Load bins_r (" << xbinnum << ") " << std::endl;
        } else {
            bins_r3 = new TAxis(xbinnum, 0, 5600);
            bins_r = new TAxis(xbinnum, 0, 17.7);
        }

        // Y

        if (f_main->Get("bins_theta")) {
            bins_theta = dynamic_cast<TAxis*>(f_main->Get("bins_theta"));
            ymode = kTheta;
            ybinnum = bins_theta->GetNbins();
            std::cout << "===== Load bins_theta (" << ybinnum << ") " << std::endl;
        } else if (f_main->Get("bins_costheta")) {
            bins_costheta = dynamic_cast<TAxis*>(f_main->Get("bins_costheta"));
            ymode = kCosTheta;
            ybinnum = bins_costheta->GetNbins();
            std::cout << "===== Load bins_costheta (" << ybinnum << ") " << std::endl;
        } else {
            bins_costheta = new TAxis(ybinnum, -1, 1.01);
            bins_theta = new TAxis(ybinnum, 0, 180.01*TMath::Pi()/180.);
        }

        // Load histograms
        for (int i = 0; i < xbinnum; ++i) {
            for (int j = 0; j < ybinnum; ++j) {
                TH1* h = 0;
                if (!m_lazy_loading) {
                    TString th_name = TString::Format("%d", i*ybinnum+j);
                    h = dynamic_cast<TH1*>(f_main->Get(th_name));
                    assert(h);
                }
                ths_main.push_back(h);
            }
        }
        std::cout << "==== Total load " << ths_main.size() << std::endl;
    }

    virtual Int_t get_bin_x(Float_t r) { // unit is m
        Int_t binx = 1;
        if (xmode == kR3) {
            binx = bins_r3->FindBin(TMath::Power(r, 3));
        } else if (xmode == kR) {
            binx = bins_r->FindBin(r);
        } else {
            std::cerr << "unknown mode" << std::endl;
        }
        return binx;
    }

    virtual Int_t get_bin_y(Float_t theta) { // unit is rad
        Int_t biny = 1;
        if (ymode == kCosTheta) {
            biny = bins_costheta->FindBin(TMath::Cos(theta));
        } else if (ymode == kTheta) {
            biny = bins_theta->FindBin(theta);
        } else {
            std::cerr << "unknown mode" << std::endl;
        }
        return biny;
    }

    virtual TH1* get_hist(Int_t binx, Int_t biny) {
        // binx: [1, xbinnum]
        // biny: [1, ybinnum]
        if (binx<1) { binx = 1; }
        else if (binx > xbinnum) { binx = xbinnum;}
        if (biny<1) { biny = 1; }
        else if (biny > ybinnum) { biny = ybinnum;}

        int idx = (binx-1)*ybinnum+(biny-1);

        TH1* h = ths_main[idx];
        if (!h) {
            TString th_name = TString::Format("%d", idx);
            h = dynamic_cast<TH1*>(f_main->Get(th_name));
            ths_main[idx] = h;
        }
        return h;
    }


    Mode xmode; // kR3 or kR
    Mode ymode; // kCosTheta or kTheta

    Int_t xbinnum;
    Int_t ybinnum;

    TAxis* bins_r3;
    TAxis* bins_r;
    TAxis* bins_costheta;
    TAxis* bins_theta;

    bool m_lazy_loading;

    TFile* f_main;

    std::vector<TH1*> ths_main; // the histograms

    
};

struct HitTimeLoaderV2: BaseDataLoader {
    HitTimeLoaderV2()
        : BaseDataLoader(kR, 200, kTheta, 180, false) {
    }

    void load() {
        load_from(m_single_filename);
    }

    Double_t get_hittime(Float_t r, Float_t theta) {
        Int_t binx = get_bin_x(r);
        Int_t biny = get_bin_y(theta);
        return get_hittime(binx, biny);
    }

    Double_t get_hittime(Int_t binx, Int_t biny, int mode) {

        TH1* h = get_hist(binx, biny);
        Double_t hittime = h->GetRandom();

        return hittime;
    }

    TString m_single_filename;
};

struct NPELoaderV2: BaseDataLoader {
    NPELoaderV2()
        : BaseDataLoader(kR, 200, kTheta, 180, false) {
    }

    void load() {
        load_from(m_filename_npe_single);
    }

    Int_t get_npe(Float_t r, Float theta) {
        Int_t binx = get_bin_x(r);
        Int_t biny = get_bin_y(theta);
        return get_npe(binx, biny);
    }

    Int_t get_npe(Int_t binx, Int_t biny) {
        Int_t npe_from_single = 0;
        if (1<=binx and binx<=xbinnum and 1<=biny and biny<=ybinnum) {
            TH1* th = get_hist(binx, biny);
            npe_from_single = th->GetRandom();
        } else if (binx==1 and (biny<1 or biny>ybinnum)) {
            biny = gRandom->Uniform(1,ybinnum);
            TH1* th = get_hist(binx, biny);
            npe_from_single = th->GetRandom();
        } else if (binx>1 and (biny<1 or biny>ybinnum)) {
            // std::cerr << "npe maybe lost: " << binx << "/" << biny << std::endl;
            // FIXME how to handle such situation.
            // biny = gRandom->Uniform(1,100);
            if (biny>ybinnum) { biny = ybinnum; }
            else if (biny<1){ biny = 1; }

            TH1* th = get_hist(binx, biny);
            npe_from_single = th->GetRandom();
        } else {
            static long warning = 0;
            ++warning;
            if (warning < 10) {
                std::cerr << "npe lost: " << binx << "/" << biny << std::endl;
            } else if (warning == 10) {
                std::cerr << "too many npe lost complains." << std::endl;
            }
        }

        return npe_from_single;
    }

    TString m_filename_npe_single;
};

}

#endif
