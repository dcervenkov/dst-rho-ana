#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TH3F.h"
#include "TROOT.h"
#include "TStyle.h"

double interpolate(const TH3* histo, int x_org, int y_org, int z_org, int size) {
    double new_value = 0;
    int points = 0;
    for (int x = x_org - size; x <= x_org + size; x++) {
        for (int y = y_org - size; y <= y_org + size; y++) {
            for (int z = z_org - size; z <= z_org + size; z++) {
                if (x == x_org && y == y_org && z == z_org) continue;
                int bin = histo->GetBin(x, y, z);
                new_value += histo->GetBinContent(bin);
                // printf("delta = %f\n", histo->GetBinContent(bin));
                points++;
            }
        }
    }
    if (points == 0) return 0;
    return new_value/points;
}

int check_efficiency(const char* evtgen_filename, const char* gsim_filename, const char* output_filename) {
    gStyle->SetPalette(kViridis);
    gStyle->SetOptStat(0);

    TFile* gsim_file = new TFile(gsim_filename, "read");
    TH3F* gsim_histo = (TH3F*)((TH3F*)gsim_file->Get("R__TF3"))->Clone("gsim_histo");
    gsim_histo->Rebin3D();
    TFile* evtgen_file = new TFile(evtgen_filename, "read");
    TH3F* evtgen_histo = (TH3F*)((TH3F*)evtgen_file->Get("R__TF3"))->Clone("evtgen_histo");
    evtgen_histo->Rebin3D();

    TH3* eff_histo = dynamic_cast<TH3*>(gsim_histo->Clone("eff_histo"));
    eff_histo->Divide(eff_histo, evtgen_histo, 1.0, 1.0, "B");

    // TH3F* eff_histo = new TH3F("eff_histo", "eff_histo", 100, -1, 1, 100, -1, 1, 100, -1, 1);

    printf("bins: x %i, y %i, z %i\n", gsim_histo->GetXaxis()->GetNbins(), gsim_histo->GetYaxis()->GetNbins(), gsim_histo->GetZaxis()->GetNbins());

    // for (int i = 10000; i < 10010; i++) {
    //     double gsim = gsim_histo->GetBinContent(i);
    //     double evtgen = evtgen_histo->GetBinContent(i);
    //     double eff = eff_histo->GetBinContent(i);
    //     printf("gsim = %f, evtgen = %f, eff = %f\n", gsim, evtgen, eff);
    // }

    int total_bins = 0;
    int fixed_bins = 0;
    for (int x = 0; x <= eff_histo->GetNbinsX() + 1; x++) {
        for (int y = 0; y <= eff_histo->GetNbinsY() + 1; y++) {
            for (int z = 0; z <= eff_histo->GetNbinsZ() + 1; z++) {
                total_bins++;
                int bin = eff_histo->GetBin(x, y, z);
                if (eff_histo->GetBinContent(bin) > 1 || eff_histo->GetBinContent(bin) < 0) {
                    fixed_bins++;
                    int size = 1;
                    double interpolation = -1;
                    do {
                        interpolation = interpolate(eff_histo, x, y, z, size++);
                    } while (interpolation < 0 || interpolation > 1);
                    // printf("bin %i: value = %f, interpolation = %f\n", bin, eff_histo->GetBinContent(bin), interpolation);
                    eff_histo->SetBinContent(bin, interpolation);
                }
            }
        }
    }

    printf("Fixed %i/%i (%.2f%%) bins.\n", fixed_bins, total_bins, (double)fixed_bins/total_bins * 100);

    TFile output_file(output_filename, "RECREATE");
    eff_histo->Write();
    output_file.Close();

    TH1* gsim_1D[3];
    TH1* evtgen_1D[3];
    TH1* eff_1D[3];
    Option_t* proj[] = {"x", "y", "z"};

    TCanvas* c1 = new TCanvas("c1", "c1", 1200, 900);
    c1->Divide(4, 3);

    TH2* eff_2D[3];
    eff_2D[0] = (TH2*)eff_histo->Project3D("xy");
    eff_2D[1] = (TH2*)eff_histo->Project3D("xz");
    eff_2D[2] = (TH2*)eff_histo->Project3D("yz");

    int nbins = eff_histo->GetNbinsX();

    for (int i = 0; i < 3; i++) {
        gsim_1D[i] = gsim_histo->Project3D(proj[i]);
        evtgen_1D[i] = evtgen_histo->Project3D(proj[i]);
        eff_1D[i] = eff_histo->Project3D(proj[i]);

        eff_1D[i]->Scale(1. / nbins / nbins);

        c1->cd(1 + i * 4);
        gsim_1D[i]->Draw();
        c1->cd(2 + i * 4);
        evtgen_1D[i]->Draw();
        c1->cd(3 + i * 4);
        eff_1D[i]->SetAxisRange(0,0.3,"y");
        eff_1D[i]->Draw();

        c1->cd(4 + i * 4);
        eff_2D[i]->Scale(1. / nbins);
        eff_2D[i]->Draw("colz");
        eff_2D[i]->SetAxisRange(0,0.3,"z");
    }

    c1->Modified();
    c1->Update();

    return 0;
}
