/**
 *  @file    fittershe.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2018-06-05
 *
 *  @brief A class to calculate a spherical harmonics expansion of an
 *  arbitrary 3D function.
 * 
 * This class can calculate a spherical harmonics expansion (SHE) of an
 * arbitrary 3D function with the appropriate domain. This can be used to model
 * efficiency/acceptance in a particle physics analysis. Most of this code was
 * donated by Pavel Reznicek.
 *
 */

#include "fittershe.h"

// Standard includes
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

// ROOT includes
#include "Math/SpecFunc.h"
#include "Rtypes.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TF3.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3.h"
#include "TH3D.h"
#include "TMath.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"

// Local includes
#include "cksum.h"
#include "colors.h"
#include "constants.h"
#include "gitversion.h"
#include "log.h"
#include "tools.h"

// int main(int argc, char const *argv[]) {
//     colors::setColors();
//     // tools::SetupPlotStyle();
//     gStyle->SetOptStat(0);
//     const char* dataname = argv[1];
//     printf("Processing %s...\n", dataname);
//     AnalyzeDataset(dataname, -1, 20, 50, 2, 1, 5, 5, "fast");
//     return 0;
// }
FitterSHE::FitterSHE() {
}

FitterSHE::~FitterSHE() {
    if (output_file_) {
        if (output_file_->IsOpen()) {
            output_file_->Close();
        }
    }
}

/**
 * Reads in data from ROOT file(s). Constructs separate datasets for the 4 categories.
 * Binds the variables to the dataset, so that dataset->get(i) changes values of, e.g., expno_
 *
 * @param file_names Vector of paths to the ROOT files
 * @param num_events [optional] Maximum number of events to use (0 to read all)
 */
void FitterSHE::ReadInFile(std::vector<const char*> file_names, const int& num_events) {
    TChain* input_chain = new TChain("h2000");
    for (auto file_name : file_names) {
        input_chain->Add(file_name);
    }

    TTree* input_tree;
    if (num_events) {
        input_tree = input_chain->CloneTree(num_events);
    } else {
        input_tree = input_chain->CloneTree();
    }

    TString cuts = tools::GetCommonCutsString();
    cuts += "&&evmcflag!=1";

    TTree* temp_tree = input_tree->CopyTree(cuts.Data());
    delete input_tree;
    input_data_ = temp_tree;
}

/**
 * Set file in which output will be saved.
 */
void FitterSHE::SetOutputFile(const char* filename) {
    output_file_ = new TFile(filename, "RECREATE");
}

Double_t FitterSHE::Ylm(UInt_t l, Int_t m, Double_t thetat, Double_t phi) {
    // https://en.wikipedia.org/wiki/Spherical_harmonics (real form)
    // https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
    if (m > 0)
        return TMath::Power(-1, m) * TMath::Sqrt2() * ROOT::Math::sph_legendre(l, m, thetat) *
               TMath::Cos(m * phi);  // integral is 1
    else if (m < 0)
        return TMath::Power(-1, m) * TMath::Sqrt2() * ROOT::Math::sph_legendre(l, -m, thetat) *
               TMath::Sin(-m * phi);  // integral is 1
    else
        return ROOT::Math::sph_legendre(l, 0, thetat);  // integral is 1
}

Double_t FitterSHE::Pk(UInt_t k, Double_t thetab) {
    // https://en.wikipedia.org/wiki/Legendre_polynomials
    Double_t cosThetab = TMath::Cos(thetab);
    return ROOT::Math::legendre(k, cosThetab) * TMath::Sqrt(k + 0.5);  // integral is 1
}

Double_t FitterSHE::SumAklmYlmPk(double *x, double *par) {
    Double_t phi = x[0];
    Double_t thetat = x[1];
    Double_t thetab = x[2];

    // while ( phi < -TMath::Pi() ) phi += 2*TMath::Pi();
    // while ( phi > +TMath::Pi() ) phi -= 2*TMath::Pi();

    Double_t &norm = par[0];
    Int_t max_k = (Int_t)(par[1] + 0.1);
    Int_t max_l = (Int_t)(par[2] + 0.1);

    Double_t result = 0;

    Int_t i = 3;  // par[0],[1],[2] are constants: normalization, max_k and max_l
    for (Int_t k = 0; k <= max_k; k++) {
        for (Int_t l = 0; l <= max_l; l++) {
            for (Int_t m = -l; m <= l; m++) {
                Double_t &aklm = par[i];
                result += aklm * Ylm(l, m, thetat, phi) * Pk(k, thetab);
                i += 1;
            }
        }
    }

    return norm * result;
}

std::map<TString, Double_t> FitterSHE::GetChi2NDF(TH3D *h3_phi_cosTheta1_cosTheta2, TF3 *f3, TString title,
                                       TString options, Int_t smoothScale) {
    static Int_t counter = 0;

    Bool_t opt_useIntegral = options.Contains("integ", TString::kIgnoreCase);
    Bool_t opt_controlPlots = options.Contains("plot", TString::kIgnoreCase);
    Bool_t opt_textOutput = options.Contains("text", TString::kIgnoreCase);

    Int_t nbins_phi = h3_phi_cosTheta1_cosTheta2->GetNbinsX();
    Int_t nbins_cosTheta1 = h3_phi_cosTheta1_cosTheta2->GetNbinsY();
    Int_t nbins_cosTheta2 = h3_phi_cosTheta1_cosTheta2->GetNbinsZ();

    // Create 3D histogram from the TF3 with binning "smoothScale"-times finer than the data
    // histogram
    f3->SetNpx(h3_phi_cosTheta1_cosTheta2->GetNbinsX() * smoothScale);
    f3->SetNpy(h3_phi_cosTheta1_cosTheta2->GetNbinsY() * smoothScale);
    f3->SetNpz(h3_phi_cosTheta1_cosTheta2->GetNbinsZ() * smoothScale);
    TH3 *hf3_phi_cosTheta1_cosTheta2 = (TH3 *)(f3->GetHistogram());
    Double_t hf3_nBinsX = hf3_phi_cosTheta1_cosTheta2->GetNbinsX();
    Double_t hf3_nBinsY = hf3_phi_cosTheta1_cosTheta2->GetNbinsY();
    Int_t hf3_nBinsZ = hf3_phi_cosTheta1_cosTheta2->GetNbinsZ();
    for (Int_t ix = 1; ix <= hf3_phi_cosTheta1_cosTheta2->GetNbinsX(); ix++) {
        for (Int_t iy = 1; iy <= hf3_phi_cosTheta1_cosTheta2->GetNbinsY(); iy++) {
            if ((iy % ((100 / hf3_nBinsZ) + 1)) == 0) {
                printf("Evaluating function: %.5g%%          \r",
                       100. * ((ix - 1) * hf3_nBinsY + (iy - 1)) / (hf3_nBinsX * hf3_nBinsY));
                std::flush(std::cout);
            }
            for (Int_t iz = 1; iz <= hf3_phi_cosTheta1_cosTheta2->GetNbinsZ(); iz++) {
                double xpos = hf3_phi_cosTheta1_cosTheta2->GetXaxis()->GetBinCenter(ix);
                double ypos = hf3_phi_cosTheta1_cosTheta2->GetYaxis()->GetBinCenter(iy);
                double zpos = hf3_phi_cosTheta1_cosTheta2->GetZaxis()->GetBinCenter(iz);
                if (opt_useIntegral)
                    hf3_phi_cosTheta1_cosTheta2->SetBinContent(
                        ix, iy, iz,
                        f3->Eval(xpos, ypos, zpos));  // TODO: replace eval by integral in the bin
                else
                    hf3_phi_cosTheta1_cosTheta2->SetBinContent(ix, iy, iz,
                                                               f3->Eval(xpos, ypos, zpos));
            }
        }
    }
    printf("                                     \r");
    std::flush(std::cout);

    TFile output_file(title + ".root", "RECREATE");
    hf3_phi_cosTheta1_cosTheta2->Write();
    output_file.Close();

    // Create 1D and 2D projections from the data histogram
    TH1D *h1_phi = (TH1D *)h3_phi_cosTheta1_cosTheta2->Project3D("x");
    TH1D *h1_cosTheta1 = (TH1D *)h3_phi_cosTheta1_cosTheta2->Project3D("y");
    TH1D *h1_cosTheta2 = (TH1D *)h3_phi_cosTheta1_cosTheta2->Project3D("z");
    TH2D *h2_phi_cosTheta1 = (TH2D *)h3_phi_cosTheta1_cosTheta2->Project3D("yx");
    TH2D *h2_phi_cosTheta2 = (TH2D *)h3_phi_cosTheta1_cosTheta2->Project3D("zx");
    TH2D *h2_cosTheta1_cosTheta2 = (TH2D *)h3_phi_cosTheta1_cosTheta2->Project3D("zy");

    // Create 1D and 2D projections from the TF3 histogram
    TH1D *hf1_phi = (TH1D *)hf3_phi_cosTheta1_cosTheta2->Project3D("x");
    TH1D *hf1_cosTheta1 = (TH1D *)hf3_phi_cosTheta1_cosTheta2->Project3D("y");
    TH1D *hf1_cosTheta2 = (TH1D *)hf3_phi_cosTheta1_cosTheta2->Project3D("z");
    TH2D *hf2_phi_cosTheta1 = (TH2D *)hf3_phi_cosTheta1_cosTheta2->Project3D("yx");
    TH2D *hf2_phi_cosTheta2 = (TH2D *)hf3_phi_cosTheta1_cosTheta2->Project3D("zx");
    TH2D *hf2_cosTheta1_cosTheta2 = (TH2D *)hf3_phi_cosTheta1_cosTheta2->Project3D("zy");
    hf1_phi->Scale(1. / smoothScale / smoothScale);
    hf1_cosTheta1->Scale(1. / smoothScale / smoothScale);
    hf1_cosTheta2->Scale(1. / smoothScale / smoothScale);
    hf2_phi_cosTheta1->Scale(1. / smoothScale);
    hf2_phi_cosTheta2->Scale(1. / smoothScale);
    hf2_cosTheta1_cosTheta2->Scale(1. / smoothScale);

    // Calculate chi2/NDF for all the projections and fill-up the histograms of residuals
    TH1D hres_phi("hres_phi", "hres_phi", 200, -10, 10);
    TH1D hres_cosTheta1("hres_cosTheta1", "hres_cosTheta1", 200, -10, 10);
    TH1D hres_cosTheta2("hres_cosTheta2", "hres_cosTheta2", 200, -10, 10);
    TH1D hres_phi_cosTheta1("hres_phi_cosTheta1", "hres_phi_cosTheta1", 200, -10, 10);
    TH1D hres_phi_cosTheta2("hres_phi_cosTheta2", "hres_phi_cosTheta2", 200, -10, 10);
    TH1D hres_cosTheta1_cosTheta2("hres_cosTheta1_cosTheta2", "hres_cosTheta1_cosTheta2", 200, -10,
                                  10);
    TH1D hres_phi_cosTheta1_cosTheta2("hres_phi_cosTheta1_cosTheta2",
                                      "hres_phi_cosTheta1_cosTheta2", 200, -10, 10);
    Int_t NDF_phi = nbins_phi /*- f3->GetNumberFreeParameters()*/;
    Int_t NDF_cosTheta1 = nbins_cosTheta1 /*- f3->GetNumberFreeParameters()*/;
    Int_t NDF_cosTheta2 = nbins_cosTheta2 /*- f3->GetNumberFreeParameters()*/;
    Int_t NDF_phi_cosTheta1 = nbins_phi * nbins_cosTheta1 /*- f3->GetNumberFreeParameters()*/;
    Int_t NDF_phi_cosTheta2 = nbins_phi * nbins_cosTheta2 /*- f3->GetNumberFreeParameters()*/;
    Int_t NDF_cosTheta1_cosTheta2 =
        nbins_cosTheta1 * nbins_cosTheta2 /*- f3->GetNumberFreeParameters()*/;
    Int_t NDF_phi_cosTheta1_cosTheta2 =
        nbins_phi * nbins_cosTheta1 * nbins_cosTheta2 /*- f3->GetNumberFreeParameters()*/;
    Double_t chi2_phi = 0;
    Double_t chi2_cosTheta1 = 0;
    Double_t chi2_cosTheta2 = 0;
    Double_t chi2_phi_cosTheta1 = 0;
    Double_t chi2_phi_cosTheta2 = 0;
    Double_t chi2_cosTheta1_cosTheta2 = 0;
    Double_t chi2_phi_cosTheta1_cosTheta2 = 0;
    Double_t res = 0;
    // 1D
    for (Int_t ix = 1; ix <= nbins_phi; ix++) {
        res = 0;
        for (Int_t ifx = 1 + (ix - 1) * smoothScale; ifx < 1 + ix * smoothScale; ifx++) {
            res -= hf1_phi->GetBinContent(ifx);
        }
        res /= smoothScale;
        res += h1_phi->GetBinContent(ix);
        if (h1_phi->GetBinError(ix) > 0) {
            chi2_phi += res * res / TMath::Sq(h1_phi->GetBinError(ix));
            hres_phi.Fill(res / h1_phi->GetBinError(ix));
        }
    }
    for (Int_t ix = 1; ix <= nbins_cosTheta1; ix++) {
        res = 0;
        for (Int_t ifx = 1 + (ix - 1) * smoothScale; ifx < 1 + ix * smoothScale; ifx++) {
            res -= hf1_cosTheta1->GetBinContent(ifx);
        }
        res /= smoothScale;
        res += h1_cosTheta1->GetBinContent(ix);
        if (h1_cosTheta1->GetBinError(ix) > 0) {
            chi2_cosTheta1 += res * res / TMath::Sq(h1_cosTheta1->GetBinError(ix));
            hres_cosTheta1.Fill(res / h1_cosTheta1->GetBinError(ix));
        }
    }
    for (Int_t ix = 1; ix <= nbins_cosTheta2; ix++) {
        res = 0;
        for (Int_t ifx = 1 + (ix - 1) * smoothScale; ifx < 1 + ix * smoothScale; ifx++) {
            res -= hf1_cosTheta2->GetBinContent(ifx);
        }
        res /= smoothScale;
        res += h1_cosTheta2->GetBinContent(ix);
        if (h1_cosTheta2->GetBinError(ix) > 0) {
            chi2_cosTheta2 += res * res / TMath::Sq(h1_cosTheta2->GetBinError(ix));
            hres_cosTheta2.Fill(res / h1_cosTheta2->GetBinError(ix));
        }
    }
    // 2D
    for (Int_t ix = 1; ix <= nbins_phi; ix++) {
        for (Int_t iy = 1; iy <= nbins_cosTheta1; iy++) {
            res = 0;
            for (Int_t ifx = 1 + (ix - 1) * smoothScale; ifx < 1 + ix * smoothScale; ifx++) {
                for (Int_t ify = 1 + (iy - 1) * smoothScale; ify < 1 + iy * smoothScale; ify++) {
                    res -= hf2_phi_cosTheta1->GetBinContent(ifx, ify);
                }
            }
            res /= smoothScale * smoothScale;
            res += h2_phi_cosTheta1->GetBinContent(ix, iy);
            if (h2_phi_cosTheta1->GetBinError(ix, iy) > 0) {
                chi2_phi_cosTheta1 += res * res / TMath::Sq(h2_phi_cosTheta1->GetBinError(ix, iy));
                hres_phi_cosTheta1.Fill(res / h2_phi_cosTheta1->GetBinError(ix, iy));
            }
        }
    }
    for (Int_t ix = 1; ix <= nbins_phi; ix++) {
        for (Int_t iy = 1; iy <= nbins_cosTheta2; iy++) {
            res = 0;
            for (Int_t ifx = 1 + (ix - 1) * smoothScale; ifx < 1 + ix * smoothScale; ifx++) {
                for (Int_t ify = 1 + (iy - 1) * smoothScale; ify < 1 + iy * smoothScale; ify++) {
                    res -= hf2_phi_cosTheta2->GetBinContent(ifx, ify);
                }
            }
            res /= smoothScale * smoothScale;
            res += h2_phi_cosTheta2->GetBinContent(ix, iy);
            if (h2_phi_cosTheta2->GetBinError(ix, iy) > 0) {
                chi2_phi_cosTheta2 += res * res / TMath::Sq(h2_phi_cosTheta2->GetBinError(ix, iy));
                hres_phi_cosTheta2.Fill(res / h2_phi_cosTheta2->GetBinError(ix, iy));
            }
        }
    }
    for (Int_t ix = 1; ix <= nbins_cosTheta1; ix++) {
        for (Int_t iy = 1; iy <= nbins_cosTheta2; iy++) {
            res = 0;
            for (Int_t ifx = 1 + (ix - 1) * smoothScale; ifx < 1 + ix * smoothScale; ifx++) {
                for (Int_t ify = 1 + (iy - 1) * smoothScale; ify < 1 + iy * smoothScale; ify++) {
                    res -= hf2_cosTheta1_cosTheta2->GetBinContent(ifx, ify);
                }
            }
            res /= smoothScale * smoothScale;
            res += h2_cosTheta1_cosTheta2->GetBinContent(ix, iy);
            if (h2_cosTheta1_cosTheta2->GetBinError(ix, iy) > 0) {
                chi2_cosTheta1_cosTheta2 +=
                    res * res / TMath::Sq(h2_cosTheta1_cosTheta2->GetBinError(ix, iy));
                hres_cosTheta1_cosTheta2.Fill(res / h2_cosTheta1_cosTheta2->GetBinError(ix, iy));
            }
        }
    }
    // 3D
    for (Int_t ix = 1; ix <= nbins_phi; ix++) {
        for (Int_t iy = 1; iy <= nbins_cosTheta1; iy++) {
            for (Int_t iz = 1; iz <= nbins_cosTheta2; iz++) {
                res = 0;
                for (Int_t ifx = 1 + (ix - 1) * smoothScale; ifx < 1 + ix * smoothScale; ifx++) {
                    for (Int_t ify = 1 + (iy - 1) * smoothScale; ify < 1 + iy * smoothScale;
                         ify++) {
                        for (Int_t ifz = 1 + (iz - 1) * smoothScale; ifz < 1 + iz * smoothScale;
                             ifz++) {
                            res -= hf3_phi_cosTheta1_cosTheta2->GetBinContent(ifx, ify, ifz);
                        }
                    }
                }
                res /= smoothScale * smoothScale * smoothScale;
                res += h3_phi_cosTheta1_cosTheta2->GetBinContent(ix, iy, iz);
                if (h3_phi_cosTheta1_cosTheta2->GetBinError(ix, iy, iz) > 0) {
                    chi2_phi_cosTheta1_cosTheta2 +=
                        res * res / TMath::Sq(h3_phi_cosTheta1_cosTheta2->GetBinError(ix, iy, iz));
                    hres_phi_cosTheta1_cosTheta2.Fill(
                        res / h3_phi_cosTheta1_cosTheta2->GetBinError(ix, iy, iz));
                }
            }
        }
    }

    // Print the chi2/NDF
    TString text_phi =
        TString::Format("chi2/NDF = %.3g = %.5g / %d", chi2_phi / NDF_phi, chi2_phi, NDF_phi);
    TString text_cosTheta1 =
        TString::Format("chi2/NDF = %.3g = %.5g / %d", chi2_cosTheta1 / NDF_cosTheta1,
                        chi2_cosTheta1, NDF_cosTheta1);
    TString text_cosTheta2 =
        TString::Format("chi2/NDF = %.3g = %.5g / %d", chi2_cosTheta2 / NDF_cosTheta2,
                        chi2_cosTheta2, NDF_cosTheta2);
    TString text_phi_cosTheta1 =
        TString::Format("chi2/NDF = %.3g = %.5g / %d", chi2_phi_cosTheta1 / NDF_phi_cosTheta1,
                        chi2_phi_cosTheta1, NDF_phi_cosTheta1);
    TString text_phi_cosTheta2 =
        TString::Format("chi2/NDF = %.3g = %.5g / %d", chi2_phi_cosTheta2 / NDF_phi_cosTheta2,
                        chi2_phi_cosTheta2, NDF_phi_cosTheta2);
    TString text_cosTheta1_cosTheta2 = TString::Format(
        "chi2/NDF = %.3g = %.5g / %d", chi2_cosTheta1_cosTheta2 / NDF_cosTheta1_cosTheta2,
        chi2_cosTheta1_cosTheta2, NDF_cosTheta1_cosTheta2);
    TString text_phi_cosTheta1_cosTheta2 = TString::Format(
        "chi2/NDF = %.3g = %.5g / %d", chi2_phi_cosTheta1_cosTheta2 / NDF_phi_cosTheta1_cosTheta2,
        chi2_phi_cosTheta1_cosTheta2, NDF_phi_cosTheta1_cosTheta2);

    if (opt_textOutput) {
        printf("1D phi      :  %s\n", text_phi.Data());
        printf("1D cosTheta1:  %s\n", text_cosTheta1.Data());
        printf("1D cosTheta2:  %s\n", text_cosTheta2.Data());
        printf("2D phi       x cosTheta1:  %s\n", text_phi_cosTheta1.Data());
        printf("2D phi       x cosTheta2:  %s\n", text_phi_cosTheta2.Data());
        printf("2D cosTheta1 x cosTheta2:  %s\n", text_cosTheta1_cosTheta2.Data());
        printf("3D phi x cosTheta1 x cosTheta2:  %s\n", text_phi_cosTheta1_cosTheta2.Data());
    }

    hres_phi.SetXTitle(text_phi);
    hres_cosTheta1.SetXTitle(text_cosTheta1);
    hres_cosTheta2.SetXTitle(text_cosTheta2);
    hres_phi_cosTheta1.SetXTitle(text_phi_cosTheta1);
    hres_phi_cosTheta2.SetXTitle(text_phi_cosTheta2);
    hres_cosTheta1_cosTheta2.SetXTitle(text_cosTheta1_cosTheta2);
    hres_phi_cosTheta1_cosTheta2.SetXTitle(text_phi_cosTheta1_cosTheta2);

    // Draw the projections and residuals
    if (opt_controlPlots) {
        Double_t max_phi = h1_phi->GetMaximum() * 1.2;
        Double_t max_cosTheta1 = h1_cosTheta1->GetMaximum() * 1.2;
        Double_t max_cosTheta2 = h1_cosTheta2->GetMaximum() * 1.2;
        Double_t max_phi_cosTheta1 = h2_phi_cosTheta1->GetMaximum() * 1.2;
        Double_t max_phi_cosTheta2 = h2_phi_cosTheta2->GetMaximum() * 1.2;
        Double_t max_cosTheta1_cosTheta2 = h2_cosTheta1_cosTheta2->GetMaximum() * 1.2;
        TCanvas *c = new TCanvas("c_chi2ndf_" + TString::Format("%d", counter),
                                 "Chi2/NDF and projections: " + title, 7 * 400, 3 * 400);
        c->Draw();
        c->Divide(7, 3);
        c->cd(1);
        h1_phi->DrawCopy("hist")->GetYaxis()->SetRangeUser(0, max_phi);
        c->cd(1 + 7);
        hf1_phi->DrawCopy("hist")->GetYaxis()->SetRangeUser(0, max_phi);
        c->cd(1 + 14);
        hres_phi.DrawCopy("");
        c->cd(2);
        h1_cosTheta1->DrawCopy("hist")->GetYaxis()->SetRangeUser(0, max_cosTheta1);
        c->cd(2 + 7);
        hf1_cosTheta1->DrawCopy("hist")->GetYaxis()->SetRangeUser(0, max_cosTheta1);
        c->cd(2 + 14);
        hres_cosTheta1.DrawCopy("");
        c->cd(3);
        h1_cosTheta2->DrawCopy("hist")->GetYaxis()->SetRangeUser(0, max_cosTheta2);
        c->cd(3 + 7);
        hf1_cosTheta2->DrawCopy("hist")->GetYaxis()->SetRangeUser(0, max_cosTheta2);
        c->cd(3 + 14);
        hres_cosTheta2.DrawCopy("");
        c->cd(4);
        h2_phi_cosTheta1->DrawCopy("colz")->GetZaxis()->SetRangeUser(0, max_phi_cosTheta1);
        c->cd(4 + 7);
        hf2_phi_cosTheta1->DrawCopy("colz")->GetZaxis()->SetRangeUser(0, max_phi_cosTheta1);
        c->cd(4 + 14);
        hres_phi_cosTheta1.DrawCopy("");
        c->cd(5);
        h2_phi_cosTheta2->DrawCopy("colz")->GetZaxis()->SetRangeUser(0, max_phi_cosTheta2);
        c->cd(5 + 7);
        hf2_phi_cosTheta2->DrawCopy("colz")->GetZaxis()->SetRangeUser(0, max_phi_cosTheta2);
        c->cd(5 + 14);
        hres_phi_cosTheta2.DrawCopy("");
        c->cd(6);
        h2_cosTheta1_cosTheta2->DrawCopy("colz")->GetZaxis()->SetRangeUser(0,
                                                                           max_cosTheta1_cosTheta2);
        c->cd(6 + 7);
        hf2_cosTheta1_cosTheta2->DrawCopy("colz")->GetZaxis()->SetRangeUser(
            0, max_cosTheta1_cosTheta2);
        c->cd(6 + 14);
        hres_cosTheta1_cosTheta2.DrawCopy("");
        c->cd(7);
        h3_phi_cosTheta1_cosTheta2->DrawCopy("lego");
        c->cd(7 + 7);
        hf3_phi_cosTheta1_cosTheta2->DrawCopy("lego");
        c->cd(7 + 14);
        hres_phi_cosTheta1_cosTheta2.DrawCopy("");
        c->Modified();
        c->Update();
        c->SaveAs(title + ".png");
    }

    std::map<TString, Double_t> result;
    result["NDF_phi"] = NDF_phi;
    result["chi2_phi"] = chi2_phi;
    result["chi2NDF_phi"] = chi2_phi / NDF_phi;
    result["prob_phi"] = TMath::Prob(chi2_phi, NDF_phi);
    result["NDF_cosTheta1"] = NDF_cosTheta1;
    result["chi2_cosTheta1"] = chi2_cosTheta1;
    result["chi2NDF_cosTheta1"] = chi2_cosTheta1 / NDF_cosTheta1;
    result["prob_cosTheta1"] = TMath::Prob(chi2_cosTheta1, NDF_cosTheta1);
    result["NDF_cosTheta2"] = NDF_cosTheta2;
    result["chi2_cosTheta2"] = chi2_cosTheta2;
    result["chi2NDF_cosTheta2"] = chi2_cosTheta2 / NDF_cosTheta2;
    result["prob_cosTheta2"] = TMath::Prob(chi2_cosTheta2, NDF_cosTheta2);
    result["NDF_phi_cosTheta1"] = NDF_phi_cosTheta1;
    result["chi2_phi_cosTheta1"] = chi2_phi_cosTheta1;
    result["chi2NDF_phi_cosTheta1"] = chi2_phi_cosTheta1 / NDF_phi_cosTheta1;
    result["prob_phi_cosTheta1"] = TMath::Prob(chi2_phi_cosTheta1, NDF_phi_cosTheta1);
    result["NDF_phi_cosTheta2"] = NDF_phi_cosTheta2;
    result["chi2_phi_cosTheta2"] = chi2_phi_cosTheta2;
    result["chi2NDF_phi_cosTheta2"] = chi2_phi_cosTheta2 / NDF_phi_cosTheta2;
    result["prob_phi_cosTheta2"] = TMath::Prob(chi2_phi_cosTheta2, NDF_phi_cosTheta2);
    result["NDF_cosTheta1_cosTheta2"] = NDF_cosTheta1_cosTheta2;
    result["chi2_cosTheta1_cosTheta2"] = chi2_cosTheta1_cosTheta2;
    result["chi2NDF_cosTheta1_cosTheta2"] = chi2_cosTheta1_cosTheta2 / NDF_cosTheta1_cosTheta2;
    result["prob_cosTheta1_cosTheta2"] =
        TMath::Prob(chi2_cosTheta1_cosTheta2, NDF_cosTheta1_cosTheta2);
    result["NDF_phi_cosTheta1_cosTheta2"] = NDF_phi_cosTheta1_cosTheta2;
    result["chi2_phi_cosTheta1_cosTheta2"] = chi2_phi_cosTheta1_cosTheta2;
    result["chi2NDF_phi_cosTheta1_cosTheta2"] =
        chi2_phi_cosTheta1_cosTheta2 / NDF_phi_cosTheta1_cosTheta2;
    result["prob_phi_cosTheta1_cosTheta2"] =
        TMath::Prob(chi2_phi_cosTheta1_cosTheta2, NDF_phi_cosTheta1_cosTheta2);

    counter++;
    return result;
}

void FitterSHE::AnalyzeDataset(TString dataset, Int_t max_entries, Int_t nBins3D, Int_t nBins1D,
                    Int_t smooth3D, Int_t smooth1D, Int_t max_k, Int_t max_l, TString options) {
    Bool_t opt_checkNorm = options.Contains("norm", TString::kIgnoreCase);
    Bool_t opt_fastSearch = options.Contains("fast", TString::kIgnoreCase);
    Bool_t opt_noSearch = options.Contains("fullonly", TString::kIgnoreCase);

    // Read input ntuple, define angles titles for different datasets
    TChain *t_data;
    TString selection;
    TString angle_names[3];
    TString angle_title[3];

    selection = tools::GetCommonCutsString();
    selection += "&&evmcflag!=1";
    t_data = new TChain("h2000");
    t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-charged_s10.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-charged_s11.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-charged_s12.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-charged_s13.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-charged_s14.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-charged_s15.root");
    t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-charm_s10.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-charm_s11.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-charm_s12.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-charm_s13.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-charm_s14.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-charm_s15.root");
    t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-mixed_s10.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-mixed_s11.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-mixed_s12.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-mixed_s13.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-mixed_s14.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-mixed_s15.root");
    t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-uds_s10.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-uds_s11.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-uds_s12.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-uds_s13.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-uds_s14.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd1_on_resonance_evtgen-uds_s15.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-charged_s0.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-charged_s1.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-charged_s2.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-charged_s3.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-charged_s4.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-charged_s5.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-charm_s0.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-charm_s1.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-charm_s2.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-charm_s3.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-charm_s4.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-charm_s5.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-mixed_s0.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-mixed_s1.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-mixed_s2.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-mixed_s3.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-mixed_s4.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-mixed_s5.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-uds_s0.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-uds_s1.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-uds_s2.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-uds_s3.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-uds_s4.root");
    // t_data->Add("../data/Kpi/mc/DSRhoSkim_svd2_on_resonance_evtgen-uds_s5.root");

    angle_names[0] = "phit";
    angle_title[0] = "#phi_{t}";  // phi must be 1st !!!, cos theta1/2 to follow
    angle_names[1] = "thetat";
    angle_title[1] = "#theta_{t}";
    angle_names[2] = "thetab";
    angle_title[2] = "#theta_{b}";

    if (max_entries > 0) {
        if (selection != "")
            selection = "(" + selection + TString::Format(") && Entry$ < %d", max_entries);
        else
            selection = TString::Format("Entry$ < %d", max_entries);
    }

    // Calculate aklm parameters
    printf("Reading data for %s\n", dataset.Data());
    Int_t entries = t_data->Draw("1", selection, "goff");
    if (entries <= 0) {
        printf("No events in dataset !\n");
        return;
    }
    printf("Found %d events\n", entries);
    t_data->SetEstimate(entries + 1);
    t_data->Draw(angle_names[2] + ":" + angle_names[1] + ":" + angle_names[0], selection, "goff");
    Double_t *thetab = t_data->GetV1();
    Double_t *thetat = t_data->GetV2();
    Double_t *phit = t_data->GetV3();

    printf("Calculating aklm parameters\n");
    std::map<UInt_t, std::map<UInt_t, std::map<Int_t, Double_t> > > aklm;
    std::vector<Double_t> aklm_abs_vector;
    Int_t nPar = 0;
    Int_t i;
    for (Int_t k = 0; k <= max_k; k++) {
        for (Int_t l = 0; l <= max_l; l++) {
            for (Int_t m = -l; m <= l; m++) {
                aklm[k][l][m] = 0;
                for (i = 0; i < entries; i++) {
                    aklm[k][l][m] += Ylm(l, m, thetat[i], phit[i]) * Pk(k, thetab[i]);
                }
                aklm[k][l][m] /= (Double_t)entries;  // normalization of the function to 1
                aklm_abs_vector.push_back(TMath::Abs(aklm[k][l][m]));
                printf(" - aklm[%d][%d][%d] = %g\n", k, l, m, aklm[k][l][m]);
                nPar++;
            }
        }
    }
    std::sort(aklm_abs_vector.begin(), aklm_abs_vector.end(), std::greater<Double_t>());
    printf(" ... got %d parameters\n", nPar);

    // Create appropriate 3D function
    printf("Creating the 3D function\n");
    TF3 *f3 = new TF3("f_SumAklmYlmPk", SumAklmYlmPk, constants::cuts::phit_low,
                      constants::cuts::phit_high, constants::cuts::thetat_low,
                      constants::cuts::thetat_high, constants::cuts::thetab_low,
                      constants::cuts::thetab_high, nPar + 3);
    f3->SetParameter(0, 1);
    f3->FixParameter(1, max_k);
    f3->FixParameter(2, max_l);

    i = 3;
    for (Int_t k = 0; k <= max_k; k++) {
        for (Int_t l = 0; l <= max_l; l++) {
            for (Int_t m = -l; m <= l; m++) {
                f3->SetParameter(i, aklm[k][l][m]);
                i++;
            }
        }
    }

    // Create 3D projection (optimal binning for 1D/2D/3D projections)
    printf("Creating 3D histogram\n");
    if (nBins3D == -1)
        nBins3D = TMath::Min(
            (Int_t)(TMath::Max(5., TMath::Power(entries / 100., 1. / 2. /*3.*/ /*1.*/)) + 0.1), 50);
    TH3D *h3 = new TH3D("h3", "h3", nBins3D, constants::cuts::phit_low, constants::cuts::phit_high,
                        nBins3D, constants::cuts::thetat_low, constants::cuts::thetat_high, nBins3D,
                        constants::cuts::thetab_low, constants::cuts::thetab_high);
    t_data->Draw(angle_names[2] + ":" + angle_names[1] + ":" + angle_names[0] + ">>+ " +
                     TString(h3->GetName()),
                 selection, "goff");
    h3->SetXTitle(angle_title[0]);
    h3->SetYTitle(angle_title[1]);
    h3->SetZTitle(angle_title[2]);
    f3->FixParameter(0, entries * 2. * TMath::Pi() / h3->GetNbinsX() * 2. / h3->GetNbinsY() * 2. /
                            h3->GetNbinsZ());
    if (opt_checkNorm) {
        printf("Checking normalization...\n");
        printf("Normalization check: %g %g\n",
               f3->Integral(constants::cuts::phit_low, constants::cuts::phit_high,
                            constants::cuts::thetat_low, constants::cuts::thetat_high,
                            constants::cuts::thetab_low, constants::cuts::thetab_high),
               h3->Integral("Width"));
    }

    // Test plots with full aklm set
    GetChi2NDF(h3, f3, dataset, "plot text", smooth3D);
    if (opt_noSearch) return;

    std::map<TString, Double_t> chi2;

    // Scan over biggest aklm parameters and find minimal needed set
    TCanvas *c_chi2 = new TCanvas("c_chi2", "c_chi2", 2 * 300, 7 * 150);
    c_chi2->Divide(2, 7);
    c_chi2->Draw();
    TGraph *g_chi2[7];
    TGraph *g_prob[7];
    for (Int_t jj = 0; jj < 7; jj++) {
        g_chi2[jj] = new TGraph(1);
        g_prob[jj] = new TGraph(1);
        c_chi2->cd(2 * jj + 1);
        g_chi2[jj]->Draw("ALP");
        c_chi2->cd(2 * jj + 2);
        g_prob[jj]->Draw("ALP");
        g_chi2[jj]->SetLineColor(kBlue);
        g_prob[jj]->SetLineColor(kRed);
        g_chi2[jj]->SetMarkerSize(0.5);
        g_prob[jj]->SetMarkerSize(0.5);
        g_chi2[jj]->GetXaxis()->SetTitle("Number of free a_{klm} parameters");
        g_prob[jj]->GetXaxis()->SetTitle("Number of free a_{klm} parameters");
    }
    g_chi2[0]->GetYaxis()->SetTitle("#chi^2/NDF for 1D " + angle_title[0]);
    g_chi2[1]->GetYaxis()->SetTitle("#chi^2/NDF for 1D " + angle_title[1]);
    g_chi2[2]->GetYaxis()->SetTitle("#chi^2/NDF for 1D " + angle_title[2]);
    g_chi2[3]->GetYaxis()->SetTitle("#chi^2/NDF for 2D " + angle_title[0] + " vs. " +
                                    angle_title[1]);
    g_chi2[4]->GetYaxis()->SetTitle("#chi^2/NDF for 2D " + angle_title[0] + " vs. " +
                                    angle_title[2]);
    g_chi2[5]->GetYaxis()->SetTitle("#chi^2/NDF for 2D " + angle_title[1] + " vs. " +
                                    angle_title[2]);
    g_chi2[6]->GetYaxis()->SetTitle("#chi^2/NDF for 3D");
    g_prob[0]->GetYaxis()->SetTitle("Probability for 1D " + angle_title[0]);
    g_prob[1]->GetYaxis()->SetTitle("Probability for 1D " + angle_title[1]);
    g_prob[2]->GetYaxis()->SetTitle("Probability for 1D " + angle_title[2]);
    g_prob[3]->GetYaxis()->SetTitle("Probability for 2D " + angle_title[0] + " vs. " +
                                    angle_title[1]);
    g_prob[4]->GetYaxis()->SetTitle("Probability for 2D " + angle_title[0] + " vs. " +
                                    angle_title[2]);
    g_prob[5]->GetYaxis()->SetTitle("Probability for 2D " + angle_title[1] + " vs. " +
                                    angle_title[2]);
    g_prob[6]->GetYaxis()->SetTitle("Probability for 3D");
    c_chi2->Modified();
    c_chi2->Update();

    UInt_t j = 0;
    UInt_t min_j = 0;
    UInt_t max_j = f3->GetNpar() - 3;
    Bool_t condition;

    for (Int_t iter = 0; iter < f3->GetNpar() - 3; iter++) {
        i = 3;
        for (Int_t k = 0; k <= max_k; k++) {
            for (Int_t l = 0; l <= max_l; l++) {
                for (Int_t m = -l; m <= l; m++) {
                    f3->ReleaseParameter(i);
                    if (TMath::Abs(aklm[k][l][m]) >= aklm_abs_vector[j])
                        f3->SetParameter(i, aklm[k][l][m]);
                    else
                        f3->FixParameter(i, 0);
                    i++;
                }
            }
        }
        for (nPar = 0, i = 3; i < f3->GetNpar(); i++)
            if (f3->GetParameter(i) != 0) nPar++;  // fix for buggy ROOT !
        printf("Interation %d (min %g): %d free parameters\n", iter, aklm_abs_vector[j],
               nPar /*f3->GetNumberFreeParameters()*/);
        chi2 = GetChi2NDF(h3, f3, dataset + "_iterations", "" /*"plot text"*/, smooth3D);
        g_chi2[0]->SetPoint(j, nPar + 1, chi2["chi2NDF_phi"]);
        g_chi2[1]->SetPoint(j, nPar + 1, chi2["chi2NDF_cosTheta1"]);
        g_chi2[2]->SetPoint(j, nPar + 1, chi2["chi2NDF_cosTheta2"]);
        g_chi2[3]->SetPoint(j, nPar + 1, chi2["chi2NDF_phi_cosTheta1"]);
        g_chi2[4]->SetPoint(j, nPar + 1, chi2["chi2NDF_phi_cosTheta2"]);
        g_chi2[5]->SetPoint(j, nPar + 1, chi2["chi2NDF_cosTheta1_cosTheta2"]);
        g_chi2[6]->SetPoint(j, nPar + 1, chi2["chi2NDF_phi_cosTheta1_cosTheta2"]);
        g_prob[0]->SetPoint(j, nPar + 1, chi2["prob_phi"]);
        g_prob[1]->SetPoint(j, nPar + 1, chi2["prob_cosTheta1"]);
        g_prob[2]->SetPoint(j, nPar + 1, chi2["prob_cosTheta2"]);
        g_prob[3]->SetPoint(j, nPar + 1, chi2["prob_phi_cosTheta1"]);
        g_prob[4]->SetPoint(j, nPar + 1, chi2["prob_phi_cosTheta2"]);
        g_prob[5]->SetPoint(j, nPar + 1, chi2["prob_cosTheta1_cosTheta2"]);
        g_prob[6]->SetPoint(j, nPar + 1, chi2["prob_phi_cosTheta1_cosTheta2"]);
        for (Int_t jj = 0; jj < 7; jj++) {
            c_chi2->cd(2 * jj + 1);
            g_chi2[jj]->Draw("ALP");
            c_chi2->cd(2 * jj + 2);
            g_prob[jj]->Draw("ALP");
        }
        c_chi2->Modified();
        c_chi2->Update();

        condition =  // chi2["prob_phi"]                 > 0.05 &&
                     // chi2["prob_cosTheta1"]           > 0.05 &&
                     // chi2["prob_cosTheta2"]           > 0.05 &&
                     // chi2["prob_phi_cosTheta1"]       > 0.05 &&
                     // chi2["prob_phi_cosTheta2"]       > 0.05 &&
                     // chi2["prob_cosTheta1_cosTheta2"] > 0.05;
                     // chi2["prob_phi_cosTheta1_cosTheta2"] > 0.05;
                     // chi2["chi2NDF_phi"]                 < 2 &&
                     // chi2["chi2NDF_cosTheta1"]           < 2 &&
                     // chi2["chi2NDF_cosTheta2"]           < 2 &&
            chi2["chi2NDF_phi_cosTheta1"] < 3 && chi2["chi2NDF_phi_cosTheta2"] < 3 &&
            chi2["chi2NDF_cosTheta1_cosTheta2"] < 3;

        if (!opt_fastSearch) {
            // Scan starting from the smallest number of parameters and increment by 1
            j++;
            if (condition) break;
        } else {
            // Jump by 2x there and back to find the optimum
            if (condition) {
                max_j = j;
                if (min_j == max_j) break;
                j -= (j + 1 - min_j) / 2;
            } else {
                min_j = j + 1;
                j += (max_j + 1 - j) / 2;
            }
        }
    }

    printf("Non-zero aklm parameters:\n");
    i = 3;
    for (Int_t k = 0; k <= max_k; k++) {
        for (Int_t l = 0; l <= max_l; l++) {
            for (Int_t m = -l; m <= l; m++) {
                if (f3->GetParameter(i) != 0) {
                    printf(" - a[k=%d][l=%d][m=%d] = %g\n", k, l, m, aklm[k][l][m]);
                }
                i++;
            }
        }
    }

    // Draw the final result of the optimization
    GetChi2NDF(h3, f3, dataset + "_minimal", "plot text", smooth3D);

    // Draw the final result with fine binning
    printf("Creating 3D histogram for fine 1D bining\n");
    if (nBins1D == -1) nBins1D = TMath::Min((Int_t)(TMath::Max(5., entries / 100.) + 0.1), 100);
    if (nBins1D != nBins3D) {
        TH3D *h3_fine = new TH3D("h3_fine", "h3_fine", nBins1D, constants::cuts::phit_low,
                                 constants::cuts::phit_high, nBins1D, constants::cuts::thetat_low,
                                 constants::cuts::thetat_high, nBins1D, constants::cuts::thetab_low,
                                 constants::cuts::thetab_high);
        t_data->Draw(angle_names[2] + ":" + angle_names[1] + ":" + angle_names[0] + ">>+ " +
                         TString(h3_fine->GetName()),
                     selection, "goff");
        h3_fine->SetXTitle(angle_title[0]);
        h3_fine->SetYTitle(angle_title[1]);
        h3_fine->SetZTitle(angle_title[2]);
        f3->FixParameter(0, entries * 2. * TMath::Pi() / h3_fine->GetNbinsX() * 2. /
                                h3_fine->GetNbinsY() * 2. / h3_fine->GetNbinsZ());
        if (opt_checkNorm) {
            printf("Checking normalization...\n");
            printf("Normalization check: %g %g\n",
                   f3->Integral(constants::cuts::phit_low, constants::cuts::phit_high,
                                constants::cuts::thetat_low, constants::cuts::thetat_high,
                                constants::cuts::thetab_low, constants::cuts::thetab_high),
                   h3->Integral("Width"));
        }
        GetChi2NDF(h3_fine, f3, dataset + "_minimal1Dopt", "plot text",
                   smooth1D);
    } else {
        printf(" - not needed, same binning defined for 3D and 1D\n");
    }
}

const void FitterSHE::LogTextFromFile(const char* field_name, const char* filename) {
    // Set the current directory back to the one for plots (ugly ROOT stuff)
    if (output_file_) {
        output_file_->cd();
    }
    std::ifstream file;
    file.open(filename);
    std::stringstream buffer;
    buffer << file.rdbuf();
    TNamed text(field_name, buffer.str());
    text.Write();
}

const void FitterSHE::LogFileCRC(const char* field_name, const char* filename) {
    char buffer[100];
    snprintf(buffer, 100, "%lu", cksum(filename, true));
    TNamed crc(field_name, buffer);
    crc.Write();
}

const void FitterSHE::LogText(const char* field_name, const char* text) {
    TNamed text_field(field_name, text);
    text_field.Write();
}

const void FitterSHE::LogCLIArguments(int argc, char* argv[]) {
    // Set the current directory back to the one for plots (ugly ROOT stuff)
    if (output_file_) {
        output_file_->cd();
    }

    std::string str;
    for (int i = 0; i < argc; i++) {
        str += argv[i];
        str += " ";
    }
    // Remove the final space
    str.pop_back();
    TNamed cli_arguments("cli_arguments", str.c_str());
    cli_arguments.Write();
}

const void FitterSHE::LogEnvironmentMetadata() {
    // Set the current directory back to the one for plots (ugly ROOT stuff)
    if (output_file_) {
        output_file_->cd();
    }

    TNamed root_version("root_version", ROOT_RELEASE);
    root_version.Write();

    char buffer[100];
    gethostname(buffer, 100);
    TNamed hostname("hostname", buffer);
    hostname.Write();

    const time_t now = time(0);
    const char* local_time_string = ctime(&now);
    TNamed local_date("local_date", local_time_string);

    tm* gmtm = gmtime(&now);
    const char* utc_time_string = asctime(gmtm);
    TNamed utc_date("utc_date", utc_time_string);

    TNamed git_version("git_version", gitversion);
    git_version.Write();
}

/**
 * Set the directory to which to ouput plots
 */
void FitterSHE::SetPlotDir(const char* plot_dir) {
    make_plots_ = true;
    tools::SetPlotDir(plot_dir);
}