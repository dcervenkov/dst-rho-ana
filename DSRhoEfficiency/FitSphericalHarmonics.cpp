#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TF3.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3.h"
#include "TH3D.h"
#include "TString.h"

#include "TMath.h"
#include "Math/SpecFunc.h"

#include "Rtypes.h"

#include <algorithm>
#include <map>
#include <vector>



Double_t Ylm(UInt_t l, Int_t m, Double_t cosTheta1, Double_t phi) {

  // https://en.wikipedia.org/wiki/Spherical_harmonics (real form)
  // https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
  Double_t theta1 = TMath::ACos(cosTheta1);
  if      ( m > 0 ) return TMath::Power(-1, m) * TMath::Sqrt2() * ROOT::Math::sph_legendre(l,  m, theta1) * TMath::Cos( m*phi); // integral is 1
  else if ( m < 0 ) return TMath::Power(-1, m) * TMath::Sqrt2() * ROOT::Math::sph_legendre(l, -m, theta1) * TMath::Sin(-m*phi); // integral is 1
  else              return ROOT::Math::sph_legendre(l, 0, theta1); // integral is 1
}

Double_t Pk(UInt_t k, Double_t cosTheta2) {

  // https://en.wikipedia.org/wiki/Legendre_polynomials
  return ROOT::Math::legendre(k, cosTheta2) * TMath::Sqrt(k+0.5); // integral is 1
}



Double_t SumAklmYlmPk(double *x, double *par) {

  Double_t phi       = x[0];
  Double_t cosTheta1 = x[1];
  Double_t cosTheta2 = x[2];

  //while ( phi < -TMath::Pi() ) phi += 2*TMath::Pi();
  //while ( phi > +TMath::Pi() ) phi -= 2*TMath::Pi();
  if ( cosTheta1 < -1 ) cosTheta1 = -1;
  if ( cosTheta1 > +1 ) cosTheta1 = +1;
  if ( cosTheta2 < -1 ) cosTheta2 = -1;
  if ( cosTheta2 > +1 ) cosTheta2 = +1;

  Double_t &norm  = par[0];
  Int_t max_k     = (Int_t)(par[1]+0.1);
  Int_t max_l     = (Int_t)(par[2]+0.1);

  Double_t result = 0;

  Int_t i = 3; // par[0],[1],[2] are constants: normalization, max_k and max_l
  for ( Int_t k=0; k<=max_k; k++ ) {
    for ( Int_t l=0; l<=max_l; l++ ) {
      for ( Int_t m=-l; m<=l; m++ ) {
        Double_t &aklm = par[i];
        result += aklm * Ylm(l, m, cosTheta1, phi) * Pk(k, cosTheta2);
        i += 1;
      }
    }
  }

  return norm*result;
}



std::map<TString,Double_t> GetChi2NDF(TH3D *h3_phi_cosTheta1_cosTheta2, TF3 *f3, TString title = "", TString options = "plot text", Int_t smoothScale = 100) {

  static Int_t counter = 0;

  Bool_t opt_useIntegral  = options.Contains("integ", TString::kIgnoreCase);
  Bool_t opt_controlPlots = options.Contains("plot" , TString::kIgnoreCase);
  Bool_t opt_textOutput   = options.Contains("text" , TString::kIgnoreCase);

  Int_t nbins_phi       = h3_phi_cosTheta1_cosTheta2->GetNbinsX();
  Int_t nbins_cosTheta1 = h3_phi_cosTheta1_cosTheta2->GetNbinsY();
  Int_t nbins_cosTheta2 = h3_phi_cosTheta1_cosTheta2->GetNbinsZ();

  // Create 3D histogram from the TF3 with binning "smoothScale"-times finer than the data histogram
  f3->SetNpx(h3_phi_cosTheta1_cosTheta2->GetNbinsX()*smoothScale);
  f3->SetNpy(h3_phi_cosTheta1_cosTheta2->GetNbinsY()*smoothScale);
  f3->SetNpz(h3_phi_cosTheta1_cosTheta2->GetNbinsZ()*smoothScale);
  TH3 *hf3_phi_cosTheta1_cosTheta2 = (TH3*)(f3->GetHistogram());
  Double_t hf3_nBinsX = hf3_phi_cosTheta1_cosTheta2->GetNbinsX();
  Double_t hf3_nBinsY = hf3_phi_cosTheta1_cosTheta2->GetNbinsY();
  Int_t    hf3_nBinsZ = hf3_phi_cosTheta1_cosTheta2->GetNbinsZ();
  for ( Int_t ix=1; ix<=hf3_phi_cosTheta1_cosTheta2->GetNbinsX(); ix++ ) {
    for ( Int_t iy=1; iy<=hf3_phi_cosTheta1_cosTheta2->GetNbinsY(); iy++ ) {
      if ( (iy%((100/hf3_nBinsZ)+1)) == 0 ) {
        printf("Evaluating function: %.5g%%          \r", 100.*((ix-1)*hf3_nBinsY+(iy-1))/(hf3_nBinsX*hf3_nBinsY));
        flush(cout);
      }
      for ( Int_t iz=1; iz<=hf3_phi_cosTheta1_cosTheta2->GetNbinsZ(); iz++ ) {
        double xpos = hf3_phi_cosTheta1_cosTheta2->GetXaxis()->GetBinCenter(ix);
        double ypos = hf3_phi_cosTheta1_cosTheta2->GetYaxis()->GetBinCenter(iy);
        double zpos = hf3_phi_cosTheta1_cosTheta2->GetZaxis()->GetBinCenter(iz);
        if ( opt_useIntegral ) hf3_phi_cosTheta1_cosTheta2->SetBinContent(ix, iy, iz, f3->Eval(xpos, ypos, zpos)); // TODO: replace eval by integral in the bin
        else                   hf3_phi_cosTheta1_cosTheta2->SetBinContent(ix, iy, iz, f3->Eval(xpos, ypos, zpos));
      }
    }
  }
  printf("                                     \r");
  flush(cout);

  // Create 1D and 2D projections from the data histogram
  TH1D *h1_phi                 = (TH1D*)h3_phi_cosTheta1_cosTheta2->Project3D("x");
  TH1D *h1_cosTheta1           = (TH1D*)h3_phi_cosTheta1_cosTheta2->Project3D("y");
  TH1D *h1_cosTheta2           = (TH1D*)h3_phi_cosTheta1_cosTheta2->Project3D("z");
  TH2D *h2_phi_cosTheta1       = (TH2D*)h3_phi_cosTheta1_cosTheta2->Project3D("yx");
  TH2D *h2_phi_cosTheta2       = (TH2D*)h3_phi_cosTheta1_cosTheta2->Project3D("zx");
  TH2D *h2_cosTheta1_cosTheta2 = (TH2D*)h3_phi_cosTheta1_cosTheta2->Project3D("zy");

  // Create 1D and 2D projections from the TF3 histogram
  TH1D *hf1_phi                 = (TH1D*)hf3_phi_cosTheta1_cosTheta2->Project3D("x");
  TH1D *hf1_cosTheta1           = (TH1D*)hf3_phi_cosTheta1_cosTheta2->Project3D("y");
  TH1D *hf1_cosTheta2           = (TH1D*)hf3_phi_cosTheta1_cosTheta2->Project3D("z");
  TH2D *hf2_phi_cosTheta1       = (TH2D*)hf3_phi_cosTheta1_cosTheta2->Project3D("yx");
  TH2D *hf2_phi_cosTheta2       = (TH2D*)hf3_phi_cosTheta1_cosTheta2->Project3D("zx");
  TH2D *hf2_cosTheta1_cosTheta2 = (TH2D*)hf3_phi_cosTheta1_cosTheta2->Project3D("zy");
  hf1_phi                ->Scale(1./smoothScale/smoothScale);
  hf1_cosTheta1          ->Scale(1./smoothScale/smoothScale);
  hf1_cosTheta2          ->Scale(1./smoothScale/smoothScale);
  hf2_phi_cosTheta1      ->Scale(1./smoothScale);
  hf2_phi_cosTheta2      ->Scale(1./smoothScale);
  hf2_cosTheta1_cosTheta2->Scale(1./smoothScale);

  // Calculate chi2/NDF for all the projections and fill-up the histograms of residuals
  TH1D hres_phi                    ("hres_phi"                    , "hres_phi"                    , 200, -10, 10);
  TH1D hres_cosTheta1              ("hres_cosTheta1"              , "hres_cosTheta1"              , 200, -10, 10);
  TH1D hres_cosTheta2              ("hres_cosTheta2"              , "hres_cosTheta2"              , 200, -10, 10);
  TH1D hres_phi_cosTheta1          ("hres_phi_cosTheta1"          , "hres_phi_cosTheta1"          , 200, -10, 10);
  TH1D hres_phi_cosTheta2          ("hres_phi_cosTheta2"          , "hres_phi_cosTheta2"          , 200, -10, 10);
  TH1D hres_cosTheta1_cosTheta2    ("hres_cosTheta1_cosTheta2"    , "hres_cosTheta1_cosTheta2"    , 200, -10, 10);
  TH1D hres_phi_cosTheta1_cosTheta2("hres_phi_cosTheta1_cosTheta2", "hres_phi_cosTheta1_cosTheta2", 200, -10, 10);
  Int_t NDF_phi                      = nbins_phi                                 /*- f3->GetNumberFreeParameters()*/;
  Int_t NDF_cosTheta1                = nbins_cosTheta1                           /*- f3->GetNumberFreeParameters()*/;
  Int_t NDF_cosTheta2                = nbins_cosTheta2                           /*- f3->GetNumberFreeParameters()*/;
  Int_t NDF_phi_cosTheta1            = nbins_phi*nbins_cosTheta1                 /*- f3->GetNumberFreeParameters()*/;
  Int_t NDF_phi_cosTheta2            = nbins_phi*nbins_cosTheta2                 /*- f3->GetNumberFreeParameters()*/;
  Int_t NDF_cosTheta1_cosTheta2      = nbins_cosTheta1*nbins_cosTheta2           /*- f3->GetNumberFreeParameters()*/;
  Int_t NDF_phi_cosTheta1_cosTheta2  = nbins_phi*nbins_cosTheta1*nbins_cosTheta2 /*- f3->GetNumberFreeParameters()*/;
  Double_t chi2_phi                      = 0;
  Double_t chi2_cosTheta1                = 0;
  Double_t chi2_cosTheta2                = 0;
  Double_t chi2_phi_cosTheta1            = 0;
  Double_t chi2_phi_cosTheta2            = 0;
  Double_t chi2_cosTheta1_cosTheta2      = 0;
  Double_t chi2_phi_cosTheta1_cosTheta2  = 0;
  Double_t res = 0;
  // 1D
  for ( Int_t ix=1 ; ix<=nbins_phi ; ix++ ) {
    res = 0;
    for ( Int_t ifx=1+(ix-1)*smoothScale; ifx<1+ix*smoothScale; ifx++ ) {
      res -= hf1_phi->GetBinContent(ifx);
    }
    res /= smoothScale;
    res += h1_phi->GetBinContent(ix);
    if ( h1_phi->GetBinError(ix) > 0 ) {
      chi2_phi += res*res / TMath::Sq(h1_phi->GetBinError(ix));
      hres_phi.Fill(res/h1_phi->GetBinError(ix));
    }
  }
  for ( Int_t ix=1 ; ix<=nbins_cosTheta1 ; ix++ ) {
    res = 0;
    for ( Int_t ifx=1+(ix-1)*smoothScale; ifx<1+ix*smoothScale; ifx++ ) {
      res -= hf1_cosTheta1->GetBinContent(ifx);
    }
    res /= smoothScale;
    res += h1_cosTheta1->GetBinContent(ix);
    if ( h1_cosTheta1->GetBinError(ix) > 0 ) {
      chi2_cosTheta1 += res*res / TMath::Sq(h1_cosTheta1->GetBinError(ix));
      hres_cosTheta1.Fill(res/h1_cosTheta1->GetBinError(ix));
    }
  }
  for ( Int_t ix=1 ; ix<=nbins_cosTheta2 ; ix++ ) {
    res = 0;
    for ( Int_t ifx=1+(ix-1)*smoothScale; ifx<1+ix*smoothScale; ifx++ ) {
      res -= hf1_cosTheta2->GetBinContent(ifx);
    }
    res /= smoothScale;
    res += h1_cosTheta2->GetBinContent(ix);
    if ( h1_cosTheta2->GetBinError(ix) > 0 ) {
      chi2_cosTheta2 += res*res / TMath::Sq(h1_cosTheta2->GetBinError(ix));
      hres_cosTheta2.Fill(res/h1_cosTheta2->GetBinError(ix));
    }
  }
  // 2D
  for (   Int_t ix=1 ; ix<=nbins_phi       ; ix++ ) {
    for ( Int_t iy=1 ; iy<=nbins_cosTheta1 ; iy++ ) {
      res = 0;
      for (   Int_t ifx=1+(ix-1)*smoothScale; ifx<1+ix*smoothScale; ifx++ ) {
        for ( Int_t ify=1+(iy-1)*smoothScale; ify<1+iy*smoothScale; ify++ ) {
          res -= hf2_phi_cosTheta1->GetBinContent(ifx, ify);
        }
      }
      res /= smoothScale*smoothScale;
      res += h2_phi_cosTheta1->GetBinContent(ix, iy);
      if ( h2_phi_cosTheta1->GetBinError(ix, iy) > 0 ) {
        chi2_phi_cosTheta1 += res*res / TMath::Sq(h2_phi_cosTheta1->GetBinError(ix, iy));
        hres_phi_cosTheta1.Fill(res/h2_phi_cosTheta1->GetBinError(ix, iy));
      }
    }
  }
  for (   Int_t ix=1 ; ix<=nbins_phi       ; ix++ ) {
    for ( Int_t iy=1 ; iy<=nbins_cosTheta2 ; iy++ ) {
      res = 0;
      for (   Int_t ifx=1+(ix-1)*smoothScale; ifx<1+ix*smoothScale; ifx++ ) {
        for ( Int_t ify=1+(iy-1)*smoothScale; ify<1+iy*smoothScale; ify++ ) {
          res -= hf2_phi_cosTheta2->GetBinContent(ifx, ify);
        }
      }
      res /= smoothScale*smoothScale;
      res += h2_phi_cosTheta2->GetBinContent(ix, iy);
      if ( h2_phi_cosTheta2->GetBinError(ix, iy) > 0 ) {
        chi2_phi_cosTheta2 += res*res / TMath::Sq(h2_phi_cosTheta2->GetBinError(ix, iy));
        hres_phi_cosTheta2.Fill(res/h2_phi_cosTheta2->GetBinError(ix, iy));
      }
    }
  }
  for (   Int_t ix=1 ; ix<=nbins_cosTheta1 ; ix++ ) {
    for ( Int_t iy=1 ; iy<=nbins_cosTheta2 ; iy++ ) {
      res = 0;
      for (   Int_t ifx=1+(ix-1)*smoothScale; ifx<1+ix*smoothScale; ifx++ ) {
        for ( Int_t ify=1+(iy-1)*smoothScale; ify<1+iy*smoothScale; ify++ ) {
          res -= hf2_cosTheta1_cosTheta2->GetBinContent(ifx, ify);
        }
      }
      res /= smoothScale*smoothScale;
      res += h2_cosTheta1_cosTheta2->GetBinContent(ix, iy);
      if ( h2_cosTheta1_cosTheta2->GetBinError(ix, iy) > 0 ) {
        chi2_cosTheta1_cosTheta2 += res*res / TMath::Sq(h2_cosTheta1_cosTheta2->GetBinError(ix, iy));
        hres_cosTheta1_cosTheta2.Fill(res/h2_cosTheta1_cosTheta2->GetBinError(ix, iy));
      }
    }
  }
  // 3D
  for (     Int_t ix=1 ; ix<=nbins_phi       ; ix++ ) {
    for (   Int_t iy=1 ; iy<=nbins_cosTheta1 ; iy++ ) {
      for ( Int_t iz=1 ; iz<=nbins_cosTheta2 ; iz++ ) {
        res = 0;
        for (     Int_t ifx=1+(ix-1)*smoothScale; ifx<1+ix*smoothScale; ifx++ ) {
          for (   Int_t ify=1+(iy-1)*smoothScale; ify<1+iy*smoothScale; ify++ ) {
            for ( Int_t ifz=1+(iz-1)*smoothScale; ifz<1+iz*smoothScale; ifz++ ) {
              res -= hf3_phi_cosTheta1_cosTheta2->GetBinContent(ifx, ify, ifz);
            }
          }
        }
        res /= smoothScale*smoothScale*smoothScale;
        res += h3_phi_cosTheta1_cosTheta2->GetBinContent(ix, iy, iz);
        if ( h3_phi_cosTheta1_cosTheta2->GetBinError(ix, iy, iz) > 0 ) {
          chi2_phi_cosTheta1_cosTheta2 += res*res / TMath::Sq(h3_phi_cosTheta1_cosTheta2->GetBinError(ix, iy, iz));
          hres_phi_cosTheta1_cosTheta2.Fill(res/h3_phi_cosTheta1_cosTheta2->GetBinError(ix, iy, iz));
        }
      }
    }
  }

  // Print the chi2/NDF
  TString text_phi                     = TString::Format("chi2/NDF = %.3g = %.5g / %d", chi2_phi                    /NDF_phi,
                                                                                        chi2_phi                    ,NDF_phi);
  TString text_cosTheta1               = TString::Format("chi2/NDF = %.3g = %.5g / %d", chi2_cosTheta1              /NDF_cosTheta1,
                                                                                        chi2_cosTheta1              ,NDF_cosTheta1);
  TString text_cosTheta2               = TString::Format("chi2/NDF = %.3g = %.5g / %d", chi2_cosTheta2              /NDF_cosTheta2,
                                                                                        chi2_cosTheta2              ,NDF_cosTheta2);
  TString text_phi_cosTheta1           = TString::Format("chi2/NDF = %.3g = %.5g / %d", chi2_phi_cosTheta1          /NDF_phi_cosTheta1,
                                                                                        chi2_phi_cosTheta1          ,NDF_phi_cosTheta1);
  TString text_phi_cosTheta2           = TString::Format("chi2/NDF = %.3g = %.5g / %d", chi2_phi_cosTheta2          /NDF_phi_cosTheta2,
                                                                                        chi2_phi_cosTheta2          ,NDF_phi_cosTheta2);
  TString text_cosTheta1_cosTheta2     = TString::Format("chi2/NDF = %.3g = %.5g / %d", chi2_cosTheta1_cosTheta2    /NDF_cosTheta1_cosTheta2,
                                                                                        chi2_cosTheta1_cosTheta2    ,NDF_cosTheta1_cosTheta2);
  TString text_phi_cosTheta1_cosTheta2 = TString::Format("chi2/NDF = %.3g = %.5g / %d", chi2_phi_cosTheta1_cosTheta2/NDF_phi_cosTheta1_cosTheta2,
                                                                                        chi2_phi_cosTheta1_cosTheta2,NDF_phi_cosTheta1_cosTheta2);

  if ( opt_textOutput ) {
    printf("1D phi      :  %s\n"                  , text_phi                    .Data());
    printf("1D cosTheta1:  %s\n"                  , text_cosTheta1              .Data());
    printf("1D cosTheta2:  %s\n"                  , text_cosTheta2              .Data());
    printf("2D phi       x cosTheta1:  %s\n"      , text_phi_cosTheta1          .Data());
    printf("2D phi       x cosTheta2:  %s\n"      , text_phi_cosTheta2          .Data());
    printf("2D cosTheta1 x cosTheta2:  %s\n"      , text_cosTheta1_cosTheta2    .Data());
    printf("3D phi x cosTheta1 x cosTheta2:  %s\n", text_phi_cosTheta1_cosTheta2.Data());
  }

  hres_phi                    .SetXTitle(text_phi);
  hres_cosTheta1              .SetXTitle(text_cosTheta1);
  hres_cosTheta2              .SetXTitle(text_cosTheta2);
  hres_phi_cosTheta1          .SetXTitle(text_phi_cosTheta1);
  hres_phi_cosTheta2          .SetXTitle(text_phi_cosTheta2);
  hres_cosTheta1_cosTheta2    .SetXTitle(text_cosTheta1_cosTheta2);
  hres_phi_cosTheta1_cosTheta2.SetXTitle(text_phi_cosTheta1_cosTheta2);

  // Draw the projections and residuals
  if ( opt_controlPlots ) {
    Double_t max_phi                     = h1_phi                    ->GetMaximum()*1.2;
    Double_t max_cosTheta1               = h1_cosTheta1              ->GetMaximum()*1.2;
    Double_t max_cosTheta2               = h1_cosTheta2              ->GetMaximum()*1.2;
    Double_t max_phi_cosTheta1           = h2_phi_cosTheta1          ->GetMaximum()*1.2;
    Double_t max_phi_cosTheta2           = h2_phi_cosTheta2          ->GetMaximum()*1.2;
    Double_t max_cosTheta1_cosTheta2     = h2_cosTheta1_cosTheta2    ->GetMaximum()*1.2;
    TCanvas *c = new TCanvas("c_chi2ndf_"+TString::Format("%d", counter), "Chi2/NDF and projections: "+title, 7*200, 3*200);
    c->Draw();
    c->Divide(7, 3);
    c->cd(1   ); h1_phi                      ->DrawCopy("hist")->GetYaxis()->SetRangeUser(0, max_phi);
    c->cd(1+7 ); hf1_phi                     ->DrawCopy("hist")->GetYaxis()->SetRangeUser(0, max_phi);
    c->cd(1+14); hres_phi                     .DrawCopy("");
    c->cd(2   ); h1_cosTheta1                ->DrawCopy("hist")->GetYaxis()->SetRangeUser(0, max_cosTheta1);
    c->cd(2+7 ); hf1_cosTheta1               ->DrawCopy("hist")->GetYaxis()->SetRangeUser(0, max_cosTheta1);
    c->cd(2+14); hres_cosTheta1               .DrawCopy("");
    c->cd(3   ); h1_cosTheta2                ->DrawCopy("hist")->GetYaxis()->SetRangeUser(0, max_cosTheta2);
    c->cd(3+7 ); hf1_cosTheta2               ->DrawCopy("hist")->GetYaxis()->SetRangeUser(0, max_cosTheta2);
    c->cd(3+14); hres_cosTheta2               .DrawCopy("");
    c->cd(4   ); h2_phi_cosTheta1            ->DrawCopy("colz")->GetZaxis()->SetRangeUser(0, max_phi_cosTheta1);
    c->cd(4+7 ); hf2_phi_cosTheta1           ->DrawCopy("colz")->GetZaxis()->SetRangeUser(0, max_phi_cosTheta1);
    c->cd(4+14); hres_phi_cosTheta1           .DrawCopy("");
    c->cd(5   ); h2_phi_cosTheta2            ->DrawCopy("colz")->GetZaxis()->SetRangeUser(0, max_phi_cosTheta2);
    c->cd(5+7 ); hf2_phi_cosTheta2           ->DrawCopy("colz")->GetZaxis()->SetRangeUser(0, max_phi_cosTheta2);
    c->cd(5+14); hres_phi_cosTheta2           .DrawCopy("");
    c->cd(6   ); h2_cosTheta1_cosTheta2      ->DrawCopy("colz")->GetZaxis()->SetRangeUser(0, max_cosTheta1_cosTheta2);
    c->cd(6+7 ); hf2_cosTheta1_cosTheta2     ->DrawCopy("colz")->GetZaxis()->SetRangeUser(0, max_cosTheta1_cosTheta2);
    c->cd(6+14); hres_cosTheta1_cosTheta2     .DrawCopy("");
    c->cd(7   ); h3_phi_cosTheta1_cosTheta2  ->DrawCopy("lego");
    c->cd(7+7 ); hf3_phi_cosTheta1_cosTheta2 ->DrawCopy("lego");
    c->cd(7+14); hres_phi_cosTheta1_cosTheta2 .DrawCopy("");
    c->Modified();
    c->Update();
  }

  std::map<TString,Double_t> result;
  result[    "NDF_phi"]                     =  NDF_phi;
  result[   "chi2_phi"]                     = chi2_phi;
  result["chi2NDF_phi"]                     = chi2_phi / NDF_phi;
  result[   "prob_phi"]                     = TMath::Prob(chi2_phi, NDF_phi);
  result[    "NDF_cosTheta1"]               =  NDF_cosTheta1;
  result[   "chi2_cosTheta1"]               = chi2_cosTheta1;
  result["chi2NDF_cosTheta1"]               = chi2_cosTheta1 / NDF_cosTheta1;
  result[   "prob_cosTheta1"]               = TMath::Prob(chi2_cosTheta1, NDF_cosTheta1);
  result[    "NDF_cosTheta2"]               =  NDF_cosTheta2;
  result[   "chi2_cosTheta2"]               = chi2_cosTheta2;
  result["chi2NDF_cosTheta2"]               = chi2_cosTheta2 / NDF_cosTheta2;
  result[   "prob_cosTheta2"]               = TMath::Prob(chi2_cosTheta2, NDF_cosTheta2);
  result[    "NDF_phi_cosTheta1"]           =  NDF_phi_cosTheta1;
  result[   "chi2_phi_cosTheta1"]           = chi2_phi_cosTheta1;
  result["chi2NDF_phi_cosTheta1"]           = chi2_phi_cosTheta1 / NDF_phi_cosTheta1;
  result[   "prob_phi_cosTheta1"]           = TMath::Prob(chi2_phi_cosTheta1, NDF_phi_cosTheta1);
  result[    "NDF_phi_cosTheta2"]           =  NDF_phi_cosTheta2;
  result[   "chi2_phi_cosTheta2"]           = chi2_phi_cosTheta2;
  result["chi2NDF_phi_cosTheta2"]           = chi2_phi_cosTheta2 / NDF_phi_cosTheta2;
  result[   "prob_phi_cosTheta2"]           = TMath::Prob(chi2_phi_cosTheta2, NDF_phi_cosTheta2);
  result[    "NDF_cosTheta1_cosTheta2"]     =  NDF_cosTheta1_cosTheta2;
  result[   "chi2_cosTheta1_cosTheta2"]     = chi2_cosTheta1_cosTheta2;
  result["chi2NDF_cosTheta1_cosTheta2"]     = chi2_cosTheta1_cosTheta2 / NDF_cosTheta1_cosTheta2;
  result[   "prob_cosTheta1_cosTheta2"]     = TMath::Prob(chi2_cosTheta1_cosTheta2, NDF_cosTheta1_cosTheta2);
  result[    "NDF_phi_cosTheta1_cosTheta2"] =  NDF_phi_cosTheta1_cosTheta2;
  result[   "chi2_phi_cosTheta1_cosTheta2"] = chi2_phi_cosTheta1_cosTheta2;
  result["chi2NDF_phi_cosTheta1_cosTheta2"] = chi2_phi_cosTheta1_cosTheta2 / NDF_phi_cosTheta1_cosTheta2;
  result[   "prob_phi_cosTheta1_cosTheta2"] = TMath::Prob(chi2_phi_cosTheta1_cosTheta2, NDF_phi_cosTheta1_cosTheta2);

  counter++;
  return result;
}



void AnalyzeDataset(TString dataset = "BsJpsiPhi", Int_t max_entries = -1,
                    Int_t nBins3D = 20, Int_t nBins1D = 50, Int_t smooth3D = 1, Int_t smooth1D = 1,
                    Int_t max_k = 10, Int_t max_l = 10, TString options = "") {

  Bool_t opt_checkNorm  = options.Contains("norm"    , TString::kIgnoreCase);
  Bool_t opt_fastSearch = options.Contains("fast"    , TString::kIgnoreCase);
  Bool_t opt_noSearch   = options.Contains("fullonly", TString::kIgnoreCase);

  // Read input ntuple, define angles titles for different datasets
  TChain  *t_data;
  TString  selection;
  TString  angle_names[3];
  TString  angle_title[3];

  if ( dataset == "BsJpsiPhi" ) {

    selection = "";
    t_data = new TChain("resolvedDup");
    t_data->Add("/home/reznicek/data/BsJpsiPhi/MC11/Acceptance/mu4mu4MoreStatsforAccpbb/*.root");
    angle_names[0] = "CDFPhi"     ; angle_title[0] = "#phi";       // phi must be 1st !!!, cos theta1/2 to follow
    angle_names[1] = "CDFCostheta"; angle_title[1] = "cos #theta";
    angle_names[2] = "CDFCosPsi"  ; angle_title[2] = "cos #psi";
    //nBins3D = 20;
    //nBins1D = 100;

  } else if ( dataset == "BdKstarMuMu_Bd_SMangles_208446" ) {

    selection = "";
    t_data = new TChain("Bd");
    t_data->Add("/home/reznicek/data/SemileptonicRareB/Data2012/MCntuples/cut/skim_208446_sel_red_cut.root");
    angle_names[0] = "Extra_helicityAnglePhi" ; angle_title[0] = "#phi";       // phi must be 1st !!!, cos theta1/2 to follow
    angle_names[1] = "Extra_helicityCosThetaL"; angle_title[1] = "cos #theta_{L}";
    angle_names[2] = "Extra_helicityCosThetaK"; angle_title[2] = "cos #theta_{K}";

  } else if ( dataset == "BsPhiMuMu_SMangles_208441" ) {

    selection = "";
    t_data = new TChain("Bd");
    t_data->Add("/home/reznicek/data/SemileptonicRareB/Data2012/MCntuples/cut/skim_208441_sel_red_cut.root");
    angle_names[0] = "Extra_helicityAnglePhi" ; angle_title[0] = "#phi";       // phi must be 1st !!!, cos theta1/2 to follow
    angle_names[1] = "Extra_helicityCosThetaL"; angle_title[1] = "cos #theta_{L}";
    angle_names[2] = "Extra_helicityCosThetaK"; angle_title[2] = "cos #theta_{K}";
    nBins3D = 15;
    nBins1D = 50;
    smooth3D = 3;
    smooth1D = 1;

  } else if ( dataset == "LambdabL0MuMu_SMangles_208456" ) {

    selection = "";
    t_data = new TChain("Bd");
    t_data->Add("/home/reznicek/data/SemileptonicRareB/Data2012/MCntuples/cut/skim_208456_sel_red_cut.root");
    angle_names[0] = "Extra_helicityAnglePhi" ; angle_title[0] = "#phi";       // phi must be 1st !!!, cos theta1/2 to follow
    angle_names[1] = "Extra_helicityCosThetaL"; angle_title[1] = "cos #theta_{L}";
    angle_names[2] = "Extra_helicityCosThetaK"; angle_title[2] = "cos #theta_{K}";

  } else if ( dataset == "LambdabpKMuMu_SMangles_208457" ) {

    selection = "";
    t_data = new TChain("Bd");
    t_data->Add("/home/reznicek/data/SemileptonicRareB/Data2012/MCntuples/cut/skim_208457_sel_red_cut.root");
    angle_names[0] = "Extra_helicityAnglePhi" ; angle_title[0] = "#phi";       // phi must be 1st !!!, cos theta1/2 to follow
    angle_names[1] = "Extra_helicityCosThetaL"; angle_title[1] = "cos #theta_{L}";
    angle_names[2] = "Extra_helicityCosThetaK"; angle_title[2] = "cos #theta_{K}";

  } else if ( dataset == "Data2012_Bsidebands" ) {

    selection = "Bd_VTX_mass < 5079 || Bd_VTX_mass > 5479";
    t_data = new TChain("Bd");
    t_data->Add("/home/reznicek/data/SemileptonicRareB/Data2012/cut/data12_v5.root");
    angle_names[0] = "Extra_helicityAnglePhi" ; angle_title[0] = "#phi";       // phi must be 1st !!!, cos theta1/2 to follow
    angle_names[1] = "Extra_helicityCosThetaL"; angle_title[1] = "cos #theta_{L}";
    angle_names[2] = "Extra_helicityCosThetaK"; angle_title[2] = "cos #theta_{K}";

  } else {
    printf("Unknown dataset %s !\n", dataset.Data());
    return;
  }

  if ( max_entries > 0 ) {
    if ( selection != "" ) selection = "("+selection+TString::Format(") && Entry$ < %d", max_entries);
    else selection = TString::Format("Entry$ < %d", max_entries);
  }

  // Calculate aklm parameters
  printf("Reading data for %s\n", dataset.Data());
  Int_t entries = t_data->Draw("1", selection, "goff");
  if ( entries <= 0 ) {
    printf("No events in dataset !\n");
    return;
  }
  printf("Found %d events\n", entries);
  t_data->SetEstimate(entries+1);
  t_data->Draw(angle_names[2]+":"+angle_names[1]+":"+angle_names[0], selection, "goff");
  Double_t *cosTheta2 = t_data->GetV1();
  Double_t *cosTheta1 = t_data->GetV2();
  Double_t *phi       = t_data->GetV3();

  printf("Calculating aklm parameters\n");
  std::map<UInt_t,std::map<UInt_t,std::map<Int_t, Double_t> > > aklm;
  std::vector<Double_t> aklm_abs_vector;
  Int_t nPar = 0;
  Int_t i;
  for ( Int_t k=0; k<=max_k; k++ ) {
    for ( Int_t l=0; l<=max_l; l++ ) {
      for ( Int_t m=-l; m<=l; m++ ) {
        aklm[k][l][m] = 0;
        for ( i=0; i<entries; i++ ) {
          aklm[k][l][m] += Ylm(l, m, cosTheta1[i], phi[i]) * Pk(k, cosTheta2[i]);
        }
        aklm[k][l][m] /= (Double_t)entries; // normalization of the function to 1
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
  TF3 *f3 = new TF3("f_SumAklmYlmPk", SumAklmYlmPk, -TMath::Pi(),TMath::Pi(), -1, 1, -1, 1, nPar + 3);
  f3->SetParameter(0, 1);
  f3->FixParameter(1, max_k);
  f3->FixParameter(2, max_l);

  i = 3;
  for ( Int_t k=0; k<=max_k; k++ ) {
    for ( Int_t l=0; l<=max_l; l++ ) {
      for ( Int_t m=-l; m<=l; m++ ) {
        f3->SetParameter(i, aklm[k][l][m]);
        i++;
      }
    }
  }

  // Create 3D projection (optimal binning for 1D/2D/3D projections)
  printf("Creating 3D histogram\n");
  if ( nBins3D == -1 ) nBins3D = TMath::Min((Int_t)(TMath::Max(5., TMath::Power(entries/100., 1./2./*3.*//*1.*/)) + 0.1), 50);
  TH3D *h3 = new TH3D("h3", "h3", nBins3D, -TMath::Pi(), TMath::Pi(), nBins3D, -1, 1, nBins3D, -1, 1);
  t_data->Draw(angle_names[2]+":"+angle_names[1]+":"+angle_names[0] + ">>+ " + TString(h3->GetName()), selection, "goff");
  h3->SetXTitle(angle_title[0]);
  h3->SetYTitle(angle_title[1]);
  h3->SetZTitle(angle_title[2]);
  f3->FixParameter(0, entries * 2.*TMath::Pi()/h3->GetNbinsX() * 2./h3->GetNbinsY() * 2./h3->GetNbinsZ());
  if ( opt_checkNorm ) {
    printf("Checking normalization...\n");
    printf("Normalization check: %g %g\n", f3->Integral(-TMath::Pi(), TMath::Pi(), -1, 1, -1, 1), h3->Integral("Width"));
  }

  // Test plots with full aklm set
  GetChi2NDF(h3, f3, "full aklm set", "plot text", smooth3D);
  if ( opt_noSearch ) return;

  std::map<TString,Double_t> chi2;

  // Scan over biggest aklm parameters and find minimal needed set
  TCanvas *c_chi2 = new TCanvas("c_chi2", "c_chi2", 2*300, 7*150);
  c_chi2->Divide(2, 7);
  c_chi2->Draw();
  TGraph *g_chi2[7];
  TGraph *g_prob[7];
  for ( Int_t jj=0; jj<7; jj++ ) {
    g_chi2[jj] = new TGraph(1);
    g_prob[jj] = new TGraph(1);
    c_chi2->cd(2*jj+1); g_chi2[jj]->Draw("ALP");
    c_chi2->cd(2*jj+2); g_prob[jj]->Draw("ALP");
    g_chi2[jj]->SetLineColor(kBlue);
    g_prob[jj]->SetLineColor(kRed);
    g_chi2[jj]->SetMarkerSize(0.5);
    g_prob[jj]->SetMarkerSize(0.5);
    g_chi2[jj]->GetXaxis()->SetTitle("Number of free a_{klm} parameters");
    g_prob[jj]->GetXaxis()->SetTitle("Number of free a_{klm} parameters");
  }
  g_chi2[0]->GetYaxis()->SetTitle("#chi^2/NDF for 1D "+angle_title[0]);
  g_chi2[1]->GetYaxis()->SetTitle("#chi^2/NDF for 1D "+angle_title[1]);
  g_chi2[2]->GetYaxis()->SetTitle("#chi^2/NDF for 1D "+angle_title[2]);
  g_chi2[3]->GetYaxis()->SetTitle("#chi^2/NDF for 2D "+angle_title[0]+" vs. "+angle_title[1]);
  g_chi2[4]->GetYaxis()->SetTitle("#chi^2/NDF for 2D "+angle_title[0]+" vs. "+angle_title[2]);
  g_chi2[5]->GetYaxis()->SetTitle("#chi^2/NDF for 2D "+angle_title[1]+" vs. "+angle_title[2]);
  g_chi2[6]->GetYaxis()->SetTitle("#chi^2/NDF for 3D");
  g_prob[0]->GetYaxis()->SetTitle("Probability for 1D "+angle_title[0]);
  g_prob[1]->GetYaxis()->SetTitle("Probability for 1D "+angle_title[1]);
  g_prob[2]->GetYaxis()->SetTitle("Probability for 1D "+angle_title[2]);
  g_prob[3]->GetYaxis()->SetTitle("Probability for 2D "+angle_title[0]+" vs. "+angle_title[1]);
  g_prob[4]->GetYaxis()->SetTitle("Probability for 2D "+angle_title[0]+" vs. "+angle_title[2]);
  g_prob[5]->GetYaxis()->SetTitle("Probability for 2D "+angle_title[1]+" vs. "+angle_title[2]);
  g_prob[6]->GetYaxis()->SetTitle("Probability for 3D");
  c_chi2->Modified();
  c_chi2->Update();

  UInt_t j     = 0;
  UInt_t min_j = 0;
  UInt_t max_j = f3->GetNpar()-3;
  Bool_t condition;

  for ( Int_t iter=0; iter < f3->GetNpar()-3; iter++ ) {

    i = 3;
    for ( Int_t k=0; k<=max_k; k++ ) {
      for ( Int_t l=0; l<=max_l; l++ ) {
        for ( Int_t m=-l; m<=l; m++ ) {
          f3->ReleaseParameter(i);
          if ( TMath::Abs(aklm[k][l][m]) >= aklm_abs_vector[j] ) f3->SetParameter(i, aklm[k][l][m]);
          else                                                   f3->FixParameter(i, 0);
          i++;
        }
      }
    }
    for ( nPar=0, i=3; i<f3->GetNpar(); i++ ) if ( f3->GetParameter(i) != 0 ) nPar++; // fix for buggy ROOT !
    printf("Interation %d (min %g): %d free parameters\n", iter, aklm_abs_vector[j], nPar/*f3->GetNumberFreeParameters()*/);
    chi2 = GetChi2NDF(h3, f3, "iterations", ""/*"plot text"*/, smooth3D);
    g_chi2[0]->SetPoint(j, nPar+1, chi2["chi2NDF_phi"]);
    g_chi2[1]->SetPoint(j, nPar+1, chi2["chi2NDF_cosTheta1"]);
    g_chi2[2]->SetPoint(j, nPar+1, chi2["chi2NDF_cosTheta2"]);
    g_chi2[3]->SetPoint(j, nPar+1, chi2["chi2NDF_phi_cosTheta1"]);
    g_chi2[4]->SetPoint(j, nPar+1, chi2["chi2NDF_phi_cosTheta2"]);
    g_chi2[5]->SetPoint(j, nPar+1, chi2["chi2NDF_cosTheta1_cosTheta2"]);
    g_chi2[6]->SetPoint(j, nPar+1, chi2["chi2NDF_phi_cosTheta1_cosTheta2"]);
    g_prob[0]->SetPoint(j, nPar+1, chi2["prob_phi"]);
    g_prob[1]->SetPoint(j, nPar+1, chi2["prob_cosTheta1"]);
    g_prob[2]->SetPoint(j, nPar+1, chi2["prob_cosTheta2"]);
    g_prob[3]->SetPoint(j, nPar+1, chi2["prob_phi_cosTheta1"]);
    g_prob[4]->SetPoint(j, nPar+1, chi2["prob_phi_cosTheta2"]);
    g_prob[5]->SetPoint(j, nPar+1, chi2["prob_cosTheta1_cosTheta2"]);
    g_prob[6]->SetPoint(j, nPar+1, chi2["prob_phi_cosTheta1_cosTheta2"]);
    for ( Int_t jj=0; jj<7; jj++ ) {
      c_chi2->cd(2*jj+1); g_chi2[jj]->Draw("ALP");
      c_chi2->cd(2*jj+2); g_prob[jj]->Draw("ALP");
    }
    c_chi2->Modified();
    c_chi2->Update();

    condition = //chi2["prob_phi"]                 > 0.05 &&
                //chi2["prob_cosTheta1"]           > 0.05 &&
                //chi2["prob_cosTheta2"]           > 0.05 &&
                //chi2["prob_phi_cosTheta1"]       > 0.05 &&
                //chi2["prob_phi_cosTheta2"]       > 0.05 &&
                //chi2["prob_cosTheta1_cosTheta2"] > 0.05;
                //chi2["prob_phi_cosTheta1_cosTheta2"] > 0.05;
                //chi2["chi2NDF_phi"]                 < 2 &&
                //chi2["chi2NDF_cosTheta1"]           < 2 &&
                //chi2["chi2NDF_cosTheta2"]           < 2 &&
                chi2["chi2NDF_phi_cosTheta1"]       < 2 &&
                chi2["chi2NDF_phi_cosTheta2"]       < 2 &&
                chi2["chi2NDF_cosTheta1_cosTheta2"] < 2;

    if ( !opt_fastSearch ) {
      // Scan starting from the smallest number of parameters and increment by 1
      j++;
      if ( condition ) break;
    } else {
      // Jump by 2x there and back to find the optimum
      if ( condition ) {
        max_j = j;
        if ( min_j == max_j ) break;
        j -= (j+1 - min_j)/2;
      } else {
        min_j = j+1;
        j += (max_j+1 - j)/2;
      }
    }
  }

  printf("Non-zero aklm parameters:\n");
  i = 3;
  for ( Int_t k=0; k<=max_k; k++ ) {
    for ( Int_t l=0; l<=max_l; l++ ) {
      for ( Int_t m=-l; m<=l; m++ ) {
        if ( f3->GetParameter(i) != 0 ) {
          printf(" - a[k=%d][l=%d][m=%d] = %g\n", k, l, m, aklm[k][l][m]);
        }
        i++;
      }
    }
  }

  // Draw the final result of the optimization
  GetChi2NDF(h3, f3, "minimal aklm set", "plot text", smooth3D);

  // Draw the final result with fine binning
  printf("Creating 3D histogram for fine 1D bining\n");
  if ( nBins1D == -1 ) nBins1D = TMath::Min((Int_t)(TMath::Max(5., entries/100.) + 0.1), 100);
  if ( nBins1D != nBins3D ) {
    TH3D *h3_fine = new TH3D("h3_fine", "h3_fine", nBins1D, -TMath::Pi(), TMath::Pi(), nBins1D, -1, 1, nBins1D, -1, 1);
    t_data->Draw(angle_names[2]+":"+angle_names[1]+":"+angle_names[0] + ">>+ " + TString(h3_fine->GetName()), selection, "goff");
    h3_fine->SetXTitle(angle_title[0]);
    h3_fine->SetYTitle(angle_title[1]);
    h3_fine->SetZTitle(angle_title[2]);
    f3->FixParameter(0, entries * 2.*TMath::Pi()/h3_fine->GetNbinsX() * 2./h3_fine->GetNbinsY() * 2./h3_fine->GetNbinsZ());
    if ( opt_checkNorm ) {
      printf("Checking normalization...\n");
      printf("Normalization check: %g %g\n", f3->Integral(-TMath::Pi(), TMath::Pi(), -1, 1, -1, 1), h3->Integral("Width"));
    }
    GetChi2NDF(h3_fine, f3, "minimal aklm set, optimized for 1D projections", "plot text", smooth1D);
  } else {
    printf (" - not needed, same binning defined for 3D and 1D\n");
  }
}
