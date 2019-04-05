/**
 *  @file    fit_spherical_harmonics.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2018-06-05
 *
 *  @brief Main header
 *
 */

#ifndef FIT_SPHERICAL_HARMONICS_H_
#define FIT_SPHERICAL_HARMONICS_H_

// Standard includes
#include <map>

// ROOT includes
#include "TF3.h"
#include "TH3D.h"
#include "TROOT.h"
#include "TString.h"

Double_t Ylm(UInt_t l, Int_t m, Double_t cosTheta1, Double_t phi);
Double_t Pk(UInt_t k, Double_t cosTheta2);
Double_t SumAklmYlmPk(double *x, double *par);

std::map<TString, Double_t> GetChi2NDF(TH3D *h3_phi_cosTheta1_cosTheta2, TF3 *f3,
                                       TString title = "", TString options = "plot text",
                                       Int_t smoothScale = 100);

void AnalyzeDataset(TString dataset = "BsJpsiPhi", Int_t max_entries = -1, Int_t nBins3D = 20,
                    Int_t nBins1D = 50, Int_t smooth3D = 1, Int_t smooth1D = 1, Int_t max_k = 10,
                    Int_t max_l = 10, TString options = "");

#endif /* FIT_SPHERICAL_HARMONICS_H_ */