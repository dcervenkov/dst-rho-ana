/**
 *  @file    fittershe.h
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

#pragma once

// Standard includes
#include <map>

// ROOT includes
#include "TF3.h"
#include "TH3D.h"
#include "TROOT.h"
#include "TString.h"
#include "TTree.h"

class FitterSHE {
   public:
    FitterSHE();
    ~FitterSHE();

    void AnalyzeDataset(TString dataset = "BsJpsiPhi", Int_t max_entries = -1, Int_t nBins3D = 20,
                        Int_t nBins1D = 50, Int_t smooth3D = 1, Int_t smooth1D = 1,
                        Int_t max_k = 10, Int_t max_l = 10, TString options = "");

    void ReadInFile(std::vector<const char*> file_names, const int& num_events = 0);
    void SetPlotDir(const char* output_dir);
    void SetOutputFile(const char* filename);
    const void LogEnvironmentMetadata();
    const void LogCLIArguments(int argc, char* argv[]);
    const void LogTextFromFile(const char* field_name, const char* filename);
    const void LogFileCRC(const char* field_name, const char* filename);
    const void LogText(const char* field_name, const char* text);

   private:
    static Double_t Ylm(UInt_t l, Int_t m, Double_t cosTheta1, Double_t phi);
    static Double_t Pk(UInt_t k, Double_t cosTheta2);
    static Double_t SumAklmYlmPk(double *x, double *par);

    std::map<TString, Double_t> GetChi2NDF(TH3D *h3_phi_cosTheta1_cosTheta2, TF3 *f3,
                                           TString title = "", TString options = "plot text",
                                           Int_t smoothScale = 100);

    TFile* output_file_ = nullptr;
    TTree* input_data_ = nullptr;

    bool make_plots_ = false;
};
