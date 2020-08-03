/**
 *  @file    fitterbkg.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2017-12-15
 *
 *  @brief This class performs the fitting itself as well as plotting
 *
 */

#ifndef FITTERBKG_H_
#define FITTERBKG_H_

// Standard includes
#include <array>

// ROOT includes
#include "RooAddPdf.h"
#include "RooBreitWigner.h"
#include "RooBifurGauss.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "RooVoigtian.h"
#include "RooRealVar.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TPaveText.h"

// Meerkat include
#include "AdaptiveKernelDensity.hh"

// Local includes
#include "constants.h"
#include "dtbkg.h"
#include "nlohmann/json.hpp"

class FitterBKG {
   public:
    FitterBKG();
    virtual ~FitterBKG();

    void SetNumCPUs(const int& numCPUs) { num_CPUs_ = numCPUs; };
    int GetNumCPUs() { return num_CPUs_; };
    void SetNoDeltaPDF();
    void SetNoTailPDF();

    void PlotVar(RooRealVar& var, const RooDataSet* data) const;
    void PlotWithPull(RooRealVar& var, const RooDataSet& data, const RooAbsPdf& pdf) const;
    void PlotKDE(AdaptiveKernelDensity kde) const;

    void ReadInFile(std::vector<const char*> file_names, const int& num_events = 0);
    void SetPlotDir(const char* output_dir);
    void Fit(RooAbsPdf* pdf, RooDataSet* data);
    void CreateHistoPDF(RooDataSet* data, const std::string results_file, const bool randomize);
    AdaptiveKernelDensity CreateKDEPDF(RooDataSet* data, const std::string results_file);
    void PrintResultsJSON() const;

    nlohmann::json FitDt(RooAbsPdf* model, std::string prefix, bool plot, bool randomize);
    nlohmann::json FitAngular(bool plot, bool randomize);

    RooRealVar dt_{"dt", "#Deltat [ps]", constants::cuts::dt_low, constants::cuts::dt_high};
    RooRealVar thetat_{"thetat", "#theta_{t} [rad]", constants::cuts::thetat_low, constants::cuts::thetat_high};
    RooRealVar thetab_{"thetab", "#theta_{b} [rad]", constants::cuts::thetab_low, constants::cuts::thetab_high};
    RooRealVar phit_{"phit", "#phi_{t} [rad]", constants::cuts::phit_low, constants::cuts::phit_high};

    RooRealVar* expno_;
    RooRealVar* expmc_;

    RooRealVar* evmcflag_;
    RooRealVar* brecflav_;
    RooRealVar* btagmcli_;
    RooRealVar* tagqr_;
    RooRealVar* tagwtag_;

    RooRealVar* benergy_;
    RooRealVar* mbc_;
    RooRealVar* de_;
    RooRealVar* csbdtg_;

    RooRealVar* shcosthb_;

    RooRealVar* vrusable_;
    RooRealVar* vrvtxz_;
    RooRealVar* vrerr6_;
    RooRealVar* vrchi2_;
    //  RooRealVar* vreffxi_;
    RooRealVar* vrndf_;
    //  RooRealVar* vreffndf_;
    RooRealVar* vrntrk_;

    RooRealVar* vtusable_;
    RooRealVar* vtvtxz_;
    RooRealVar* vterr6_;
    RooRealVar* vtchi2_;
    RooRealVar* vtndf_;
    RooRealVar* vtntrk_;
    RooRealVar* vtistagl_;

    RooRealVar* nocand_;

    RooCategory* decaytype_;

    RooDataSet* dataset_ = nullptr;
    RooDataSet* dataset_FB_ = nullptr;
    RooDataSet* dataset_FA_ = nullptr;
    RooDataSet* dataset_SB_ = nullptr;
    RooDataSet* dataset_SA_ = nullptr;

   private:
    TH3F* ConvertDensityToHisto(AdaptiveKernelDensity pdf) const;
    TH3F* Create3DHisto(const RooDataSet* dataset) const;

    TH3F* NormalizePDF(const TH3F* pdf, const double low, const double high);
    double Interpolate(const TH3F* histo, int x_org, int y_org, int z_org, int size);

    std::vector<RooRealVar**> conditional_vars_;
    std::vector<RooRealVar**> dataset_vars_;
    RooArgSet conditional_vars_argset_;
    RooArgSet dataset_vars_argset_;

    int num_CPUs_;

    RooFitResult* result_ = nullptr;

    TFile* output_file_ = nullptr;
    TChain* input_tree = nullptr;
    TTree* data_tree = nullptr;

    // Physics-based background dt model
    RooRealVar dt_tau_{"dt_tau", "#tau", 1.5, 0.1, 10};
    RooRealVar dt_f_delta_{"dt_f_delta", "f_{d}", 0.1, 0, 1};
    RooRealVar dt_mu_delta_{"dt_mu_delta", "#mu_{d}", 0, -1, 1};
    RooRealVar dt_mu_lifetime_{"dt_mu_lifetime", "#mu_{l}", 0, -1, 1};
    RooRealVar dt_f_tail_{"dt_f_tail", "f_{t}", 0.1, 0, 1};
    RooRealVar dt_S_main_{"dt_S_main", "S_{m}", 1, 0, 1000};
    RooRealVar dt_S_tail_{"dt_S_tail", "S_{t}", 1, 0, 1000};

    // Background dt model
    RooRealVar dt_voigt_mu_{"dt_voigt_mu", "v_{#mu}", -0.303, -1, 1};
    RooRealVar dt_voigt_sigma_{"dt_voigt_sigma", "v_{#sigma}", 2.323, 0, 10};
    RooRealVar dt_voigt_width_{"dt_voigt_width", "v_{w}", 0.851, 0, 10};
    RooVoigtian dt_voigt_{"dt_voigt", "dt_voigt", dt_, dt_voigt_mu_, dt_voigt_width_, dt_voigt_sigma_};

    RooRealVar dt_gaus_mu_{"dt_gaus_mu", "g_{#mu}", -0.161 -1, 1};
    RooRealVar dt_gaus_sigma_{"dt_gaus_sigma", "g_{#sigma}", 1.096, 0, 10};
    RooGaussian dt_gaus_{"dt_gaus", "dt_gaus", dt_, dt_gaus_mu_, dt_gaus_sigma_};

    RooRealVar dt_f_{"dt_f", "f_{v/g}", 0.631, 0, 1};

    // Background phit model
    RooRealVar phit_poly_p2_{"phit_poly_p2", "p_{2}", 0.856, -0.1, 2};
    RooRealVar phit_f_{"phit_f", "f_{poly}", 0.147, 0.1, 0.9};
    RooPolynomial phit_poly_{"phit_poly", "phit_poly", phit_, phit_poly_p2_, 2};
    RooRealVar phit_offset_{"phit_offset", "#phi_{t}^{offset}", 0.056, -0.2, 0.2};
    RooFormulaVar phit_phit_{"phit_phit", "phit_phit", "phit - phit_offset",
                                 RooArgList(phit_, phit_offset_)};
    RooGenericPdf phit_cos_{"phit_cos", "phit_cos", "cos(phit_phit)^2",
                                RooArgList(phit_phit_)};

    // Background thetat model
    RooRealVar thetat_p1_{"thetat_p1", "p_{1}", 0};
    RooRealVar thetat_p2_{"thetat_p2", "p_{2}", -1.147, -10, 10};
    RooRealVar thetat_p3_{"thetat_p3", "p_{3}", 0};
    RooRealVar thetat_p4_{"thetat_p4", "p_{4}", 0.174, -1, 1};
    RooRealVar thetat_p5_{"thetat_p5", "p_{5}", 0};
    RooRealVar thetat_p6_{"thetat_p6", "p_{6}", -0.029, -1, 1};
    RooArgList thetat_pars_ {thetat_p1_, thetat_p2_, thetat_p3_, thetat_p4_, thetat_p5_, thetat_p6_};

    // Background thetab model
    RooRealVar thetab_gaus_mu_{"thetab_gaus_mu", "#mu", 2.885, 1.5, 3};
    RooRealVar thetab_gaus_sigma_l_{"thetab_gaus_sigma_l", "#sigma_{L}", 0.411, 0, 3};
    RooRealVar thetab_gaus_sigma_r_{"thetab_gaus_sigma_r", "#sigma_{R}", 0.094, 0, 3};
    RooBifurGauss thetab_gaus_{
        "thetab_gaus",  "thetab_gaus",       thetab_,
        thetab_gaus_mu_, thetab_gaus_sigma_l_, thetab_gaus_sigma_r_};
    RooRealVar thetab_exp_alpha_{"thetab_exp_alpha", "#alpha", -4.63, -10, 0.0};
    RooExponential thetab_exp_{"thetab_exp", "thetab_exp", thetab_,
                                   thetab_exp_alpha_};
    RooRealVar thetab_f_{"thetab_f", "f_{exp}", 0.625, 0, 0.8};

   public:
    RooAddPdf dt_model_{"dt_model", "dt_model",
                              RooArgList(dt_voigt_, dt_gaus_), RooArgList(dt_f_)};

    RooAddPdf phit_model_{"phit_model", "phit_model",
                              RooArgList(phit_poly_, phit_cos_), RooArgList(phit_f_)};

    RooChebychev thetat_model_{"thetat_model", "thetat_model", thetat_, thetat_pars_};

    RooAddPdf thetab_model_{"thetab_model", "thetab_model",
                              RooArgList(thetab_exp_, thetab_gaus_), RooArgList(thetab_f_)};

    RooProdPdf model_{"model", "model", RooArgList(phit_model_, thetab_model_, thetat_model_)};

    DtBKG* physics_dt_model_;
}

;

#endif /* FITTERBKG_H_ */
