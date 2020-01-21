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

class FitterBKG {
   public:
    FitterBKG();
    virtual ~FitterBKG();

    void SetNumCPUs(const int& numCPUs) { num_CPUs_ = numCPUs; };
    int GetNumCPUs() { return num_CPUs_; };

    void PlotVar(RooRealVar& var, const RooDataSet* data) const;
    void PlotWithPull(const RooRealVar& var, const RooDataSet&, const RooAbsPdf& pdf,
                      const std::vector<RooAbsPdf*> components = std::vector<RooAbsPdf*>(),
                      const char* title = "") const;
    void PlotKDE(AdaptiveKernelDensity kde) const;

    void ReadInFile(std::vector<const char*> file_names, const int& num_events = 0);
    void SetPlotDir(const char* output_dir);
    void Fit(RooAbsPdf* pdf, RooDataSet* data);
    void CreateHistoPDF(RooDataSet* data, const std::string results_file);
    AdaptiveKernelDensity CreateKDEPDF(RooDataSet* data, const std::string results_file);
    void PrintResultsJSON() const;

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

    RooRealVar* tau_;
    RooRealVar* dm_;

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

    RooCategory* decaytype_;

    RooDataSet* dataset_ = nullptr;
    RooDataSet* dataset_a_ = nullptr;
    RooDataSet* dataset_ab_ = nullptr;
    RooDataSet* dataset_b_ = nullptr;
    RooDataSet* dataset_bb_ = nullptr;

   private:
    TPaveText* CreateStatBox(const double chi2, const int ndof, const bool position_top,
                             const bool position_left) const;
    TH3F* ConvertDensityToHisto(AdaptiveKernelDensity pdf) const;
    TH3F* Create3DHisto(const RooDataSet* dataset) const;

    TH3F* NormalizePDF(const TH3F* pdf, const double low, const double high);
    double Interpolate(const TH3F* histo, int x_org, int y_org, int z_org, int size);

    std::vector<RooRealVar**> conditional_vars_;
    std::vector<RooRealVar**> dataset_vars_;
    RooArgSet conditional_vars_argset_;
    RooArgSet dataset_vars_argset_;
    std::vector <RooRealVar*> model_parameters_;

    int num_CPUs_;

    RooFitResult* result_ = nullptr;

    TFile* output_file_ = nullptr;
    TChain* input_tree = nullptr;
    TTree* data_tree = nullptr;

    // Physics-based background dt model
    RooRealVar bkg_dt_tau_{"bkg_dt_tau", "#tau_{bkg}", 1.5, 0.1, 10};
    RooRealVar bkg_dt_f_delta_{"bkg_dt_f_delta", "f_{d}", 0.1, 0, 1};
    RooRealVar bkg_dt_mu_delta_{"bkg_dt_mu_delta", "#mu_{d}", 0, -1, 1};
    RooRealVar bkg_dt_mu_lifetime_{"bkg_dt_mu_lifetime", "#mu_{l}", 0, -1, 1};
    RooRealVar bkg_dt_f_tail_{"bkg_dt_f_tail", "f_{t}", 0.1, 0, 1};
    RooRealVar bkg_dt_S_main_{"bkg_dt_S_main", "S_{m}", 1, 0, 1000};
    RooRealVar bkg_dt_S_tail_{"bkg_dt_S_tail", "S_{t}", 1, 0, 1000};

    // Background dt model
    RooRealVar bkg_dt_voigt_mu_{"bkg_dt_voigt_mu", "v_{#mu}", -0.303, -1, 1};
    RooRealVar bkg_dt_voigt_sigma_{"bkg_dt_voigt_sigma_", "v_{#sigma}", 2.323, 0, 10};
    RooRealVar bkg_dt_voigt_width_{"bkg_dt_voigt_width_", "v_{w}", 0.851, 0, 10};
    RooVoigtian bkg_dt_voigt_{"bkg_dt_voigt", "bkg_dt_voigt", dt_, bkg_dt_voigt_mu_, bkg_dt_voigt_width_, bkg_dt_voigt_sigma_};
    
    RooRealVar bkg_dt_gaus_mu_{"bkg_dt_gaus_mu", "g_{#mu}", -0.161 -1, 1};
    RooRealVar bkg_dt_gaus_sigma_{"bkg_dt_gaus_sigma_", "g_{#sigma}", 1.096, 0, 10};
    RooGaussian bkg_dt_gaus_{"bkg_dt_gaus", "bkg_dt_gaus", dt_, bkg_dt_gaus_mu_, bkg_dt_gaus_sigma_};

    RooRealVar bkg_dt_f_{"bkg_dt_f", "f_{v/g}", 0.631, 0, 1};

    // Background phit model
    RooRealVar bkg_phit_poly_p2_{"bkg_phit_poly_p2", "p_{2}", 0.856, -0.1, 2};
    RooRealVar bkg_phit_f_{"bkg_phit_f", "f_{poly}", 0.147, 0.1, 0.9};
    RooPolynomial bkg_phit_poly_{"bkg_phit_poly", "bkg_phit_poly", phit_, bkg_phit_poly_p2_, 2};
    RooRealVar bkg_phit_offset_{"bkg_phit_offset", "#phi_{t}^{offset}", 0.056, -0.2, 0.2};
    RooFormulaVar bkg_phit_phit_{"bkg_phit_phit", "bkg_phit_phit", "phit - bkg_phit_offset",
                                 RooArgList(phit_, bkg_phit_offset_)};
    RooGenericPdf bkg_phit_cos_{"bkg_phit_cos", "bkg_phit_cos", "cos(bkg_phit_phit)^2",
                                RooArgList(bkg_phit_phit_)};

    // Background thetat model
    RooRealVar bkg_thetat_p1_{"bkg_thetat_p1", "p_{1}", 0};
    RooRealVar bkg_thetat_p2_{"bkg_thetat_p2", "p_{2}", -1.147, -10, 10};
    RooRealVar bkg_thetat_p3_{"bkg_thetat_p3", "p_{3}", 0};
    RooRealVar bkg_thetat_p4_{"bkg_thetat_p4", "p_{4}", 0.174, -1, 1};
    RooRealVar bkg_thetat_p5_{"bkg_thetat_p5", "p_{5}", 0};
    RooRealVar bkg_thetat_p6_{"bkg_thetat_p6", "p_{6}", -0.029, -1, 1};
    RooArgList bkg_thetat_pars_ {bkg_thetat_p1_, bkg_thetat_p2_, bkg_thetat_p3_, bkg_thetat_p4_, bkg_thetat_p5_, bkg_thetat_p6_};

    // Background thetab model
    RooRealVar bkg_thetab_gaus_mu_{"bkg_thetab_gaus_mu", "#mu", 2.885, 1.5, 3};
    RooRealVar bkg_thetab_gaus_sigma_l_{"bkg_thetab_gaus_sigma_l", "#sigma_{L}", 0.411, 0, 3};
    RooRealVar bkg_thetab_gaus_sigma_r_{"bkg_thetab_gaus_sigma_r", "#sigma_{R}", 0.094, 0, 3};
    RooBifurGauss bkg_thetab_gaus_{
        "bkg_thetab_gaus",  "bkg_thetab_gaus",       thetab_,
        bkg_thetab_gaus_mu_, bkg_thetab_gaus_sigma_l_, bkg_thetab_gaus_sigma_r_};
    RooRealVar bkg_thetab_exp_alpha_{"bkg_thetab_exp_alpha", "#alpha", -4.63, -10, 0.0};
    RooExponential bkg_thetab_exp_{"bkg_thetab_exp", "bkg_thetab_exp", thetab_,
                                   bkg_thetab_exp_alpha_};
    RooRealVar bkg_thetab_f_{"bkg_thetab_f", "f_{exp}", 0.625, 0, 0.8};

   public:
    RooAddPdf bkg_dt_model_{"bkg_dt_model", "bkg_dt_model",
                              RooArgList(bkg_dt_voigt_, bkg_dt_gaus_), RooArgList(bkg_dt_f_)};

    RooAddPdf bkg_phit_model_{"bkg_phit_model", "bkg_phit_model",
                              RooArgList(bkg_phit_poly_, bkg_phit_cos_), RooArgList(bkg_phit_f_)};

    RooChebychev bkg_thetat_model_{"bkg_thetat_model", "bkg_thetat_model", thetat_, bkg_thetat_pars_};

    RooAddPdf bkg_thetab_model_{"bkg_thetab_model", "bkg_thetab_model",
                              RooArgList(bkg_thetab_exp_, bkg_thetab_gaus_), RooArgList(bkg_thetab_f_)};

    RooProdPdf model_{"model", "model", RooArgList(bkg_phit_model_, bkg_thetab_model_, bkg_thetat_model_)};

    DtBKG* bkg_physics_dt_model_;
}

;

#endif /* FITTERBKG_H_ */
