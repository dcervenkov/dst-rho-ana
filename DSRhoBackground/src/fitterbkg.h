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
#include "RooVoigtian.h"
#include "RooRealVar.h"
#include "TCanvas.h"
#include "TPaveText.h"

// Local includes
#include "constants.h"

class FitterBKG {
   public:
    FitterBKG();
    virtual ~FitterBKG();

    void PlotVar(RooRealVar& var, const RooDataSet* data) const;
    void PlotWithPull(const RooRealVar& var, const RooDataSet*, const RooAbsPdf* pdf,
                      const char* title = "") const;

    void ReadInFile(std::vector<const char*> file_names, const int& num_events = 0);
    void SetPlotDir(const char* output_dir);
    void Fit(RooAbsPdf* pdf, RooDataSet* data);

    RooRealVar dt_{"dt", "dt", constants::fit_range_dt_low, constants::fit_range_dt_high};
    RooRealVar thetat_{"thetat", "thetat", 0, constants::pi};
    RooRealVar thetab_{"thetab", "thetab", constants::cuts::thetab_low, constants::cuts::thetab_high};
    RooRealVar phit_{"phit", "phit", -constants::pi, constants::pi};

    RooRealVar* expno_;
    RooRealVar* expmc_;

    RooRealVar* evmcflag_;
    RooRealVar* brecflav_;
    RooRealVar* btagmcli_;
    RooRealVar* tagqr_;
    RooRealVar* tagwtag_;

    RooRealVar* benergy_;
    RooRealVar* mbc_;

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

    RooDataSet* dataset_ = NULL;
    RooDataSet* dataset_a_ = NULL;
    RooDataSet* dataset_ab_ = NULL;
    RooDataSet* dataset_b_ = NULL;
    RooDataSet* dataset_bb_ = NULL;

   private:
    TPaveText* CreateStatBox(const double chi2, const bool position_top,
                             const bool position_left) const;
    TString GetCommonCutsString() const;

    std::vector<RooRealVar**> conditional_vars_;
    std::vector<RooRealVar**> dataset_vars_;
    RooArgSet conditional_vars_argset_;
    RooArgSet dataset_vars_argset_;

    RooFitResult* result_ = NULL;

    TFile* output_file_ = NULL;

    // Self-cross-feed dt model
    RooRealVar scf_dt_voigt_mu_{"scf_dt_voigt_mu", "v_{#mu}", -0.303, -1, 1};
    RooRealVar scf_dt_voigt_sigma_{"scf_dt_voigt_sigma_", "v_{#sigma}", 2.323, 0, 10};
    RooRealVar scf_dt_voigt_width_{"scf_dt_voigt_width_", "v_{w}", 0.851, 0, 10};
    RooVoigtian scf_dt_voigt_{ "scf_dt_voigt",  "scf_dt_voigt", dt_, scf_dt_voigt_mu_, scf_dt_voigt_width_, scf_dt_voigt_sigma_};
    
    RooRealVar scf_dt_gaus_mu_{"scf_dt_gaus_mu", "g_{#mu}", -0.161 -1, 1};
    RooRealVar scf_dt_gaus_sigma_{"scf_dt_gaus_sigma_", "g_{#sigma}", 1.096, 0, 10};
    RooGaussian scf_dt_gaus_{ "scf_dt_gaus",  "scf_dt_gaus", dt_, scf_dt_gaus_mu_, scf_dt_gaus_sigma_};

    RooRealVar scf_dt_f_{"scf_dt_f", "f_{v/g}", 0.631, 0, 1};

    // Self-cross-feed phit model
    RooRealVar scf_phit_poly_p2_{"scf_phit_poly_p2", "p_{2}", 0.856, -0.1, 2};
    RooRealVar scf_phit_f_{"scf_phit_f", "f_{poly}", 0.147, 0.1, 0.9};
    RooPolynomial scf_phit_poly_{"scf_phit_poly", "scf_phit_poly", phit_, scf_phit_poly_p2_, 2};
    RooRealVar scf_phit_offset_{"scf_phit_offset", "#phi_{t}^{offset}", 0.056, -0.1, 0.1};
    RooFormulaVar scf_phit_phit_{"scf_phit_phit", "scf_phit_phit", "phit - scf_phit_offset",
                                 RooArgList(phit_, scf_phit_offset_)};
    RooGenericPdf scf_phit_cos_{"scf_phit_cos", "scf_phit_cos", "cos(scf_phit_phit)^2",
                                RooArgList(scf_phit_phit_)};

    // Self-cross-feed thetat model
    RooRealVar scf_thetat_f_{"scf_thetat_f", "#theta_{t}^{w}", -0.051, -0.5, 0.5};
    RooFormulaVar scf_thetat_thetat_{"scf_thetat_thetat", "scf_thetat_thetat",
                                     "(thetat - 1.5708)*(1+scf_thetat_f) + 1.5708",
                                     RooArgList(thetat_, scf_thetat_f_)};

    // Self-cross-feed thetab model
    RooRealVar scf_thetab_gaus_mu_{"scf_thetab_gaus_mu", "#mu", 2.885, 1.5, 3};
    RooRealVar scf_thetab_gaus_sigma_l_{"scf_thetab_gaus_sigma_l", "#sigma_{L}", 0.411, 0, 3};
    RooRealVar scf_thetab_gaus_sigma_r_{"scf_thetab_gaus_sigma_r", "#sigma_{R}", 0.094, 0, 3};
    RooBifurGauss scf_thetab_gaus_{
        "scf_thetab_gaus",  "scf_thetab_gaus",       thetab_,
        scf_thetab_gaus_mu_, scf_thetab_gaus_sigma_l_, scf_thetab_gaus_sigma_r_};
    RooRealVar scf_thetab_exp_alpha_{"scf_thetab_exp_alpha", "#alpha", -4.63, -10, 0.0};
    RooExponential scf_thetab_exp_{"scf_thetab_exp", "scf_thetab_exp", thetab_,
                                   scf_thetab_exp_alpha_};
    RooRealVar scf_thetab_f_{"scf_thetab_f", "f_{exp}", 0.625, 0, 1};

   public:
    // RooAddPdf scf_dt_model_{"scf_dt_model", "scf_dt_model",
    //                           RooArgList(scf_dt_voigt_, scf_dt_gaus_), RooArgList(scf_dt_f_)};
    RooAddPdf scf_dt_model_{"scf_dt_model", "scf_dt_model",
                              RooArgList(scf_dt_voigt_, scf_dt_gaus_), RooArgList(scf_dt_f_)};

    RooAddPdf scf_phit_model_{"scf_phit_model", "scf_phit_model",
                              RooArgList(scf_phit_poly_, scf_phit_cos_), RooArgList(scf_phit_f_)};

    RooGenericPdf scf_thetat_model_{"scf_thetat_model", "scf_thetat_model",
                                    "sin(scf_thetat_thetat)^3", RooArgList(scf_thetat_thetat_)};

    RooAddPdf scf_thetab_model_{"scf_thetab_model", "scf_thetab_model",
                              RooArgList(scf_thetab_exp_, scf_thetab_gaus_), RooArgList(scf_thetab_f_)};
}

;

#endif /* FITTERBKG_H_ */
