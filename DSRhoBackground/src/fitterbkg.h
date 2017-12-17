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
#include "RooChebychev.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
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

    void ReadInFile(const char* file_path, const int& num_events = 0);
    void SetPlotDir(const char* output_dir);
    void Fit(RooAbsPdf* pdf, RooDataSet* data);

    RooRealVar thetat_{"thetat", "thetat", 0, constants::pi};
    RooRealVar thetab_{"thetab", "thetab", 0, constants::pi};
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

    RooRealVar* dt_;
    RooFormulaVar* dt_formula_;

    RooRealVar* vrzerr_;
    RooFormulaVar* vrzerr_formula_;

    RooRealVar* vtzerr_;
    RooFormulaVar* vtzerr_formula_;

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

    // Self-cross-feed phit model
    RooRealVar scf_phit_poly_p2_{"scf_phit_poly_p2", "p_{2}", 0.4, -0.1, 2};
    RooRealVar scf_phit_f_{"scf_phit_f", "f_{poly}", 0.25, 0.1, 0.9};
    RooPolynomial scf_phit_poly_{"scf_phit_poly", "scf_phit_poly", phit_, scf_phit_poly_p2_, 2};
    RooRealVar scf_phit_offset_{"scf_phit_offset", "#phi_{t}^{offset}", 0, -0.1, 0.1};
    RooFormulaVar scf_phit_phit_{"scf_phit_phit", "scf_phit_phit", "phit - scf_phit_offset",
                                 RooArgList(phit_, scf_phit_offset_)};
    RooGenericPdf scf_phit_cos_{"scf_phit_cos", "scf_phit_cos", "cos(scf_phit_phit)^2",
                                RooArgList(scf_phit_phit_)};

    // Self-cross-feed thetat model
    RooRealVar scf_thetat_poly_p1_{"scf_thetat_poly_p1", "p_{1}", 0.4, -10, 10};
    RooRealVar scf_thetat_poly_p2_{"scf_thetat_poly_p2", "p_{2}", 0.4, -10, 10};
    RooRealVar scf_thetat_poly_p3_{"scf_thetat_poly_p3", "p_{3}", 0.4, -10, 10};
    RooRealVar scf_thetat_poly_p4_{"scf_thetat_poly_p4", "p_{4}", 0.4, -10, 10};

    RooRealVar scf_thetat_f_{"scf_thetat_f", "#theta_{t}^{w}", 0, -0.1, 0.1};
    RooFormulaVar scf_thetat_thetat_{"scf_thetat_thetat", "scf_thetat_thetat",
                                     "(thetat - 1.5708)*(1+scf_thetat_f) + 1.5708",
                                     RooArgList(thetat_, scf_thetat_f_)};

   public:
    RooAddPdf scf_phit_model_{"scf_phit_model", "scf_phit_model",
                              RooArgList(scf_phit_poly_, scf_phit_cos_), RooArgList(scf_phit_f_)};

    RooGenericPdf scf_thetat_model_{"scf_thetat_model", "scf_thetat_model",
                                    "sin(scf_thetat_thetat)^3", RooArgList(scf_thetat_thetat_)};
};

#endif /* FITTERBKG_H_ */
