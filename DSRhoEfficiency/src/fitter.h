/**
 *  @file    fitter.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2015-07-31
 *
 *  @brief This class performs the yield fitting itself as well as plotting
 *
 */

#ifndef FITTER_H_
#define FITTER_H_

// ROOT includes
#include "RooAddPdf.h"
#include "RooCategory.h"
#include "RooChebychev.h"
#include "RooExtendPdf.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooRealVar.h"
#include "TCanvas.h"

// Meerkat includes
#include "AdaptiveKernelDensity.hh"
#include "BinnedKernelDensity.hh"
#include "CombinedPhaseSpace.hh"
#include "KernelDensity.hh"

// Local includes
#include "constants.h"

class Fitter {
   public:
    Fitter(const char* evtgen_filepath, const char* gsim_filepath, const char* output_dir);
    virtual ~Fitter();
    void PlotVar(const RooRealVar& var) const;
    void PlotVar(const RooRealVar& var, const RooDataHist& data1, const RooDataHist& data2,
                 const bool draw_pull, const bool draw_residual) const;
    void PlotVars2D(const RooRealVar& var1, const RooRealVar& var2) const;
    void PlotEfficiency(RooRealVar& var, bool plot_model = true, bool legend_position_top = true,
                        bool legend_position_left = true);
    void PlotEfficiency2D(RooRealVar& var1, RooRealVar& var2);
    void FitEfficiency(RooRealVar& var);
    void SetEfficiencyModel(const int model_num);

    void ProcessKDEEfficiency(const char* efficiency_file, const std::array<double, 6> bin_kde_pars,
                              const std::array<double, 6> ada_kde_pars, const double mirror_margin = 0);
    void ProcessKDEEfficiency2(const char* efficiency_file, const std::array<double, 6> bin_kde_pars,
                               const std::array<double, 6> ada_kde_pars, const double mirror_margin = 0);
    void ProcessNormalizedEfficiency(const char* efficiency_file);
    
    RooRealVar thetat_{"thetat", "#theta_{t} [rad]", constants::cuts::thetat_low, constants::cuts::thetat_high};
    RooRealVar thetab_{"thetab", "#theta_{b} [rad]", constants::cuts::thetab_low, constants::cuts::thetab_high};
    RooRealVar phit_{"phit", "#phi_{t} [rad]", constants::cuts::phit_low, constants::cuts::phit_high};
    RooRealVar dt_{"dt", "dt", constants::cuts::dt_low, constants::cuts::dt_high};
    RooCategory dec_type_{"dec_type", "dec_type"};
    RooRealVar evmcflag_{"evmcflag", "evmcflag", -100, 100};

   private:
    TH3F* GetBinned3DEfficiency();
    TTree* Histogram2TTree(TH3F* histo);
    TH3F* Create3DHisto(const RooDataSet* dataset, const char* name = nullptr) const;
    TH3F* ConvertDensityToHisto(AdaptiveKernelDensity pdf) const;
    TH3F* ConvertDensityToHisto(BinnedKernelDensity pdf) const;
    void MirrorDataAtEdges(RooDataSet* data);
    double* GetMirrorVals(double vals[], const int var_num);
    double* GetMirrorVals(double vals[], const int var_num1, const int var_num2);
    double* GetMirrorVals(double vals[], const int var_num1, const int var_num2,
                          const int var_num3);
    int CloseToEdge(double vals[], const int var_num);
    void EnlargeVarRanges(const double margin);
    TH3F* ConvertTransHisto(TH3F* trans_histo);
	bool CanUseInterpolation(const TH3F* histo, const double& phit, const double& transtht, const double& transthb) const;
    TTree* DoubleTreeToFloatTree(TTree* double_tree) const;
    TH3F* GetKDEHisto(RooDataSet* dataset, const std::array<double, 6> bin_kde_pars) const;
    TH3F* NormalizePDF(const TH3F* pdf, const double low, const double high);
    double Interpolate(const TH3F* histo, int x_org, int y_org, int z_org, int size);

    RooRealVar vrusable_{"vrusable", "vrusable", 0};
    RooRealVar vrvtxz_{"vrvtxz", "vrvtxz", 0};
    RooRealVar vrerr6_{"vrerr6", "vrerr6", 0};
    RooRealVar vrchi2_{"vrchi2", "vrchi2", 0};
    RooRealVar vrndf_{"vrndf", "vrndf", 0};
    RooRealVar vrntrk_{"vrntrk", "vrntrk", 0};

    RooRealVar vtusable_{"vtusable", "vtusable", 0};
    RooRealVar vtvtxz_{"vtvtxz", "vtvtxz", 0};
    RooRealVar vterr6_{"vterr6", "vterr6", 0};
    RooRealVar vtchi2_{"vtchi2", "vtchi2", 0};
    RooRealVar vtndf_{"vtndf", "vtndf", 0};
    RooRealVar vtntrk_{"vtntrk", "vtntrk", 0};

    RooFormulaVar dt_formula_{"dt", "#Deltat [ps]", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)",
                              RooArgSet(vrvtxz_, vtvtxz_)};

    TPaveText* CreateStatBox(double chi2, bool position_top = true, bool position_left = true);

    RooDataSet* evtgen_dataset_ = NULL;
    RooDataSet* gsim_dataset_ = NULL;

    RooDataHist* efficiency_thetat_ = NULL;
    RooDataHist* efficiency_thetab_ = NULL;
    RooDataHist* efficiency_phit_ = NULL;

    RooDataHist** efficiency_ = NULL;

    TCanvas* canvas_var_ = NULL;
    RooPlot* plot_var_ = NULL;

    TCanvas* canvas_eff_ = NULL;
    RooPlot* plot_eff_ = NULL;

    RooExtendPdf* thetat_model_e_;
    RooExtendPdf* thetab_model_e_;
    RooExtendPdf* phit_model_e_;

    RooRealVar* vars_[3] = {&thetat_, &thetab_, &phit_};
    RooRealVar orig_vars_[3] = {thetat_, thetab_, phit_};

    // Model 1
    RooRealVar model1_thetat_p0_{"model1_thetat_p0", "p_{0}", -0.056, -10, +10};
    RooRealVar model1_thetat_p1_{"model1_thetat_p1", "p_{1}", 0.119, -10, +10};
    RooRealVar model1_thetat_p2_{"model1_thetat_p2", "p_{2}", -0.038, -10, +10};
    RooRealVar model1_thetat_p3_{"model1_thetat_p3", "p_{3}", 0.025, -10, +10};
    RooChebychev model1_thetat_{
        "thetat_model", "thetat_model", thetat_,
        RooArgList(model1_thetat_p0_, model1_thetat_p1_, model1_thetat_p2_, model1_thetat_p3_)};
    RooRealVar n_thetat_model1_{"n_thetat_model1", "n_{#theta_{t}}", 0, 100000};
    RooExtendPdf thetat_model1_e_{"thetat_model1_e", "thetat_model1_e", model1_thetat_,
                                  n_thetat_model1_};

    RooRealVar model1_thetab_p0_{"model1_thetab_p0", "p_{0}", 0.089, -10, +10};
    RooRealVar model1_thetab_p1_{"model1_thetab_p1", "p_{1}", -0.358, -10, +10};
    RooRealVar model1_thetab_p2_{"model1_thetab_p2", "p_{2}", -0.024, -10, +10};
    RooChebychev model1_thetab_{
        "model1_thetab", "model1_thetab", thetab_,
        RooArgList(model1_thetab_p0_, model1_thetab_p1_, model1_thetab_p2_)};
    RooRealVar n_thetab_model1_{"n_thetab_model1", "n_{#theta_{b}}", 0, 100000};
    RooExtendPdf thetab_model1_e_{"thetab_model1_e", "thetab_model1_e", model1_thetab_,
                                  n_thetab_model1_};

    RooRealVar model1_phit_p0_{"model1_phit_p0", "p_{0}", 0.020, -10, +10};
    RooRealVar model1_phit_p1_{"model1_phit_p1", "p_{1}", 0.045, -10, +10};
    RooRealVar model1_phit_p2_{"model1_phit_p2", "p_{2}", 0.003, -10, +10};
    RooRealVar model1_phit_p3_{"model1_phit_p3", "p_{3}", -0.060, -10, +10};
    RooRealVar model1_phit_p4_{"model1_phit_p4", "p_{4}", -0.043, -10, +10};
    RooRealVar model1_phit_p5_{"model1_phit_p5", "p_{5}", 0.058, -10, +10};
    RooChebychev model1_phit_{"phit_model", "phit_model", phit_,
                              RooArgList(model1_phit_p0_, model1_phit_p1_, model1_phit_p2_,
                                         model1_phit_p3_, model1_phit_p4_, model1_phit_p5_)};
    RooRealVar n_phit_model1_{"n_phit_model1", "n_{#phit_{b}}", 0, 100000};
    RooExtendPdf phit_model1_e_{"phit_model1_e", "phit_model1_e", model1_phit_, n_phit_model1_};

    // Model 2
    RooPolynomial model2_thetat_{"thetat_model", "thetat_model", thetat_, RooArgList()};
    RooRealVar n_thetat_model2_{"n_thetat_model2", "n_{#theta_{t}}", 0, 100000};
    RooExtendPdf thetat_model2_e_{"thetat_model2_e", "thetat_model2_e", model2_thetat_,
                                  n_thetat_model2_};

    RooRealVar model2_thetab_p0_{"model2_thetab_p0", "p_{0}", 0.089, -10, +10};
    RooRealVar model2_thetab_p1_{"model2_thetab_p1", "p_{1}", -0.358, -10, +10};
    RooRealVar model2_thetab_p2_{"model2_thetab_p2", "p_{2}", -0.024, -10, +10};
    RooChebychev model2_thetab_{
        "model2_thetab", "model2_thetab", thetab_,
        RooArgList(model2_thetab_p0_, model2_thetab_p1_, model2_thetab_p2_)};
    RooRealVar n_thetab_model2_{"n_thetab_model2", "n_{#theta_{b}}", 0, 100000};
    RooExtendPdf thetab_model2_e_{"thetab_model2_e", "thetab_model2_e", model2_thetab_,
                                  n_thetab_model2_};

    RooRealVar model2_phit_p1_{"model2_phit_p1", "p_{1}", 0.007, -10, +10};
    RooRealVar model2_phit_p2_{"model2_phit_p2", "p_{2}", 0.088, -10, +10};
    RooRealVar model2_phit_gaus1_mean_{"model2_phit_gaus1_mean", "#mu_{1}", -1.467, -10, +10};
    RooRealVar model2_phit_gaus1_sigma_{"model2_phit_gaus1_sigma", "#sigma_{1}", 0.429, -10, +10};
    RooRealVar model2_phit_gaus2_mean_{"model2_phit_gaus2_mean", "#mu_{2}", 1.851, -10, +10};
    RooRealVar model2_phit_gaus2_sigma_{"model2_phit_gaus2_sigma", "#sigma_{2}", 0.368, -10, +10};
    RooChebychev model2_phit_cheb_{"model2_phit_cheb", "model2_phit_cheb", phit_,
                                   RooArgList(model2_phit_p1_, model2_phit_p2_)};
    RooGaussian model2_phit_gaus1_{"model2_phit_gaus1", "model2_phit_gaus1", phit_,
                                   model2_phit_gaus1_mean_, model2_phit_gaus1_sigma_};
    RooGaussian model2_phit_gaus2_{"model2_phit_gaus2", "model2_phit_gaus2", phit_,
                                   model2_phit_gaus2_mean_, model2_phit_gaus2_sigma_};
    RooRealVar model2_f_phit_gaus1_{"model2_f_phit_gaus1", "f_{G_{1}}", 0.039, -10, +10};
    RooRealVar model2_f_phit_gaus2_{"model2_f_phit_gaus2", "f_{G_{2}}", 0.038, -10, +10};
    RooAddPdf model2_phit_{"model2_phit", "model2_phit",
                           RooArgList(model2_phit_gaus1_, model2_phit_gaus2_, model2_phit_cheb_),
                           RooArgList(model2_f_phit_gaus1_, model2_f_phit_gaus2_)};
    RooRealVar n_phit_model2_{"n_phit_model2", "n_{#phit_{b}}", 0, 100000};
    RooExtendPdf phit_model2_e_{"phit_model2_e", "phit_model2_e", model2_phit_, n_phit_model2_};

    // Model 3
    RooChebychev thetat_model3_{"thetat_model3", "thetat_model3", thetat_, RooArgList()};
    RooRealVar n_thetat_model3_{"n_thetat_model3", "n_{#theta_{t}}", 0, 100000};
    RooExtendPdf thetat_model3_e_{"thetat_model3_e", "thetat_model3_e", thetat_model3_,
                                  n_thetat_model3_};

    RooRealVar thetab_model3_p0_{"thetab_model3_p0", "p_{0}", 0.089, -10, 10};
    RooRealVar thetab_model3_p1_{"thetab_model3_p1", "p_{1}", -0.358, -10, 10};
    RooRealVar thetab_model3_p2_{"thetab_model3_p2", "p_{2}", -0.024, -10, 10};
    RooChebychev thetab_model3_{
        "thetab_model3", "thetab_model3", thetab_,
        RooArgList(thetab_model3_p0_, thetab_model3_p1_, thetab_model3_p2_)};
    RooRealVar n_thetab_model3_{"n_thetab_model3", "n_{#theta}_{b}", 0, 100000};
    RooExtendPdf thetab_model3_e_{"thetab_model3_e", "thetab_model3_e", thetab_model3_,
                                  n_thetab_model3_};

    RooRealVar phit_gaus_mean_{"phit_gaus_mean", "#mu_{G}", -1, 1};
    RooRealVar phit_gaus_sigma_{"phit_gaus_sigma", "#sigma_{G}", 0, 1.5};
    RooRealVar f_phit_gaus_{"f_phit_gaus", "f_{G}", -1.0, 0.0};
    RooGenericPdf phit_model3_{
        "phit_model3", "phit_model3",
        "1 + f_phit_gaus*TMath::Gaus(phit, phit_gaus_mean, phit_gaus_sigma, 1)",
        RooArgSet(phit_, f_phit_gaus_, phit_gaus_mean_, phit_gaus_sigma_)};
    RooRealVar n_phit_model3_{"n_phit_model3", "n_{#phi}_{t}", 0, 100000};
    RooExtendPdf phit_model3_e_{"phit_model3_e", "phit_model3_e", phit_model3_, n_phit_model3_};

    // Model 4
    RooRealVar model4_thetat_p0_{"model4_thetat_p0", "p_{0}", -0.056, -10, +10};
    RooRealVar model4_thetat_p1_{"model4_thetat_p1", "p_{1}", 0};
    RooRealVar model4_thetat_p2_{"model4_thetat_p2", "p_{2}", -0.038, -10, +10};
    RooFormulaVar thetat_plus_pi2{"thetat_plus_pi2", "thetat + pi/2", "thetat - 1.571",
                                  RooArgSet(thetat_)};
    RooPolynomial thetat_model4_{"thetat_model4", "thetat_model4", thetat_plus_pi2,
                                 RooArgList{model4_thetat_p1_, model4_thetat_p2_}};
    RooRealVar n_thetat_model4_{"n_thetat_model4", "n_{#theta_{t}}", 0, 100000};
    RooExtendPdf thetat_model4_e_{"thetat_model4_e", "thetat_model4_e", thetat_model4_,
                                  n_thetat_model4_};

    RooAbsPdf* model_ = NULL;

    RooFitResult* result_ = NULL;

    TFile* output_file_;
};

#endif /* FITTER_H_ */
