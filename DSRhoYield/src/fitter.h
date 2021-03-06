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
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooBifurGauss.h"
#include "RooCBShape.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooExtendPdf.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooNovosibirsk.h"
#include "RooPolynomial.h"
#include "RooRealVar.h"
#include "TChain.h"
#include "TString.h"

// Local includes
#include "constants.h"

enum class Components { signal, crossfeed, signal_plus_crossfeed, background, all };

class Fitter {
   public:
    Fitter(const char* output_dir);
    virtual ~Fitter();

    void Setup(Components component);
    void FixShape(Components component);
    void FitTo(TTree* chain);
    void Plot();
    TPaveText* CreateStatBox(double chi2);
    void WriteFitResults();
    void CloseOutput();
    void SPlotFull(TChain* chain);
    void SPlotSC(TChain* chain);
    void SPlotSB(TChain* chain);
    void SPlotCB(TChain* chain);
    double GetCorrelation(TChain* chain, const RooRealVar& var1, const RooRealVar& var2,
                          bool save_plot = false);
    void PrintResults();
    void SetNumCPUs(int numCPUs) { numCPUs_ = numCPUs; };
    void SetMC(bool MC) { MC_ = MC; };

   private:
    int numCPUs_ = 1;
    bool MC_ = false;

    RooArgSet dataset_vars_;

    RooRealVar thetat_{"thetat", "#theta_{t} [rad]", constants::cuts::thetat_low,
                       constants::cuts::thetat_high};
    RooRealVar thetab_{"thetab", "#theta_{b} [rad]", constants::cuts::thetab_low,
                       constants::cuts::thetab_high};
    RooRealVar phit_{"phit", "#phi_{t} [rad]", constants::cuts::phit_low,
                     constants::cuts::phit_high};
    RooRealVar dt_{"dt", "dt", constants::cuts::dt_low, constants::cuts::dt_high};

    RooRealVar de_{"de", "#DeltaE [GeV]", constants::cuts::de_low, constants::cuts::de_high};
    RooRealVar candsel_{"candsel", "candsel", 0, 4};
    RooRealVar evmcflag_{"evmcflag", "evmcflag", -100, 100};
    RooRealVar csbdtg_{"csbdtg", "csbdtg", -1, 1};

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

    RooRealVar nocand_{"nocand", "nocand", 1, 6};

    RooRealVar width_factor_{"width_factor", "F_{width}", 1, 0.5, 1.5};

    // Signal Gaussian
    RooRealVar sig_gaus_mean_{"sig_gaus_mean", "#mu_{Gauss}", -0.1, 0.1};
    RooRealVar sig_gaus_sigma_l_{"sig_gaus_sigma_l", "#sigma_{Gauss}^{l}", 0.001, 0.1};
    RooFormulaVar sig_gaus_sigma_l_w_{"sig_gaus_sigma_l_w", "w#sigma_{Gauss}^{l}",
                                      "width_factor*sig_gaus_sigma_l",
                                      RooArgList(width_factor_, sig_gaus_sigma_l_)};
    RooRealVar sig_gaus_sigma_r_{"sig_gaus_sigma_r", "#sigma_{Gauss}^{r}", 0.001, 0.1};
    RooFormulaVar sig_gaus_sigma_r_w_{"sig_gaus_sigma_r_w", "w#sigma_{Gauss}^{r}",
                                      "width_factor*sig_gaus_sigma_r",
                                      RooArgList(width_factor_, sig_gaus_sigma_r_)};
    RooBifurGauss sig_gaus_{"sig_gaus",     "sig_gaus",          de_,
                            sig_gaus_mean_, sig_gaus_sigma_l_w_, sig_gaus_sigma_r_w_};

    // Signal crystal-ball
    RooRealVar sig_cb_t0_{"sig_cb_t0", "t0_{cb}", -0.1, 0.0};
    RooRealVar sig_cb_sigma_{"sig_cb_sigma", "#sigma_{cb}", 0.03, 0.001, 0.3};
    RooFormulaVar sig_cb_sigma_w_{"sig_cb_sigma_w", "w#sigma_{cb}", "width_factor*sig_cb_sigma",
                                  RooArgList(width_factor_, sig_cb_sigma_)};
    RooRealVar sig_cb_cut_{"sig_cb_cut", "cut_{cb}", 0.85, 0.1, 2};
    RooRealVar sig_cb_power_{"sig_cb_power", "power_{cb}", 3.67, 0.1, 10};
    RooCBShape sig_cb_{"sig_cb",        "sig_cb",    de_,          sig_cb_t0_,
                       sig_cb_sigma_w_, sig_cb_cut_, sig_cb_power_};

    // Signal model
    RooRealVar sig_f_{"sig_f", "f_{cb/Gauss}", 0.1, 0.9};
    RooAddPdf signal_model_{"signal_model", "sig_cb + sig_gaus", sig_cb_, sig_gaus_, sig_f_};

    // Cross-feed Gaussian
    RooRealVar cross_gaus_mean_{"cross_gaus_mean", "#mu_{Gauss}", -0.1, 0.1};
    RooRealVar cross_gaus_sigma_{"cross_gaus_sigma", "#sigma_{Gauss}", 0.01, 1};
    RooGaussian cross_gaus_{"cross_gaus", "cross_gaus", de_, cross_gaus_mean_, cross_gaus_sigma_};

    // Cross-feed exponential
    RooRealVar cross_tau_{"cross_tau", "#tau_{exp}", -8, 10};
    RooExponential cross_exponential_{"cross_exponential", "cross_exponential", de_, cross_tau_};

    // Cross-feed model
    RooRealVar cross_f_{"cross_f", "f_{Gauss/exp}", 0.1, 0.9};
    RooAddPdf cross_model_{"cross_model", "cross_gaus + cross_exponential", cross_gaus_,
                           cross_exponential_, cross_f_};

    // Signal + cross-feed model
    RooRealVar signal_plus_cross_f_{"signal_plus_cross_f", "f_{sig/cf}", 0.1, 0.99};
    RooAddPdf signal_plus_cross_model_{"signal_plus_cross_model", "signal_model + cross_model",
                                       signal_model_, cross_model_, signal_plus_cross_f_};
    RooRealVar n_signal_plus_cross_model_{"n_signal_plus_cross_model", "n_{sig+cf}", 20000, 100,
                                          1000000};
    RooExtendPdf signal_plus_cross_model_e_{"signal_plus_cross_model_e",
                                            "signal_plus_cross_model_e", signal_plus_cross_model_,
                                            n_signal_plus_cross_model_};

    // Background model
    RooRealVar bkg_poly_p1_{"bkg_poly_p1", "p_{1}", -100, 100};
    RooRealVar bkg_poly_p2_{"bkg_poly_p2", "p_{2}", -100, 100};
    RooChebychev bkg_model_{"bkg_model", "bkg_model", de_, RooArgList(bkg_poly_p1_, bkg_poly_p2_)};
    RooRealVar n_bkg_{"n_bkg", "n_{bkg}", 20000, 10, 1000000};
    RooExtendPdf bkg_model_e_{"bkg_model_e", "bkg_model_e", bkg_model_, n_bkg_};

    // Whole model
    RooAddPdf model_{"model", "signal_plus_cross_model_e + bkg_model_e",
                     RooArgList(signal_plus_cross_model_e_, bkg_model_e_)};

    RooAbsPdf* active_model_ = nullptr;
    Components active_component_{Components::all};
    TString active_model_name_ = {"all"};

    RooDataSet* data_set_ = nullptr;
    RooFitResult* fit_result_ = nullptr;

    TString data_cut_;

    TFile* output_file_ = nullptr;

    // Helper variables which are used to calculate final results with their uncertainties
    RooFormulaVar f_cr_{
        "f_cr", "f_cr",
        "(n_signal_plus_cross_model * signal_plus_cross_f) / (n_signal_plus_cross_model + n_bkg)",
        RooArgList(n_signal_plus_cross_model_, signal_plus_cross_f_, n_bkg_)};

    RooFormulaVar f_scf_{"f_scf", "f_scf",
                         "(n_signal_plus_cross_model * (1 - signal_plus_cross_f)) / "
                         "(n_signal_plus_cross_model + n_bkg)",
                         RooArgList(n_signal_plus_cross_model_, signal_plus_cross_f_, n_bkg_)};

    RooFormulaVar f_bkg_{"f_bkg", "f_bkg", "(n_bkg) / (n_signal_plus_cross_model + n_bkg)",
                         RooArgList(n_signal_plus_cross_model_, n_bkg_)};

    RooFormulaVar n_cr_{"n_cr", "n_cr", "n_signal_plus_cross_model * signal_plus_cross_f",
                        RooArgList(n_signal_plus_cross_model_, signal_plus_cross_f_)};

    RooFormulaVar n_scf_{"n_scf", "n_scf", "n_signal_plus_cross_model * (1 - signal_plus_cross_f)",
                         RooArgList(n_signal_plus_cross_model_, signal_plus_cross_f_)};
};

#endif /* FITTER_H_ */
