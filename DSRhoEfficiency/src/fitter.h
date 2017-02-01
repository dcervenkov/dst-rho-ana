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
#include "TCanvas.h"
#include "RooRealVar.h"
#include "RooChebychev.h"
#include "RooGaussian.h"
#include "RooCategory.h"
#include "RooExtendPdf.h"
#include "RooGenericPdf.h"
#include "RooFormulaVar.h"

// Local includes
#include "constants.h"

class Fitter {
public:
	Fitter(const char* evtgen_filepath, const char* gsim_filepath, const char* output_dir);
	virtual ~Fitter();
	void PlotVar(RooRealVar& var);
	void PlotVars2D(RooRealVar& var1, RooRealVar& var2);
	void PlotEfficiency(RooRealVar& var, bool plot_model = true, bool legend_position_top = true, bool legend_position_left = true);
	void PlotEfficiency2D(RooRealVar& var1, RooRealVar& var2);
	void FitEfficiency(RooRealVar& var);

	RooRealVar thetat_ { "thetat", "#theta_{t} [rad]", 0, kPi };
	RooRealVar thetab_ { "thetab", "#theta_{b} [rad]", 0.5, 2.95 };
	RooRealVar phit_ { "phit", "#phi_{t} [rad]", -kPi, kPi };
	RooRealVar dt_ { "dt", "dt", -15, +15 };
	RooCategory dec_type_ { "dec_type", "dec_type" };
	RooRealVar evmcflag_ { "evmcflag", "evmcflag", -100, 100 };

private:
	RooRealVar vrvtxz_ { "vrvtxz", "vrvtxz", 0 };
	RooRealVar vtvtxz_ { "vtvtxz", "vtvtxz", 0 };
	RooFormulaVar dt_formula_ {"dt", "#Deltat [ps]", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)", RooArgSet(vrvtxz_, vtvtxz_)};

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

	RooChebychev thetat_model_ { "thetat_model", "thetat_model", thetat_, RooArgList() };
	RooRealVar n_thetat_model_ { "n_thetat_model", "n_{#theta_{t}}", 0, 100000 };
	RooExtendPdf thetat_model_e_ { "thetat_model_e", "thetat_model_e", thetat_model_, n_thetat_model_ };


	RooRealVar thetab_model_p0_ { "thetab_model_p0", "p_{0}", 0.089, -100, 100};
	RooRealVar thetab_model_p1_ { "thetab_model_p1", "p_{1}", -0.358, -100, 100};
	RooRealVar thetab_model_p2_ { "thetab_model_p2", "p_{2}", -0.024, -100, 100};
	RooChebychev thetab_model_ { "thetab_model", "thetab_model", thetab_,
		RooArgList(thetab_model_p0_, thetab_model_p1_, thetab_model_p2_)};
	RooRealVar n_thetab_model_ { "n_thetab_model", "n_{#theta}_{b}", 0, 100000 };
	RooExtendPdf thetab_model_e_ { "thetab_model_e", "thetab_model_e", thetab_model_, n_thetab_model_ };


	RooRealVar phit_gaus_mean_ { "phit_gaus_mean", "#mu_{G}", -1, 1};
	RooRealVar phit_gaus_sigma_ { "phit_gaus_sigma", "#sigma_{G}", 0, 1.5};
	RooRealVar f_phit_gaus_ { "f_phit_gaus", "f_{G}", -1.0, 0.0 };
	RooGenericPdf phit_model_ { "phit_model", "phit_model", "1 + f_phit_gaus*TMath::Gaus(phit, phit_gaus_mean, phit_gaus_sigma, 1)",
		RooArgSet(phit_, f_phit_gaus_, phit_gaus_mean_, phit_gaus_sigma_)};
	RooRealVar n_phit_model_ { "n_phit_model", "n_{#phi}_{t}", 0, 100000 };
	RooExtendPdf phit_model_e_ { "phit_model_e", "phit_model_e", phit_model_, n_phit_model_ };

	RooAbsPdf* model_ = NULL;

	RooFitResult* result_ = NULL;

	TFile* output_file_;
};

#endif /* FITTER_H_ */
