/**
 *  @file    efficiency.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2016-01-22
 *
 *  @brief Class that calculates efficiency (detector acceptance) from a model for a set of parameters
 *
 */

#pragma once

// ROOT includes
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "TH3F.h"

// Meerkat includes
#include "BinnedDensity.hh"
#include "OneDimPhaseSpace.hh"
#include "CombinedPhaseSpace.hh"

// Local includes
#include "constants.h"

class Efficiency {
public:
	Efficiency();
	Efficiency(const char* filename);
	virtual ~Efficiency();
	double GetEfficiency(double thetat, double thetab, double phit, int efficiency_model) const;
	double EfficiencyInterface(double* x, double* p) const;
	void ReadInFile(const char* filename);

protected:
	// void RescaleVars(double& thetat, double& thetab, double& phit, const double margin) const;
	RooRealVar* thetat_ = new RooRealVar{"thetat", "#theta_{t}", constants::cuts::thetat_low, constants::cuts::thetat_high };
	RooRealVar* thetab_ = new RooRealVar{"thetab", "#theta_{b}", constants::cuts::thetab_low, constants::cuts::thetab_high };
	RooRealVar* phit_ = new RooRealVar{ "phit", "#phi_{t}", constants::cuts::phit_low, constants::cuts::phit_high };

	// Model 1
	RooRealVar model1_thetat_p0_ { "model1_thetat_p0", "p_{0}", 0.004249 };
	RooRealVar model1_thetat_p1_ { "model1_thetat_p1", "p_{1}", 0.138054 };
	RooRealVar model1_thetat_p2_ { "model1_thetat_p2", "p_{2}", 0.000536 };
	RooRealVar model1_thetat_p3_ { "model1_thetat_p3", "p_{3}", 0.031963 };
	RooChebychev model1_thetat_ { "thetat_model", "thetat_model", *thetat_,
			RooArgList(model1_thetat_p0_, model1_thetat_p1_, model1_thetat_p2_, model1_thetat_p3_) };

	RooRealVar model1_thetab_p0_ { "model1_thetab_p0", "p_{0}", 0.089 };
	RooRealVar model1_thetab_p1_ { "model1_thetab_p1", "p_{1}", -0.358 };
	RooRealVar model1_thetab_p2_ { "model1_thetab_p2", "p_{2}", -0.024 };
	RooChebychev model1_thetab_ {"model1_thetab", "model1_thetab", *thetab_,
			RooArgList(model1_thetab_p0_, model1_thetab_p1_, model1_thetab_p2_)};

	RooRealVar model1_phit_p0_ { "model1_phit_p0", "p_{0}", 0.020};
	RooRealVar model1_phit_p1_ { "model1_phit_p1", "p_{1}", 0.045};
	RooRealVar model1_phit_p2_ { "model1_phit_p2", "p_{2}", 0.003};
	RooRealVar model1_phit_p3_ { "model1_phit_p3", "p_{3}", -0.060};
	RooRealVar model1_phit_p4_ { "model1_phit_p4", "p_{4}", -0.043};
	RooRealVar model1_phit_p5_ { "model1_phit_p5", "p_{5}", 0.058};
	RooChebychev model1_phit_ { "phit_model", "phit_model", *phit_,
		RooArgList(model1_phit_p0_, model1_phit_p1_, model1_phit_p2_, model1_phit_p3_,
				model1_phit_p4_, model1_phit_p5_) };

	double GetModel1Efficiency() const { return model1_thetat_.getVal() * model1_thetab_.getVal() * model1_phit_.getVal(); }

	// Model 2
	RooPolynomial model2_thetat_ {"thetat_model", "thetat_model", *thetat_, RooArgList() };

	RooRealVar model2_thetab_p0_ { "model2_thetab_p0", "p_{0}", 0.075 };
	RooRealVar model2_thetab_p1_ { "model2_thetab_p1", "p_{1}", -0.369 };
	RooRealVar model2_thetab_p2_ { "model2_thetab_p2", "p_{2}", -0.036 };
	RooChebychev model2_thetab_ {"model2_thetab", "model2_thetab", *thetab_,
			RooArgList(model2_thetab_p0_, model2_thetab_p1_, model2_thetab_p2_)};

	RooRealVar model2_phit_p1_ { "model2_phit_p1", "p_{1}", 0.008};
	RooRealVar model2_phit_p2_ { "model2_phit_p2", "p_{2}", 0.093};
	RooRealVar model2_phit_gaus1_mean_ { "model2_phit_gaus1_mean", "#mu_{1}", -1.488 };
	RooRealVar model2_phit_gaus1_sigma_ { "model2_phit_gaus1_sigma", "#sigma_{1}", 0.475 };
	RooRealVar model2_phit_gaus2_mean_ { "model2_phit_gaus2_mean", "#mu_{2}", 1.841 };
	RooRealVar model2_phit_gaus2_sigma_ { "model2_phit_gaus2_sigma", "#sigma_{2}", 0.390 };
	RooChebychev model2_phit_cheb_ { "model2_phit_cheb", "model2_phit_cheb", *phit_, RooArgList(model2_phit_p1_, model2_phit_p2_) };
	RooGaussian model2_phit_gaus1_ { "model2_phit_gaus1", "model2_phit_gaus1", *phit_, model2_phit_gaus1_mean_, model2_phit_gaus1_sigma_ };
	RooGaussian model2_phit_gaus2_ { "model2_phit_gaus2", "model2_phit_gaus2", *phit_, model2_phit_gaus2_mean_, model2_phit_gaus2_sigma_ };
	RooRealVar model2_f_phit_gaus1_ { "model2_f_phit_gaus1", "f_{G_{1}}", 0.042 };
	RooRealVar model2_f_phit_gaus2_ { "model2_f_phit_gaus2", "f_{G_{2}}", 0.041 };
	RooAddPdf model2_phit_ { "model2_phit", "model2_phit", RooArgList(model2_phit_gaus1_, model2_phit_gaus2_, model2_phit_cheb_),
		RooArgList(model2_f_phit_gaus1_, model2_f_phit_gaus2_) };

	double GetModel2Efficiency() const { return model2_thetat_.getVal() * model2_thetab_.getVal() * model2_phit_.getVal(); }

	RooChebychev model3_thetat_ { "model3_thetat", "model3_thetat", *thetat_, RooArgList() };

	RooRealVar model3_thetab_p0_ { "model3_thetab_p0", "p_{0}", 0.098 };
	RooRealVar model3_thetab_p1_ { "model3_thetab_p1", "p_{1}", -0.339 };
	RooRealVar model3_thetab_p2_ { "model3_thetab_p2", "p_{2}", -0.018 };
	RooChebychev model3_thetab_ { "model3_thetab", "model3_thetab", *thetab_,
		RooArgList(model3_thetab_p0_, model3_thetab_p1_, model3_thetab_p2_) };

	RooRealVar model3_phit_gaus_mean_ { "model3_phit_gaus_mean", "#mu_{G}", 0.16750 };
	RooRealVar model3_phit_gaus_sigma_ { "model3_phit_gaus_sigma", "#sigma_{G}", 0.75956 };
	RooRealVar model3_f_phit_gaus_ { "model3_f_phit_gaus", "f_{G}", -0.43365 };
	RooGenericPdf model3_phit_ { "model3_phit", "model3_phit",
		"1 + model3_f_phit_gaus*TMath::Gaus(phit, model3_phit_gaus_mean, model3_phit_gaus_sigma, 1)",
		RooArgSet(*phit_, model3_f_phit_gaus_, model3_phit_gaus_mean_, model3_phit_gaus_sigma_)};

	double GetModel3Efficiency() const { return model3_thetat_.getVal() * model3_thetab_.getVal() * model3_phit_.getVal(); }

    RooRealVar model4_thetat_p1_{"model4_thetat_p1", "p_{1}", 0};
    RooRealVar model4_thetat_p2_{"model4_thetat_p2", "p_{2}", 0.084};
    RooFormulaVar thetat_plus_pi2{"thetat_plus_pi2", "thetat + pi/2", "thetat - 1.571", RooArgSet(*thetat_)};
    RooPolynomial model4_thetat_{"model4_thetat", "model4_thetat", thetat_plus_pi2, RooArgList{model4_thetat_p1_, model4_thetat_p2_}};

	double GetModel4Efficiency() const { return model4_thetat_.getVal() * model2_thetab_.getVal() * model2_phit_.getVal(); }

    OneDimPhaseSpace phasespace_thetat{"phasespace_thetat", thetat_->getMin(), thetat_->getMax()};
    OneDimPhaseSpace phasespace_thetab{"phasespace_thetab", thetab_->getMin(), thetab_->getMax()};
    OneDimPhaseSpace phasespace_phit{"phasespace_phit", phit_->getMin(), phit_->getMax()};
    CombinedPhaseSpace phasespace{"phasespace", &phasespace_thetat, &phasespace_thetab,
                                  &phasespace_phit};
	BinnedDensity* binned_efficiency = nullptr;

	TH3F* histo_efficiency = nullptr;
	bool CanUseInterpolation(const double& phit, const double& transtht, const double& transthb) const;
	int CloseToEdge(const std::vector<Double_t> vals, const double margin) const;
	double GetHistogramEfficiency(double thetat, double thetab, double phit) const;
	double GetKDEEfficiency(double thetat, double thetab, double phit) const;

	double histo_normalization = 0;
	double GetNormalization();
};
