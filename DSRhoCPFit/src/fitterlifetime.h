/**
 *  @file    fitter.h
 *  @author  Daniel Cervenkov, cervenkov@ipnp.mff.cuni.cz
 *  @date    2015-07-31
 *
 *  @brief This class performs the yield fitting itself as well as plotting
 *
 */

#ifndef FITTERLIFETIME_H_
#define FITTERLIFETIME_H_

// ROOT includes
#include "TCanvas.h"
#include "TPaveText.h"
#include "RooRealVar.h"

// Local includes
#include "constants.h"

class FitterLifetime {
public:
	FitterLifetime();
	virtual ~FitterLifetime();

	void PlotVar(RooRealVar& var, const RooAbsData&) const;
	void PlotWithPull(const RooRealVar& var, const RooAbsData&, const RooAbsPdf& pdf, const char* title = "") const;
	void Test();
	void Toy();

	void SetNumCPUs(const int& numCPUs) {num_CPUs_ = numCPUs;};
	int GetNumCPUs() {return num_CPUs_;};

	void SetDoLifetimeFit(const bool& do_lifetime_fit) {do_lifetime_fit_ = do_lifetime_fit;};
	int GetDoLifetimeFit() const {return do_lifetime_fit_;};

	void SetDoMixingFit(const bool& do_mixing_fit) {do_mixing_fit_ = do_mixing_fit;};
	int GetDoMixingFit() const {return do_mixing_fit_;};

	void SetMakePlots(const bool& make_plots) {make_plots_ = make_plots;};
	int GetMakePlots() const {return make_plots_;};

	void SetPerfectTagging(const bool& perfect_tagging) {perfect_tagging_ = perfect_tagging;};
	int GetPerfectTagging() const {return perfect_tagging_;};

	void ReadInFile(const char* file_path, const int& num_events = 0);
	void SetOutputDir(const char* output_dir);


	RooRealVar* thetat_;
	RooRealVar* thetab_;
	RooRealVar* phit_;

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
	RooRealVar* vreffxi_;
	RooRealVar* vrndf_;
	RooRealVar* vreffndf_;
	RooRealVar* vrntrk_;

	RooRealVar* vtusable_;
	RooRealVar* vtvtxz_;
	RooRealVar* vterr6_;
	RooRealVar* vtchi2_;
	RooRealVar* vtndf_;
	RooRealVar* vtntrk_;
	RooRealVar* vtistagl_;

	RooCategory* decaytype_;


private:
	TPaveText* CreateStatBox(const double chi2, const bool position_top, const bool position_left) const;
	TString GetCommonCutsString() const;

	std::vector<RooRealVar**>conditional_vars_;
	std::vector<RooRealVar**>dataset_vars_;
	RooArgSet conditional_vars_argset_;
	RooArgSet dataset_vars_argset_;

	RooDataSet* dataset_ = NULL;
	RooFitResult* result_ = NULL;

	TFile* output_file_ = NULL;

	int num_CPUs_;
	bool do_lifetime_fit_;
	bool do_mixing_fit_;
	bool make_plots_;
	bool perfect_tagging_;

};

#endif /* FITTERLIFETIME_H_ */
