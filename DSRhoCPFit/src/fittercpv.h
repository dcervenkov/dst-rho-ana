/**
 *  @file    fittercpv.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2016-11-03
 *
 *  @brief This class performs the yield fitting itself as well as plotting
 *
 */

#ifndef FITTERCPV_H_
#define FITTERCPV_H_

// Standard includes
#include <array>

// ROOT includes
#include "RooRealVar.h"
#include "TCanvas.h"
#include "TH3D.h"
#include "TPaveText.h"

// Local includes
#include "constants.h"

class FitterCPV {
   public:
    FitterCPV(std::array<double, 16> par_input);
    virtual ~FitterCPV();

    void PlotVar(RooRealVar& var, const RooAbsData&) const;
    void PlotWithPull(const RooRealVar& var, const RooAbsData&, const RooAbsPdf& pdf,
                      const std::vector<RooAbsPdf*> components = std::vector<RooAbsPdf*>(),
                      const char* title = "") const;
    void FitSignal();
    void FitSCF();
    void FitAll();
    void FitAngularCR();
    void GenerateToys(const int num_events, const int num_toys);
    void TestEfficiency();
    void PlotEfficiency();

    void SetNumCPUs(const int& numCPUs) { num_CPUs_ = numCPUs; };
    int GetNumCPUs() { return num_CPUs_; };

    void SetEfficiencyModel(const int& efficiency_model) { efficiency_model_ = efficiency_model; };
    int GetEfficiencyModel() { return efficiency_model_; };

    void SetDoLifetimeFit(const bool& do_lifetime_fit) { do_lifetime_fit_ = do_lifetime_fit; };
    int GetDoLifetimeFit() const { return do_lifetime_fit_; };

    void SetDoMixingFit(const bool& do_mixing_fit) { do_mixing_fit_ = do_mixing_fit; };
    int GetDoMixingFit() const { return do_mixing_fit_; };

    void SetDoTimeIndependentFit(const bool& do_time_independent_fit) {
        do_time_independent_fit_ = do_time_independent_fit;
    };
    int GetDoTimeIndependentFit() const { return do_time_independent_fit_; };

    void SetMakePlots(const bool& make_plots) { make_plots_ = make_plots; };
    int GetMakePlots() const { return make_plots_; };

    void SetPerfectTagging(const bool& perfect_tagging) { perfect_tagging_ = perfect_tagging; };
    int GetPerfectTagging() const { return perfect_tagging_; };

    void ReadInFile(const char* file_path, const int& num_events = 0);
    void SetPlotDir(const char* output_dir);
    bool FixParameters(const char* pars);
    bool SaveResults(const char* file);
    std::string SetupFitRange(const char* filename);
    RooDataSet* ReduceDataToFitRange(const char* filename);

    RooRealVar* ap_;
    RooRealVar* apa_;
    RooRealVar* a0_;
    RooRealVar* a0a_;
    RooFormulaVar* at_;
    RooRealVar* ata_;

    RooRealVar* xp_;
    RooRealVar* x0_;
    RooRealVar* xt_;

    RooRealVar* yp_;
    RooRealVar* y0_;
    RooRealVar* yt_;

    RooRealVar* xpb_;
    RooRealVar* x0b_;
    RooRealVar* xtb_;

    RooRealVar* ypb_;
    RooRealVar* y0b_;
    RooRealVar* ytb_;

    RooRealVar* dt_;
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

   private:
    TPaveText* CreateStatBox(const double chi2, const bool position_top,
                             const bool position_left) const;
    TString GetCommonCutsString() const;
    TH3D* GetBinnedEfficiency();

    std::vector<RooRealVar**> conditional_vars_;
    std::vector<RooRealVar**> dataset_vars_;
    std::vector<RooRealVar**> parameters_;
    RooArgSet conditional_vars_argset_;
    RooArgSet dataset_vars_argset_;
    RooArgSet parameters_argset_;

    std::array<double, 16> par_input_;

    RooDataSet* dataset_ = NULL;
    RooFitResult* result_ = NULL;

    TFile* output_file_ = NULL;

    int num_CPUs_;
    int efficiency_model_;
    bool do_lifetime_fit_;
    bool do_mixing_fit_;
    bool do_time_independent_fit_;
    bool make_plots_;
    bool perfect_tagging_;
};

#endif /* FITTERCPV_H_ */
