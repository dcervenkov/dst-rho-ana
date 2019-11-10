/**
 *  @file    fitter.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2015-07-31
 *
 *  @brief This class performs the yield fitting itself as well as plotting
 *
 */

#ifndef FITTERLIFETIME_H_
#define FITTERLIFETIME_H_

// ROOT includes
#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "TCanvas.h"
#include "TPaveText.h"

// Local includes
#include "constants.h"
#include "nlohmann/json.hpp"

class FitterLifetime {
   public:
    FitterLifetime(const nlohmann::json config);
    virtual ~FitterLifetime();

    void PlotVar(RooRealVar& var) const;
    void PlotWithPull(const RooRealVar& var, const RooAbsData&, const RooAbsPdf& pdf,
                      const std::vector<RooAbsPdf*>& components = {}, const char* title = "") const;
    void Test();

    void SetNumCPUs(const int& numCPUs) { num_CPUs_ = numCPUs; };
    int GetNumCPUs() { return num_CPUs_; };

    void SetChannel(const std::string channel) { channel_ = channel; };
    std::string GetChannel() { return channel_; };

    void SetDoLifetimeFit(const bool& do_lifetime_fit) { do_lifetime_fit_ = do_lifetime_fit; };
    int GetDoLifetimeFit() const { return do_lifetime_fit_; };

    void SetDoMixingFit(const bool& do_mixing_fit) { do_mixing_fit_ = do_mixing_fit; };
    int GetDoMixingFit() const { return do_mixing_fit_; };

    void SetMakePlots(const bool& make_plots) { make_plots_ = make_plots; };
    int GetMakePlots() const { return make_plots_; };

    void SetPerfectTagging(const bool& perfect_tagging) { perfect_tagging_ = perfect_tagging; };
    int GetPerfectTagging() const { return perfect_tagging_; };

	void ReadInFile(const std::vector<const char*> file_names, const int& num_events = 0);
    void SaveTXTResults(const char* results_file) const;

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

    RooRealVar* dt_;

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
    TPaveText* CreateStatBox(const double chi2, const bool position_top,
                             const bool position_left) const;
    RooAddPdf* CreateVoigtGaussDtPdf(const std::string prefix) const;
    RooAbsPdf* CreateLifetimePDF(std::vector<RooAbsPdf*>& components, const bool scf, const bool bkg) const;
    nlohmann::json ReadInJSONFile(const char* filename) const;

    std::vector<RooRealVar**> conditional_vars_;
    std::vector<RooRealVar**> dataset_vars_;
    const RooArgSet* vars;
    RooArgSet conditional_argset_;
    RooArgSet dataset_argset_;

    RooDataSet* dataset_ = nullptr;
    RooFitResult* result_ = nullptr;

    TFile* output_file_ = nullptr;

    int num_CPUs_;
    bool do_lifetime_fit_;
    bool do_mixing_fit_;
    bool make_plots_;
    bool perfect_tagging_;

    nlohmann::json config_;
    std::string channel_;
};

#endif /* FITTERLIFETIME_H_ */
