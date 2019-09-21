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
#include <vector>

// ROOT includes
#include "RooHistPdf.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TH3D.h"
#include "TPaveText.h"
#include "TTree.h"

// Local includes
#include "constants.h"
#include "dtcppdf.h"
#include "dtscfpdf.h"
#include "nlohmann/json.hpp"

class FitterCPV {
   public:
    FitterCPV(nlohmann::json config);
    virtual ~FitterCPV();

    void InitVars(std::array<double, 16> par_input);

    void PlotVar(RooRealVar& var, const RooAbsData&) const;
    void PlotVar(RooRealVar& var, const RooAbsPdf&) const;
    void PlotWithPull(const RooRealVar& var, const RooAbsData&, const RooAbsPdf& pdf,
                      const std::vector<RooAbsPdf*> components = std::vector<RooAbsPdf*>(),
                      const char* title = "") const;

    void Fit(const nlohmann::json config);

    void GenerateToys(const int num_events, const int num_toys);
    void TestEfficiency();
    void PlotEfficiency();

    void SetOutputFile(TFile* file);

    void SetEfficiencyFiles(std::vector<std::string> efficiency_files) {
        efficiency_files_ = efficiency_files;
    };
    std::vector<std::string> GetEfficiencyFiles() { return efficiency_files_; };

    void SetEfficiencyModel(const int& efficiency_model) { efficiency_model_ = efficiency_model; };
    int GetEfficiencyModel() { return efficiency_model_; };

    void SetDoLifetimeFit(const bool& do_lifetime_fit) { do_lifetime_fit_ = do_lifetime_fit; };
    int GetDoLifetimeFit() const { return do_lifetime_fit_; };

    void SetDoTimeIndependentFit(const bool& do_time_independent_fit) {
        do_time_independent_fit_ = do_time_independent_fit;
    };
    int GetDoTimeIndependentFit() const { return do_time_independent_fit_; };

    void SetMakePlots(const bool& make_plots) { make_plots_ = make_plots; };
    int GetMakePlots() const { return make_plots_; };

    void SetPerfectTagging(const bool& perfect_tagging) { perfect_tagging_ = perfect_tagging; };
    int GetPerfectTagging() const { return perfect_tagging_; };

    void SetGeneratorLevel(const bool& generator_level) { generator_level_ = generator_level; };
    int GetGeneratorLevel() const { return generator_level_; };

    bool ResultExists() const { return result_ ? true : false; };

    void ReadInFile(std::vector<const char*> file_names, const int& num_events = 0);
    void SetPlotDir(const char* output_dir);
    void SetSCFKDE(const char* file);
    void SetSCFHisto(const char* file);
    bool FixParameters(const char* pars);
    const std::string CreateResultsString(const bool time_dependent) const;
    const std::string CreatePullTableString(const bool asymmetric = false);
    const std::string CreateLatexPullTableString(const bool asymmetric = false);
    static nlohmann::json ReadJSONConfig(const char* filename);
    std::string ApplyJSONConfig(const nlohmann::json& config);

    const void LogResults();

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
    RooRealVar* de_;
    RooRealVar* csbdtg_;

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
    RooSimultaneous* CreatePDF(const nlohmann::json config);
    RooSimultaneous* CreateChannelPDF(const std::string channel_name,
                                      const nlohmann::json channel_config,
                                      const nlohmann::json common_config);
    RooDataSet* GetData(const nlohmann::json config);
    RooDataSet* GetChannelData(const std::string channel_name, const nlohmann::json channel_config,
                               const nlohmann::json common_config);

    void PrepareVarArgsets();
    void ChangeFitRanges(const nlohmann::json& config);
    void ChangeModelParameters(RooAbsPdf* pdf, const std::string channel_name,
                               const nlohmann::json& config);
    TPaveText* CreateStatBox(const double chi2, const int ndof, const bool position_top,
                             const bool position_left) const;
    TH3D* GetBinnedEfficiency(std::vector<std::string> file, const int model);
    const void SaveLikelihoodScan(RooAbsPdf& pdf, RooRealVar* var, const double margin = 0);
    const void SaveLikelihoodScan(RooAbsPdf& pdf, RooRealVar* var1, RooRealVar* var2,
                                  const double margin1 = 0, const double margin2 = 0);
    const double Calculate3DChi2(const RooDataHist& data, const RooDataHist& pdf);
    const void SaveChi2Scan(RooSimultaneous& pdf, RooRealVar* var, const double margin = 0);
    int CloseToEdge(const std::vector<Double_t> vals, const double margin) const;

    // RooAbsPdf* CreateAngularSCFPDF(const std::string channel_name);
    RooAbsPdf* CreateAngularSCFBKGPDF(const std::string prefix);

    RooAbsPdf* CreateSCFPDF(const std::string channel_name,
                            const nlohmann::json channel_config) const;
    RooAbsPdf* GetHistoSCF(const std::string filename) const;
    RooSimultaneous* CreateAngularPDF(const std::string name_prefix, const bool scf, const bool bkg,
                                      const nlohmann::json channel_config);
    RooSimultaneous* CreateTimeDependentPDF(const std::string channel_name,
                                            const nlohmann::json common_config,
                                            const nlohmann::json channel_config);

    RooAddPdf* CreateVoigtGaussDtPdf(const std::string prefix);
    void CreateDtCPPDFs(DtCPPDF*& cr_pdf_FB, DtCPPDF*& cr_pdf_FA, DtCPPDF*& cr_pdf_SB,
                        DtCPPDF*& cr_pdf_SA, const std::string channel_name,
                        const nlohmann::json common_config,
                        const nlohmann::json channel_config) const;
    void CreateDtSCFPDFs(DtSCFPDF*& scf_pdf_FB, DtSCFPDF*& scf_pdf_FA, DtSCFPDF*& scf_pdf_SB,
                         DtSCFPDF*& scf_pdf_SA, const std::string channel_name) const;
    void CreateFunctionalDtSCFBKGPDFs(RooProdPdf*& pdf_FB, RooProdPdf*& pdf_FA, RooProdPdf*& pdf_SB,
                                      RooProdPdf*& pdf_SA, const std::string channel_name,
                                      const nlohmann::json channel_config, const bool scf);

    void PlotFit(RooSimultaneous* pdf, const bool scf, const bool bkg);

    RooRealVar cr_scf_f_{"cr_scf_f", "f_{cr}", constants::fraction_cr_of_crscf, 0.80, 0.99};
    RooRealVar cr_f_{"cr_f", "f_{cr}", constants::fraction_cr_of_crscfbkg, 0.10, 0.99};
    RooRealVar scf_f_{"scf_f", "f_{scf}", constants::fraction_scf_of_crscfbkg, 0.10, 0.99};

    RooDataHist* scf_angular_kde_hist_ = nullptr;
    RooHistPdf* scf_angular_kde_ = nullptr;

    std::vector<RooRealVar**> conditional_vars_;
    std::vector<RooRealVar**> dataset_vars_;
    std::vector<RooRealVar**> parameters_;
    RooArgSet conditional_vars_argset_;
    RooArgSet dataset_vars_argset_;
    RooArgSet parameters_argset_;

    std::array<double, 16> par_input_;

    RooDataSet* dataset_ = nullptr;
    RooFitResult* result_ = nullptr;

    TFile* output_file_ = nullptr;

    std::vector<std::string> efficiency_files_;
    int efficiency_model_;
    bool do_lifetime_fit_;
    bool do_time_independent_fit_;
    bool make_plots_;
    bool perfect_tagging_;
    bool generator_level_;

    RooSimultaneous* pdf_;
    RooDataSet* data_;
    nlohmann::json config_;
};

#endif /* FITTERCPV_H_ */
