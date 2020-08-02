/**
 *  @file    fitterbkg.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2017-12-15
 *
 *  @brief Class that performs the fit itself as well as plotting
 *
 */

#include "fitterbkg.h"

// Standard includes
#include <array>
#include <boost/filesystem.hpp>
#include <sstream>
#include <string>

// ROOT includes
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooHist.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooTreeDataStore.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TEnv.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TTree.h"

// Meerkat includes
#include "AbsDensity.hh"
#include "AdaptiveKernelDensity.hh"
#include "BinnedKernelDensity.hh"
#include "CombinedPhaseSpace.hh"
#include "FormulaDensity.hh"
#include "KernelDensity.hh"
#include "OneDimPhaseSpace.hh"
#include "RooMeerkatPdf.hh"
#include "UniformDensity.hh"

// Local includes
#include "constants.h"
#include "log.h"
#include "tools.h"

FitterBKG::FitterBKG() {
    vrusable_ = new RooRealVar("vrusable", "vrusable", 0, 1);
    vrvtxz_ = new RooRealVar("vrvtxz", "vrvtxz", -10, 10);
    vrerr6_ = new RooRealVar("vrerr6", "vrerr6", -1, 1);
    vrchi2_ = new RooRealVar("vrchi2", "vrchi2", 0, 10000000);
    //  vreffxi_ = new RooRealVar( "vreffxi", "vreffxi", 0, 10000000);
    vrndf_ = new RooRealVar("vrndf", "vrndf", 0, 100);
    //  vreffndf_ = new RooRealVar( "vreffndf", "vreffndf", 0, 100);
    vrntrk_ = new RooRealVar("vrntrk", "vrntrk", 0, 100);

    vtusable_ = new RooRealVar("vtusable", "vtusable", 1);
    vtvtxz_ = new RooRealVar("vtvtxz", "vtvtxz", -10, 10);
    vtchi2_ = new RooRealVar("vtchi2", "vtchi2", 1.6);
    vtndf_ = new RooRealVar("vtndf", "vtndf", 4);
    vterr6_ = new RooRealVar("vterr6", "vterr6", -1, 1);
    vtchi2_ = new RooRealVar("vtchi2", "vtchi2", 0, 10000000);
    vtndf_ = new RooRealVar("vtndf", "vtndf", 0, 100);
    vtntrk_ = new RooRealVar("vtntrk", "vtntrk", 0, 100);
    vtistagl_ = new RooRealVar("vtistagl", "vtistagl", 0, 100);

    expno_ = new RooRealVar("expno", "expno", 0, 100);
    expmc_ = new RooRealVar("expmc", "expmc", 0, 10);

    evmcflag_ = new RooRealVar("evmcflag", "evmcflag", -2, 10);
    brecflav_ = new RooRealVar("brecflav", "brecflav", -1, 1);
    btagmcli_ = new RooRealVar("btagmcli", "btagmcli", -1000, 1000);
    tagqr_ = new RooRealVar("tagqr", "tagqr", -1, 1);
    tagwtag_ = new RooRealVar("tagwtag", "tagwtag", 0, 0.5);

    benergy_ = new RooRealVar("benergy", "benergy", 0, 100);
    mbc_ = new RooRealVar("mbc", "mbc", 5, 6);
    de_ = new RooRealVar("de", "de", -1, 1);
    csbdtg_ = new RooRealVar("csbdtg", "csbdtg", -1, 1);

    shcosthb_ = new RooRealVar("shcosthb", "shcosthb", -1, 1);

    nocand_ = new RooRealVar("nocand", "nocand", 1, 6);

    decaytype_ = new RooCategory("decaytype", "decaytype");
    decaytype_->defineType("FB", 1);
    decaytype_->defineType("FA", 2);
    decaytype_->defineType("SB", 3);
    decaytype_->defineType("SA", 4);

    // All the variables present in the ntuple are added to a vector so we can later
    // iterate through them for convenience
    conditional_vars_.push_back(&expno_);
    conditional_vars_.push_back(&expmc_);

    conditional_vars_.push_back(&evmcflag_);
    conditional_vars_.push_back(&brecflav_);
    conditional_vars_.push_back(&btagmcli_);
    conditional_vars_.push_back(&tagqr_);
    conditional_vars_.push_back(&tagwtag_);

    conditional_vars_.push_back(&benergy_);
    conditional_vars_.push_back(&mbc_);
    conditional_vars_.push_back(&de_);
    conditional_vars_.push_back(&csbdtg_);

    conditional_vars_.push_back(&shcosthb_);

    conditional_vars_.push_back(&vrusable_);
    conditional_vars_.push_back(&vrvtxz_);
    conditional_vars_.push_back(&vrerr6_);
    conditional_vars_.push_back(&vrchi2_);
    //  conditional_vars_.push_back(&vreffxi_);
    conditional_vars_.push_back(&vrndf_);
    //  conditional_vars_.push_back(&vreffndf_);
    conditional_vars_.push_back(&vrntrk_);

    conditional_vars_.push_back(&vtusable_);
    conditional_vars_.push_back(&vtvtxz_);
    conditional_vars_.push_back(&vterr6_);
    conditional_vars_.push_back(&vtchi2_);
    conditional_vars_.push_back(&vtndf_);
    conditional_vars_.push_back(&vtntrk_);

    conditional_vars_.push_back(&vtistagl_);

    conditional_vars_.push_back(&nocand_);

    dataset_vars_ = conditional_vars_;
    // dataset_vars_.push_back(&thetat_);
    // dataset_vars_.push_back(&thetab_);
    // dataset_vars_.push_back(&phit_);

    // The variables present in the ntuple are added to an RooArgSet that will be needed
    // when we create a RooDataSet from the input_tree
    for (auto var : conditional_vars_) {
        conditional_vars_argset_.add(**var);
    }

    // TODO: Comments
    for (auto var : dataset_vars_) {
        dataset_vars_argset_.add(**var);
    }
    dataset_vars_argset_.add(dt_);
    dataset_vars_argset_.add(thetat_);
    dataset_vars_argset_.add(thetab_);
    dataset_vars_argset_.add(phit_);

    num_CPUs_ = 1;

    physics_dt_model_ =
        new DtBKG("physics_dt_model", "physics_dt_model", dt_, *vrerr6_, *vterr6_, dt_tau_,
                  dt_f_delta_, dt_mu_delta_, dt_mu_lifetime_, dt_f_tail_, dt_S_main_, dt_S_tail_);
}

FitterBKG::~FitterBKG() {
    if (output_file_) {
        output_file_->Close();
    }
}

/**
 * Make a simple plot of a variable and save it to a file
 * NOTE: @var can't be const (without const_cast) because of dumb RooFit design
 */
void FitterBKG::PlotVar(RooRealVar& var, const RooDataSet* data) const {
    TCanvas* canvas = new TCanvas(var.GetName(), var.GetTitle(), 500, 500);

    double range_min;
    double range_max;
    data->getRange(var, range_min, range_max, 0.05);

    RooPlot* plot = var.frame(range_min, range_max);

    data->plotOn(plot);
    plot->Draw();

    plot->SetTitle("");
    plot->GetXaxis()->SetTitle(TString(var.GetTitle()));
    plot->GetYaxis()->SetTitle("");

    canvas->Write();
    canvas->SaveAs(constants::format);
}

void FitterBKG::PlotWithPull(RooRealVar& var, const RooDataSet& data, const RooAbsPdf& pdf) const {
    tools::PlotWithPull(var, conditional_vars_argset_, data, pdf,
                        result_ ? result_->floatParsFinal() : RooArgList(),
                        std::vector<RooAbsPdf*>(), num_CPUs_);
}

/**
 * Reads in data from ROOT file(s). Constructs separate datasets for the 4 categories.
 * Binds the variables to the dataset, so that dataset->get(i) changes values of, e.g., expno_
 *
 * @param file_names Vector of paths to the ROOT files
 * @param num_events [optional] Maximum number of events to use (0 to read all)
 */
void FitterBKG::ReadInFile(std::vector<const char*> file_names, const int& num_events) {
    input_tree = new TChain("h2000");
    for (auto file_name : file_names) {
        input_tree->Add(file_name);
    }

    // if (num_events) {
    //     TTree* temp_tree = input_tree;
    //     input_tree = temp_tree->CloneTree(num_events);
    //     delete temp_tree;
    // }

    TString common_cuts = tools::GetCommonCutsString();
    common_cuts += "&&evmcflag!=1";

    TString FB_cuts;
    TString FA_cuts;
    TString SB_cuts;
    TString SA_cuts;
    FB_cuts = "brecflav==1&&tagqr<0";
    FA_cuts = "brecflav==-1&&tagqr>0";
    SB_cuts = "brecflav==1&&tagqr>0";
    SA_cuts = "brecflav==-1&&tagqr<0";

    Log::print(Log::info, "Reading %i input files...\n", file_names.size());
    dataset_ = new RooDataSet("dataset", "dataset", input_tree, dataset_vars_argset_, common_cuts);

    Log::print(Log::info, "Creating FB, FA, SB, SA labels...\n");
    // We add an identifying label to each of the 4 categories and then combine it into a single
    // dataset for RooSimultaneous fitting
    dataset_FB_ = static_cast<RooDataSet*>(dataset_->reduce(FB_cuts));
    decaytype_->setLabel("FB");
    dataset_FB_->addColumn(*decaytype_);
    dataset_FB_->SetName("dataset_FB");

    dataset_FA_ = static_cast<RooDataSet*>(dataset_->reduce(FA_cuts));
    decaytype_->setLabel("FA");
    dataset_FA_->addColumn(*decaytype_);
    dataset_FA_->SetName("dataset_FA");

    dataset_SB_ = static_cast<RooDataSet*>(dataset_->reduce(SB_cuts));
    decaytype_->setLabel("SB");
    dataset_SB_->addColumn(*decaytype_);
    dataset_SB_->SetName("dataset_SB");

    dataset_SA_ = static_cast<RooDataSet*>(dataset_->reduce(SA_cuts));
    decaytype_->setLabel("SA");
    dataset_SA_->addColumn(*decaytype_);
    dataset_SA_->SetName("dataset_SA");

    data_tree = input_tree->CopyTree(common_cuts);
    delete input_tree;
    // input_file->Close();

    // Bind the variables to the dataset, so that dataset->get(i) changes values of, e.g., expno_
    const RooArgSet* vars = dataset_->get();
    for (RooRealVar** var : dataset_vars_) {
        *var = static_cast<RooRealVar*>(vars->find((*var)->GetName()));
    }

    dataset_->get(1)->Print("v");

    // Set the current directory back to the one for writing (ugly ROOT stuff)
    if (output_file_) {
        output_file_->cd();
    }
    Log::print(Log::info, "Input dataset ready\n");
}

/**
 * Set the directory to which to ouput plots
 */
void FitterBKG::SetPlotDir(const char* plot_dir) {
    if (!boost::filesystem::is_directory(plot_dir)) {
        boost::filesystem::create_directories(plot_dir);
    }
    TString root_filename(boost::filesystem::path(plot_dir).filename().string());
    root_filename += ".root";
    gEnv->SetValue("Canvas.PrintDirectory", plot_dir);
    Log::print(Log::info, "print dir: %s\n", gEnv->GetValue("Canvas.PrintDirectory", "not found"));
    output_file_ = new TFile(TString(plot_dir) + "/" + root_filename, "RECREATE");
}

/**
 * Fit a given PDF to a given dataset
 *
 * All the variables are already defined in the PDF.
 */
void FitterBKG::Fit(RooAbsPdf* pdf, RooDataSet* data) {
    result_ =
        pdf->fitTo(*data, RooFit::Save(), RooFit::Minimizer("Minuit2"), RooFit::NumCPU(num_CPUs_));
}

void FitterBKG::CreateHistoPDF(RooDataSet* data, const std::string results_file,
                               const bool randomize) {
    thetat_.setBins(50);
    thetab_.setBins(50);
    phit_.setBins(50);
    RooDataHist data_hist("data_hist", "data_hist", RooArgSet(thetat_, thetab_, phit_), *data);
    RooDataHist* rnd_hist = 0;
    if (randomize) {
        rnd_hist = tools::RandomizeDataHist(data_hist);
    }
    RooHistPdf* scf_hist_pdf =
        new RooHistPdf("scf_hist_pdf", "scf_hist_pdf", RooArgSet(thetat_, thetab_, phit_),
                       randomize ? *rnd_hist : data_hist);

    tools::CreateDirsIfNecessary(results_file);
    TFile f(results_file.c_str(), "RECREATE");
    scf_hist_pdf->Write();
    f.Close();
    delete scf_hist_pdf;
    if (rnd_hist) delete rnd_hist;
}

/**
 * Replace bins not in (low,high) range by average of neighbours
 *
 * The neighbourhood taken for averaging increases from nearest neighbours
 * until the average is in the (low,high) range.
 *
 * @param pdf The histogram PDF to modify
 * @param low Low edge of the required range
 * @param high High edge of the required range
 *
 * @return TH3F* Modified histogram PDF with all bins in the (low,high) range
 */
TH3F* FitterBKG::NormalizePDF(const TH3F* pdf, const double low, const double high) {
    int total_bins = 0;
    int fixed_bins = 0;
    TH3F* normalized_pdf = (TH3F*)pdf->Clone();
    for (int x = 1; x <= pdf->GetNbinsX(); x++) {
        for (int y = 1; y <= pdf->GetNbinsY(); y++) {
            for (int z = 1; z <= pdf->GetNbinsZ(); z++) {
                total_bins++;
                int bin = pdf->GetBin(x, y, z);
                if (pdf->GetBinContent(bin) > high || pdf->GetBinContent(bin) < low) {
                    fixed_bins++;
                    int size = 1;
                    double interpolation = -1;
                    do {
                        interpolation = Interpolate(pdf, x, y, z, size++);
                    } while (interpolation < 0 || interpolation > 1);
                    // printf("bin %i: value = %f, interpolation = %f\n", bin,
                    // pdf->GetBinContent(bin), interpolation);
                    normalized_pdf->SetBinContent(bin, interpolation);
                }
            }
        }
    }

    Log::print(Log::info, "Fixed %i/%i (%.2f%%) bins.\n", fixed_bins, total_bins,
               (double)fixed_bins / total_bins * 100);
    return normalized_pdf;
}

/**
 * @brief Get average value of a group of neighboring cells in a histogram
 *
 * The central cell is not excluded from the calculation.
 *
 * @param histo Histogram to process
 * @param x_org Number of the cell around which the neighbors are chosen (x-axis)
 * @param y_org Number of the cell around which the neighbors are chosen (y-axis)
 * @param z_org Number of the cell around which the neighbors are chosen (z-axis)
 * @param size How many neighbours are to be chosen to each side on each axis
 *
 * @return double The average values of the neighbor group
 */
double FitterBKG::Interpolate(const TH3F* histo, int x_org, int y_org, int z_org, int size) {
    double new_value = 0;
    int points = 0;
    int num_bins_x = histo->GetXaxis()->GetNbins();
    int num_bins_y = histo->GetYaxis()->GetNbins();
    int num_bins_z = histo->GetZaxis()->GetNbins();

    for (int x = x_org - size; x <= x_org + size; x++) {
        for (int y = y_org - size; y <= y_org + size; y++) {
            for (int z = z_org - size; z <= z_org + size; z++) {
                if (x == x_org && y == y_org && z == z_org) continue;
                if (x < 1 || y < 1 || z < 1 || x > num_bins_x || y > num_bins_y || z > num_bins_z) {
                    printf("skipping\n");
                    continue;
                }
                int bin = histo->GetBin(x, y, z);
                new_value += histo->GetBinContent(bin);
                // printf("delta = %f\n", histo->GetBinContent(bin));
                points++;
            }
        }
    }
    if (points == 0) return 0;
    return new_value / points;
}

AdaptiveKernelDensity FitterBKG::CreateKDEPDF(RooDataSet* data, const std::string results_file) {
    OneDimPhaseSpace* phasespace_thetat =
        new OneDimPhaseSpace("phasespace_thetat", thetat_.getMin(), thetat_.getMax());
    OneDimPhaseSpace* phasespace_thetab =
        new OneDimPhaseSpace("phasespace_thetab", thetab_.getMin(), thetab_.getMax());
    OneDimPhaseSpace* phasespace_phit =
        new OneDimPhaseSpace("phasespace_phit", phit_.getMin(), phit_.getMax());
    CombinedPhaseSpace* phasespace =
        new CombinedPhaseSpace("phasespace", phasespace_thetat, phasespace_thetab, phasespace_phit);
    // UniformDensity uniform_density("Uniform Density", &phasespace);
    // FormulaDensity formula_density("Formula Density", &phasespace,
    //                                "(x - 1.57)^2 * (y - 1.57)^2 * (z)^2");

    // This simple step function provides better edge approximation
    char formula[1000];
    std::snprintf(formula, 1000,
                  "(x > %f ? 1 : 0) * (x < %f ? 1 : 0) * "
                  "(y > %f ? 1 : 0) * (y < %f ? 1 : 0) * "
                  "(z > %f ? 1 : 0) * (z < %f ? 1 : 0)",
                  constants::cuts::thetat_low, constants::cuts::thetat_high,
                  constants::cuts::thetab_low, constants::cuts::thetab_high,
                  constants::cuts::phit_low, constants::cuts::phit_high);

    FormulaDensity formula_density("Formula Density", phasespace, formula);

    // RooAbsData::setDefaultStorageType(RooAbsData::Tree);
    // RooDataSet* dataNew = new RooDataSet("dataNew","dataNew", data, *data->get());
    // TTree* eff_tree = const_cast<TTree*>(dataNew->tree())->CloneTree();

    // TTree* eff_tree = data->GetClonedTree();
    TTree* eff_tree = data_tree;

    // std::array<double, 6> bin_kde_pars = {50, 50, 50, 0.2, 0.2, 0.4};
    std::array<double, 6> ada_kde_pars = {100, 100, 100, 0.1, 0.05, 0.1};

    // BinnedKernelDensity bin_kde("BinKernelPDF", phasespace, eff_tree, "thetat", "thetab", "phit",
    //                             bin_kde_pars[0], bin_kde_pars[1], bin_kde_pars[2],
    //                             bin_kde_pars[3], bin_kde_pars[4], bin_kde_pars[5],
    //                             &formula_density);

    AdaptiveKernelDensity kde("KernelPDF", phasespace, eff_tree, "thetat", "thetab", "phit",
                              ada_kde_pars[0], ada_kde_pars[1], ada_kde_pars[2], ada_kde_pars[3],
                              ada_kde_pars[4], ada_kde_pars[5], &formula_density);
    //   ada_kde_pars[3], ada_kde_pars[4], ada_kde_pars[5], &bin_kde);

    tools::CreateDirsIfNecessary(results_file);
    kde.writeToTextFile(results_file.c_str());

    return kde;
}

TH3F* FitterBKG::ConvertDensityToHisto(AdaptiveKernelDensity pdf) const {
    const int num_bins = 50;

    assert(num_bins == thetat_.getBins());
    assert(num_bins == thetab_.getBins());
    assert(num_bins == phit_.getBins());

    TH3F* pdf_histo =
        new TH3F("pdf_histo", "pdf_histo", num_bins, thetat_.getMin(), thetat_.getMax(), num_bins,
                 thetab_.getMin(), thetab_.getMax(), num_bins, phit_.getMin(), phit_.getMax());

    double thetat;
    double thetab;
    double phit;
    for (int x = 1; x <= pdf_histo->GetXaxis()->GetNbins(); x++) {
        for (int y = 1; y <= pdf_histo->GetYaxis()->GetNbins(); y++) {
            for (int z = 1; z <= pdf_histo->GetZaxis()->GetNbins(); z++) {
                thetat = pdf_histo->GetXaxis()->GetBinCenter(x);
                thetab = pdf_histo->GetYaxis()->GetBinCenter(y);
                phit = pdf_histo->GetZaxis()->GetBinCenter(z);
                std::vector<double> coords;
                coords.push_back(thetat);
                coords.push_back(thetab);
                coords.push_back(phit);
                pdf_histo->SetBinContent(pdf_histo->GetBin(x, y, z), pdf.density(coords));
            }
        }
    }

    return pdf_histo;
}

TH3F* FitterBKG::Create3DHisto(const RooDataSet* dataset) const {
    const int num_bins = 50;

    assert(num_bins == thetat_.getBins());
    assert(num_bins == thetab_.getBins());
    assert(num_bins == phit_.getBins());

    TH3F* histo = new TH3F(dataset->GetName(), dataset->GetTitle(), num_bins, thetat_.getMin(),
                           thetat_.getMax(), num_bins, thetab_.getMin(), thetab_.getMax(), num_bins,
                           phit_.getMin(), phit_.getMax());
    for (int i = 0; i < dataset->numEntries(); i++) {
        const RooArgSet* row = dataset->get(i);
        double thetat = row->getRealValue("thetat");
        double thetab = row->getRealValue("thetab");
        double phit = row->getRealValue("phit");
        histo->Fill(thetat, thetab, phit);
    }

    return histo;
}

void FitterBKG::PlotKDE(AdaptiveKernelDensity kde) const {
    TH3F* model = ConvertDensityToHisto(kde);
    TH3F* data = Create3DHisto(dataset_);
    const double scale = data->GetEntries() / model->GetEntries();

    RooDataHist roo_model("roo_model", "roo_model", RooArgList(thetat_, thetab_, phit_), model,
                          scale);
    RooDataHist roo_data("roo_data", "roo_data", RooArgList(thetat_, thetab_, phit_), data);

    RooRealVar vars[] = {thetat_, thetab_, phit_};

    RooArgList list(thetat_, thetab_, phit_);
    RooMeerkatPdf meerkat_pdf("meerkat_pdf", "meerkat_pdf", list, &kde);
    RooHistPdf meerkat_histpdf("meerkat_histpdf", "meerkat_histpdf",
                               RooArgSet(thetat_, thetab_, phit_), roo_model);
    for (auto var : vars) {
        // PlotWithPull(var, *dataset_, meerkat_pdf);
        tools::PlotWithPull(var, conditional_vars_argset_, *dataset_, meerkat_histpdf,
                            result_ ? result_->floatParsFinal() : RooArgList());
    }

    for (int i = 0; i < 3; i++) {
        for (int j = i + 1; j < 3; j++) {
            tools::PlotVars2D(vars[i], vars[j], roo_model);
            tools::PlotVars2D(vars[i], vars[j], roo_data);
            tools::PlotPull2D(vars[i], vars[j], roo_model, roo_data);
        }
    }
}

/**
 * Disable certain parts of the physics-based PDF
 */
void FitterBKG::SetNoDeltaPDF() {
    dt_f_delta_ = 0;
    dt_f_delta_.setConstant();
    dt_mu_delta_ = 0;
    dt_mu_delta_.setConstant();
}

/**
 * Disable certain parts of the physics-based PDF
 */
void FitterBKG::SetNoTailPDF() {
    dt_f_tail_ = 0;
    dt_f_tail_.setConstant();
    dt_S_tail_ = 0;
    dt_S_tail_.setConstant();
}

/**
 * Fit (and plot) a Delta t model
 *
 * @param nodelta Don't include delta function term in PDF
 * @param notail Don't include resolution tail term in PDF
 * @param plot Create plots
 *
 * @return nlohmann::json Fitted parameters and values
 */
nlohmann::json FitterBKG::FitDt(RooAbsPdf* model, std::string prefix, bool plot) {
    Fit(model, dataset_);

    nlohmann::json json_results;
    json_results[prefix + "dt_model"] = tools::GetResultsJSON(model, RooArgSet(dt_), "bkg_");
    if (plot) {
        PlotWithPull(dt_, *dataset_, *model);
    }

    {
        RooDataSet* dataset_cf =
            static_cast<RooDataSet*>(dataset_->emptyClone("dataset_cf", "dataset_cf"));
        dataset_cf->append(*dataset_FB_);
        dataset_cf->append(*dataset_FA_);
        Fit(model, dataset_cf);
        json_results[prefix + "dt_cf_model"] =
            tools::GetResultsJSON(model, RooArgSet(dt_), "bkg_cf_");
        if (plot) {
            PlotWithPull(dt_, *dataset_cf, *model);
        }
        delete dataset_cf;
    }

    {
        RooDataSet* dataset_dcs =
            static_cast<RooDataSet*>(dataset_->emptyClone("dataset_dcs", "dataset_dcs"));
        dataset_dcs->append(*dataset_SB_);
        dataset_dcs->append(*dataset_SA_);
        Fit(model, dataset_dcs);
        json_results[prefix + "dt_dcs_model"] =
            tools::GetResultsJSON(model, RooArgSet(dt_), "bkg_dcs_");
        if (plot) {
            PlotWithPull(dt_, *dataset_dcs, *model);
        }
        delete dataset_dcs;
    }
    return json_results;
}

/**
 * Fit (and plot) angular distributions
 *
 * @param plot Create plots
 *
 * @return nlohmann::json Fitted parameters and values
 */
nlohmann::json FitterBKG::FitAngular(bool plot) {
    // We must plot after each fit otherwise the results printed on the plot
    // would be of the last fit
    Fit(&phit_model_, dataset_);
    if (plot) {
        PlotWithPull(phit_, *dataset_, phit_model_);
    }

    Fit(&thetat_model_, dataset_);
    if (plot) {
        PlotWithPull(thetat_, *dataset_, thetat_model_);
    }

    Fit(&thetab_model_, dataset_);
    if (plot) {
        PlotWithPull(thetab_, *dataset_, thetab_model_);
    }

    std::vector<const RooAbsPdf*> angular_pdfs = {&thetab_model_, &thetat_model_, &phit_model_};
    nlohmann::json json_results;
    json_results["angular_pdf"] =
        tools::GetResultsJSON(angular_pdfs, RooArgSet(thetat_, thetab_, phit_), "bkg_");

    if (plot) {
        thetat_.setBins(50);
        thetab_.setBins(50);
        phit_.setBins(50);

        RooProdPdf* bkg_pdf = new RooProdPdf("bkg_pdf", "bkg_pdf",
                                             RooArgList(thetat_model_, thetab_model_, phit_model_));

        RooDataHist* pdf_hist = bkg_pdf->generateBinned(RooArgSet(thetat_, thetab_, phit_),
                                                        dataset_->numEntries(), true);

        RooDataHist data_hist("data_hist", "data_hist", RooArgSet(thetat_, thetab_, phit_),
                              *dataset_);

        tools::PlotVars2D(thetat_, thetab_, data_hist, *pdf_hist, "", constants::format);
        tools::PlotVars2D(thetat_, phit_, data_hist, *pdf_hist, "", constants::format);
        tools::PlotVars2D(thetab_, phit_, data_hist, *pdf_hist, "", constants::format);

        tools::PlotPull2D(thetat_, thetab_, data_hist, *pdf_hist, "", constants::format, true);
        tools::PlotPull2D(thetat_, phit_, data_hist, *pdf_hist, "", constants::format, true);
        tools::PlotPull2D(thetab_, phit_, data_hist, *pdf_hist, "", constants::format, true);
    }
    return json_results;
}
