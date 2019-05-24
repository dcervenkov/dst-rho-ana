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
#include "RooHistPdf.h"
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooTreeDataStore.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TEnv.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
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
#include "UniformDensity.hh"
#include "RooMeerkatPdf.hh"

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
    csbdtg_ = new RooRealVar("csbdtg", "csbdtg", -1, 1);

    shcosthb_ = new RooRealVar("shcosthb", "shcosthb", -1, 1);

    decaytype_ = new RooCategory("decaytype", "decaytype");
    decaytype_->defineType("a", 1);
    decaytype_->defineType("ab", 2);
    decaytype_->defineType("b", 3);
    decaytype_->defineType("bb", 4);

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

/**
 * Make a plot of a variable with both data points and pdf projection and display a pull plot
 * beneath it. Then save to a file.
 *
 * @param var Variable to be plotted
 * @param data Dataset against which to plot
 * @param pdf PDF to use for plotting and pull calculation
 * @param components [optional] Component PDFs to plot alongside the full PDF
 * @param title [optional] Title for the y-axis
 */
void FitterBKG::PlotWithPull(const RooRealVar& var, const RooDataSet& data, const RooAbsPdf& pdf,
                             const std::vector<RooAbsPdf*> components, const char* title) const {
    TString name = pdf.GetName();
    name += "_";
    name += data.GetName();
    name += "_";
    name += var.GetName();
    TCanvas canvas(name, name, 500, 500);

    RooPlot* plot = var.frame();
    TPad* pad_var;
    TPad* pad_pull;
    pad_var = new TPad("pad_var", "pad_var", 0, 0.25, 1, 1);
    pad_pull = new TPad("pad_pull", "pad_pull", 0, 0, 1, 0.25);
    pad_var->Draw();
    pad_pull->Draw();

    pad_var->cd();
    pad_var->SetBottomMargin(0.0);
    pad_var->SetLeftMargin(0.12);

    data.plotOn(plot);

    // Plot components before the total PDF as the pull plots are made from the
    // last plotted PDF
    const int colors[] = {4, 7, 5};
    int i = 0;
    for (auto component : components) {
        pdf.plotOn(plot, RooFit::ProjWData(conditional_vars_argset_, data, kFALSE),
                   RooFit::NumCPU(num_CPUs_), RooFit::LineColor(colors[i++]),
                   RooFit::Components(*component));
    }

    pdf.plotOn(plot, RooFit::ProjWData(conditional_vars_argset_, data, kFALSE),
               RooFit::NumCPU(num_CPUs_));

    plot->GetXaxis()->SetTitle("");
    plot->GetXaxis()->SetLabelSize(0);

    // This line makes sure the 0 is not drawn as it would overlap with the lower pad
    plot->GetYaxis()->SetRangeUser(0.001, plot->GetMaximum());
    plot->SetTitle("");
    plot->GetYaxis()->SetTitle(title);
    plot->GetYaxis()->SetTitleOffset(1.60);
    plot->Draw();

    // This is not an OK way of calculating ndof - e.g., in a case where the PDF
    // is defined only on a subset of the var range, or there are empty bins,
    // etc. it will be wrong. However, RooFit doesn't give access to anything
    // like ndof itself, so this will have to do for now.
    const int num_floating_pars = result_ ? result_->floatParsFinal().getSize() : 0;
    const int ndof = var.getBinning().numBins() - num_floating_pars;
    const double chi2 = plot->chiSquare(num_floating_pars) * ndof;
    TPaveText* stat_box = CreateStatBox(chi2, ndof, true, true);
    if (stat_box) {
        stat_box->Draw();
    }

    pad_pull->cd();
    pad_pull->SetTopMargin(0.0);
    pad_pull->SetBottomMargin(0.35);
    pad_pull->SetLeftMargin(0.12);

    // Create a new frame to draw the pull distribution and add the distribution to the frame
    RooPlot* plot_pull_ =
        var.frame(RooFit::Title("Pull Distribution"));
    plot_pull_->SetTitle("");
    RooHist* hpull = plot->pullHist();
    hpull->SetFillColor(kGray);
    // The only working way to get rid of error bars; HIST draw option doesn't work with RooPlot
    for (int i = 0; i < hpull->GetN(); i++) {
        hpull->SetPointError(i, 0, 0, 0, 0);
    }
    plot_pull_->addPlotable(hpull, "B");
    // We plot again without bars, so the points are not half covered by bars as in case of "BP"
    // draw option.
    // We need to create and plot a clone, because the ROOT object ownership is transfered to the
    // RooPlot by addPlotable().
    // If we just added the hpull twice, we would get a segfault.
    RooHist* hpull_clone = dynamic_cast<RooHist*>(hpull->Clone());

    // We plot again without bars, so the points are not half covered by bars
    plot_pull_->addPlotable(hpull_clone, "P");

    plot_pull_->GetXaxis()->SetTickLength(0.03 * pad_var->GetAbsHNDC() / pad_pull->GetAbsHNDC());
    plot_pull_->GetXaxis()->SetTitle(TString(var.GetTitle()));
    plot_pull_->GetXaxis()->SetTitleOffset(4.0);
    plot_pull_->GetXaxis()->SetLabelOffset(0.01 * pad_var->GetAbsHNDC() / pad_pull->GetAbsHNDC());
    plot_pull_->GetYaxis()->SetRangeUser(-5, 5);
    plot_pull_->GetYaxis()->SetNdivisions(505);
    plot_pull_->Draw();

    canvas.Write();
    canvas.SaveAs(constants::format);
}

/**
 * Create and return statbox with fit results to overlay on plots. Returns a nullptr pointer
 * if no fit result exist.
 *
 * @param chi2 Reduced chi2 of the projection plot
 * @param ndof Number of degrees of freedom (used for p-value computation)
 * @param position_top Should the box be displayed at the top or bottom of the plot
 * @param position_left Should the box be displayed at the left or right of the plot
 */
TPaveText* FitterBKG::CreateStatBox(const double chi2, const int ndof, const bool position_top,
                                    const bool position_left) const {
    RooArgList results = result_ ? result_->floatParsFinal() : RooArgList();

    double x_left, x_right, y_bottom, y_top;
    const double line_height = 0.06;

    if (position_top) {
        y_top = 0.9;
        y_bottom = y_top - (results.getSize() + 2) * line_height;
    } else {
        y_bottom = 0.023;
        y_top = y_bottom + (results.getSize() + 2) * line_height;
    }

    if (position_left) {
        x_left = 0.30;
        x_right = 0.30;
    } else {
        x_left = 0.7;
        x_right = 0.7;
    }

    TPaveText* stat_box = new TPaveText(x_left, y_bottom, x_right, y_top, "NDC");
    stat_box->SetShadowColor(kWhite);
    stat_box->SetBorderSize(0);
    stat_box->SetFillColor(kWhite);
    stat_box->SetTextFont(43);
    stat_box->SetTextSize(14);
    stat_box->SetY1NDC(0.1);

    char line[1000];
    for (int i = 0; i < results.getSize(); i++) {
        snprintf(line, 1000, "%s = %.3f +- %.3f", results[i].GetTitle(),
                 dynamic_cast<RooRealVar&>(results[i]).getVal(),
                 dynamic_cast<RooRealVar&>(results[i]).getError());
        stat_box->AddText(line);
    }
    snprintf(line, 1000, "#chi^{2}/ndof = %.2f (%.1f/%i)\n", chi2 / ndof, chi2, ndof);
    stat_box->AddText(line);
    snprintf(line, 1000, "p = %.2f\n", TMath::Prob(chi2, ndof));
    stat_box->AddText(line);
    return stat_box;
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

    TString a_cuts;
    TString ab_cuts;
    TString b_cuts;
    TString bb_cuts;
    a_cuts = "brecflav==1&&tagqr<0";
    ab_cuts = "brecflav==-1&&tagqr>0";
    b_cuts = "brecflav==1&&tagqr>0";
    bb_cuts = "brecflav==-1&&tagqr<0";

    Log::print(Log::info, "Reading %i input files...\n", file_names.size());
    dataset_ = new RooDataSet("dataset", "dataset", input_tree, dataset_vars_argset_, common_cuts);

    Log::print(Log::info, "Creating a, ab, b, bb labels...\n");
    // We add an identifying label to each of the 4 categories and then combine it into a single
    // dataset for RooSimultaneous fitting
    dataset_a_ = static_cast<RooDataSet*>(dataset_->reduce(a_cuts));
    decaytype_->setLabel("a");
    dataset_a_->addColumn(*decaytype_);
    dataset_a_->SetName("dataset_a");

    dataset_ab_ = static_cast<RooDataSet*>(dataset_->reduce(ab_cuts));
    decaytype_->setLabel("ab");
    dataset_ab_->addColumn(*decaytype_);
    dataset_ab_->SetName("dataset_ab");

    dataset_b_ = static_cast<RooDataSet*>(dataset_->reduce(b_cuts));
    decaytype_->setLabel("b");
    dataset_b_->addColumn(*decaytype_);
    dataset_b_->SetName("dataset_b");

    dataset_bb_ = static_cast<RooDataSet*>(dataset_->reduce(bb_cuts));
    decaytype_->setLabel("bb");
    dataset_bb_->addColumn(*decaytype_);
    dataset_bb_->SetName("dataset_bb");

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
    result_ = pdf->fitTo(*data, RooFit::Save(), RooFit::Minimizer("Minuit2"), RooFit::NumCPU(num_CPUs_));
}

void FitterBKG::CreateHistoPDF(RooDataSet* data) {

    // TH3F* scf_histo = Create3DHisto(data);
    // TH3F* scf_histo_normalized = NormalizePDF(scf_histo, 0, 1);

    // const char* output_file = "scf_histo";
    // TFile f(efficiency_file, "RECREATE");
    // binned_pdf->Write();
    // f.Close();

    RooDataHist data_hist("data_hist", "data_hist", RooArgSet(thetat_, thetab_, phit_), *data);
    RooHistPdf* scf_hist_pdf = new RooHistPdf("scf_hist_pdf", "scf_hist_pdf",
                                              RooArgSet(thetat_, thetab_, phit_), data_hist);

    scf_hist_pdf->Write();
}

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
                    // printf("bin %i: value = %f, interpolation = %f\n", bin, pdf->GetBinContent(bin), interpolation);
                    normalized_pdf->SetBinContent(bin, interpolation);
                }
            }
        }
    }

    Log::print(Log::info, "Fixed %i/%i (%.2f%%) bins.\n", fixed_bins, total_bins, (double)fixed_bins/total_bins * 100);
    return normalized_pdf;
}

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
    return new_value/points;
}


AdaptiveKernelDensity FitterBKG::FitKDE(RooDataSet* data) {
    OneDimPhaseSpace* phasespace_thetat = new OneDimPhaseSpace("phasespace_thetat", thetat_.getMin(), thetat_.getMax());
    OneDimPhaseSpace* phasespace_thetab = new OneDimPhaseSpace("phasespace_thetab", thetab_.getMin(), thetab_.getMax());
    OneDimPhaseSpace* phasespace_phit = new OneDimPhaseSpace("phasespace_phit", phit_.getMin(), phit_.getMax());
    CombinedPhaseSpace* phasespace = new CombinedPhaseSpace("phasespace", phasespace_thetat, phasespace_thetab,
                                  phasespace_phit);
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
    //                             bin_kde_pars[3], bin_kde_pars[4], bin_kde_pars[5], &formula_density);

    AdaptiveKernelDensity kde("KernelPDF", phasespace, eff_tree, "thetat", "thetab", "phit",
                              ada_kde_pars[0], ada_kde_pars[1], ada_kde_pars[2],
                              ada_kde_pars[3], ada_kde_pars[4], ada_kde_pars[5], &formula_density);
                            //   ada_kde_pars[3], ada_kde_pars[4], ada_kde_pars[5], &bin_kde);

    const char* output_file = "scf_kde";
    kde.writeToTextFile(output_file);

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
    const double scale = data->GetEntries()/model->GetEntries();

    RooDataHist roo_model("roo_model", "roo_model", RooArgList(thetat_, thetab_, phit_), model, scale);
    RooDataHist roo_data("roo_data", "roo_data", RooArgList(thetat_, thetab_, phit_), data);

    RooRealVar vars[] = {thetat_, thetab_, phit_};

    RooArgList list(thetat_, thetab_, phit_);
    RooMeerkatPdf meerkat_pdf("meerkat_pdf", "meerkat_pdf", list, &kde);
    RooHistPdf meerkat_histpdf("meerkat_histpdf", "meerkat_histpdf", RooArgSet(thetat_, thetab_, phit_), roo_model);
    for (auto var : vars) {
        // PlotWithPull(var, *dataset_, meerkat_pdf);
        PlotWithPull(var, *dataset_, meerkat_histpdf);
    }

    for (int i = 0; i < 3; i++) {
        for (int j = i + 1; j < 3; j++) {
            tools::PlotVars2D(vars[i], vars[j], roo_model);
            tools::PlotVars2D(vars[i], vars[j], roo_data);
            tools::PlotPull2D(vars[i], vars[j], roo_model, roo_data);
        }
    }

}
