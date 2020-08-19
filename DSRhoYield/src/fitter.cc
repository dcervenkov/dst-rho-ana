/**
 *  @file    fitter.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2015-07-31
 *
 *  @brief This class performs the yield fitting itself as well as plotting
 *
 */

#include "fitter.h"

// ROOT includes
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TMacro.h"
#include "TPaveText.h"
#include "TROOT.h"

// Local includes
#include "log.h"
#include "tools.h"

/**
 * Constructor that takes the output directory as argument. Plots as well
 * as ROOT files are saved to that directory.
 *
 * @param output_dir Output directory
 */
Fitter::Fitter(const char* output_dir) {
    dataset_vars_.add(de_);
    dataset_vars_.add(csbdtg_);
    dataset_vars_.add(evmcflag_);
    dataset_vars_.add(candsel_);

    dataset_vars_.add(thetat_);
    dataset_vars_.add(thetab_);
    dataset_vars_.add(phit_);
    dataset_vars_.add(dt_);

    dataset_vars_.add(vrusable_);
    dataset_vars_.add(vrvtxz_);
    dataset_vars_.add(vrerr6_);
    dataset_vars_.add(vrchi2_);
    dataset_vars_.add(vrndf_);
    dataset_vars_.add(vrntrk_);

    dataset_vars_.add(vtusable_);
    dataset_vars_.add(vtvtxz_);
    dataset_vars_.add(vterr6_);
    dataset_vars_.add(vtchi2_);
    dataset_vars_.add(vtndf_);
    dataset_vars_.add(vtntrk_);

    dataset_vars_.add(nocand_);

    Setup(Components::all);
    char output_filename[1000];
    snprintf(output_filename, 1000, "%s/fit_results.root", output_dir);
    output_file_ = new TFile(output_filename, "RECREATE");
    width_factor_.setConstant();
}

Fitter::~Fitter() {
    output_file_->Close();
    delete output_file_;
}

/**
 * Close the output file and write it to disk. Then reopen it in READ mode.
 */
void Fitter::CloseOutput() {
    // ReOpen instead of Close, because Close causes a segfault
    // ReOpen allows the TApplication windows to stay open, while finishing the file write
    output_file_->ReOpen("READ");
}

/**
 * Fit the model to data in the chain of files given.
 *
 * @param chain Chain of files containing data to be fitted
 */
void Fitter::FitTo(TTree* chain) {
    if (data_set_) delete data_set_;
    data_set_ = new RooDataSet("data_set", "data_set", chain, dataset_vars_, data_cut_);
    Log::print(Log::info, "Fitting %i events\n", data_set_->numEntries());

    //	fit_result_ = active_model_->fitTo(*data_set_, RooFit::Save(), RooFit::Minimizer("Minuit2"),
    //RooFit::Hesse(1), RooFit::Minos(1), RooFit::NumCPU(4));
    fit_result_ = active_model_->fitTo(*data_set_, RooFit::Save(), RooFit::Minimizer("Minuit2"),
                                       RooFit::NumCPU(numCPUs_));
}

/**
 * Setup the model and data cuts for a certain PDF component (signal, crossfeed, background, etc.).
 *
 * @param component Which component of the PDF is to be used for fitting
 */
void Fitter::Setup(Components component) {
    active_component_ = component;
    data_cut_ = tools::GetCommonCutsString();
    switch (component) {
        case Components::signal:
            data_cut_ += "&&(evmcflag==1||evmcflag==7||evmcflag==8)";
            active_model_ = &signal_model_;
            active_model_name_ = "signal";
            width_factor_.setConstant();
            break;

        case Components::crossfeed:
            data_cut_ += "&&(candsel!=0&&(evmcflag!=1&&evmcflag!=7))";
            active_model_ = &cross_model_;
            active_model_name_ = "crossfeed";
            width_factor_.setConstant();
            break;

        case Components::signal_plus_crossfeed:
            data_cut_ +=
                "&&((evmcflag==1||evmcflag==7||evmcflag==8)||(candsel!=0&&(evmcflag!=1&&evmcflag!=7)))";
            active_model_ = &signal_plus_cross_model_;
            active_model_name_ = "signal_plus_crossfeed";
            width_factor_.setConstant();
            break;

        case Components::background:
            data_cut_ +=
                "&&(!((evmcflag==1||evmcflag==7||evmcflag==8)||(candsel!=0&&(evmcflag!=1&&evmcflag!=7))))";
            active_model_ = &bkg_model_;
            active_model_name_ = "background";
            width_factor_.setConstant();
            break;

        case Components::all:
            active_model_ = &model_;
            active_model_name_ = "all";
            if (MC_) {
                width_factor_.setConstant(true);
            } else {
                width_factor_.setConstant(false);
            }
            break;

    }  // switch
}

/**
 * Fixes the PDF shape by setting relevant parameters constant.
 *
 * @param component PDF component shape that should be fixed
 */
void Fitter::FixShape(Components component) {
    switch (component) {
        case Components::signal:
            sig_gaus_mean_.setConstant();
            sig_gaus_sigma_l_.setConstant();
            sig_gaus_sigma_r_.setConstant();
            sig_cb_t0_.setConstant();
            sig_cb_sigma_.setConstant();
            sig_cb_cut_.setConstant();
            sig_cb_power_.setConstant();
            sig_f_.setConstant();
            break;

        case Components::crossfeed:
            cross_gaus_mean_.setConstant();
            cross_gaus_sigma_.setConstant();
            cross_tau_.setConstant();
            cross_f_.setConstant();
            break;

        case Components::signal_plus_crossfeed:
            signal_plus_cross_f_.setConstant();
            break;

        case Components::background:
            bkg_poly_p1_.setConstant();
            bkg_poly_p2_.setConstant();
            break;

        default:
            break;
    }
}

/**
 * Saves fit results to the output ROOT file in a form of a TMacro.
 */
void Fitter::WriteFitResults() {
    TMacro macro(active_model_name_, active_model_name_);
    const RooArgList results = fit_result_->floatParsFinal();
    char line[1000];
    snprintf(line, 1000, "*%s*", active_model_name_.Data());
    macro.AddLine(line);
    snprintf(line, 1000, "%i entries, fit status: %i", data_set_->numEntries(),
             fit_result_->status());
    macro.AddLine(line);
    for (int i = 0; i < results.getSize(); i++) {
        snprintf(line, 1000, "%-30s % 8f % 8f % 8f", results[i].GetName(),
                 static_cast<RooRealVar&>(results[i]).getVal(),
                 static_cast<RooRealVar&>(results[i]).getErrorLo(),
                 static_cast<RooRealVar&>(results[i]).getErrorHi());
        macro.AddLine(line);
        results[i].Write();
    }

    RooPlot* frame_de = de_.frame();
    data_set_->plotOn(frame_de);
    active_model_->plotOn(frame_de);
    snprintf(line, 1000, "chi^2 = %f\n", frame_de->chiSquare());
    macro.AddLine(line);
    snprintf(line, 1000, "rel://files/DSRhoYield/ver7/K/%s.png\n", active_model_name_.Data());
    macro.AddLine(line);

    macro.Write();

    fit_result_->Write();
}

/**
 * Plot the data, full model PDFs (or a component thereof if Setup was called) as well as signal and
 * background PDFs.
 */
void Fitter::Plot() {
    RooPlot* frame_de = de_.frame();
    data_set_->plotOn(frame_de, RooFit::XErrorSize(0.001));
    active_model_->plotOn(frame_de);
    const double chi2 = frame_de->chiSquare();
    Log::print(Log::info, "chi^2 = %f\n", chi2);

    // Create a new frame to draw the pull distribution and add the distribution to the frame
    RooPlot* frame_pull = de_.frame(RooFit::Title("Pull Distribution"));
    frame_pull->SetTitle("");
    RooHist* hpull = frame_de->pullHist();
    hpull->SetFillColor(kGray);
    // The only working way to get rid of error bars; HIST draw option doesn't work with RooPlot
    for (int i = 0; i < hpull->GetN(); i++) {
        hpull->SetPointError(i, 0, 0, 0, 0);
    }
    frame_pull->addPlotable(hpull, "B");
    // We plot again without bars, so the points are not half covered by bars as in case of "BP"
    // draw option. We need to create and plot a clone, because the ROOT object ownership is
    // transfered to the RooPlot by addPlotable(). If we just added the hpull twice, we would get a
    // segfault.
    RooHist* hpull_clone = static_cast<RooHist*>(hpull->Clone());
    frame_pull->addPlotable(
        hpull_clone,
        "P");  // We plot again without bars, so the points are not half covered by bars

    TLegend* legend = 0;
    // The plotting of the components has to take place after creating the pull histogram, otherwise
    // we would have to make sure the correct (full) pdf is used for the computation
    switch (active_component_) {
        case Components::signal:
            active_model_->plotOn(frame_de, RooFit::Components(sig_cb_), RooFit::LineColor(29),
                                  RooFit::LineStyle(kDashed));
            active_model_->plotOn(frame_de, RooFit::Components(sig_gaus_), RooFit::LineColor(32),
                                  RooFit::LineStyle(kDashed));
            break;

        case Components::crossfeed:
            active_model_->plotOn(frame_de, RooFit::Components(cross_exponential_),
                                  RooFit::LineColor(40), RooFit::LineStyle(kDashed));
            active_model_->plotOn(frame_de, RooFit::Components(cross_gaus_), RooFit::LineColor(49),
                                  RooFit::LineStyle(kDashed));
            break;

        case Components::signal_plus_crossfeed:
            active_model_->plotOn(frame_de, RooFit::Components(signal_model_), RooFit::LineColor(4),
                                  RooFit::LineStyle(kDashed));
            active_model_->plotOn(frame_de, RooFit::Components(cross_model_), RooFit::LineColor(7),
                                  RooFit::LineStyle(kDashed));
            break;

        case Components::background:
            break;

        case Components::all:
            active_model_->plotOn(frame_de, RooFit::Name("signal"),
                                  RooFit::Components(signal_model_), RooFit::LineColor(4),
                                  RooFit::LineStyle(kDashed));
            active_model_->plotOn(frame_de, RooFit::Name("crossfeed"),
                                  RooFit::Components(cross_model_), RooFit::LineColor(7),
                                  RooFit::LineStyle(kDashed));
            active_model_->plotOn(frame_de, RooFit::Name("background"),
                                  RooFit::Components(bkg_model_), RooFit::LineColor(5),
                                  RooFit::LineStyle(kDashed));

            legend = new TLegend(0.77, 0.70, 0.92, 0.9);
            legend->SetFillColor(kWhite);
            legend->SetLineColor(kWhite);
            legend->SetTextFont(43);
            legend->SetTextSize(18);
            legend->AddEntry(frame_de->findObject("signal"), "CR", "L");
            legend->AddEntry(frame_de->findObject("crossfeed"), "SCF", "L");
            legend->AddEntry(frame_de->findObject("background"), "BKG", "L");
            break;

        default:
            break;
    }

    // We plot the data points again, so they are drawn on top of the pdfs
    data_set_->plotOn(frame_de, RooFit::XErrorSize(0.001));

    TCanvas* canvas = new TCanvas(active_model_name_, active_model_name_, 500, 500);

    TPad* pad_de = new TPad("pad_de", "pad_de", 0, 0.25, 1, 1);
    TPad* pad_pull = new TPad("pad_pull", "pad_pull", 0, 0, 1, 0.25);

    pad_de->Draw();
    pad_pull->Draw();

    pad_de->cd();
    pad_de->SetBottomMargin(0.0);

    frame_de->SetTitle("");
    frame_de->GetXaxis()->SetTitle("");
    frame_de->GetXaxis()->SetLabelSize(0);
    // This line makes sure the 0 is not drawn as it would overlap with the lower pad
    frame_de->GetYaxis()->SetRangeUser(0.001, frame_de->GetMaximum());
    frame_de->GetYaxis()->SetTitle("");
    frame_de->Draw();

    TPaveText* stat_box = CreateStatBox(chi2);
    stat_box->Draw();

    if (legend) {
        legend->Draw();
    }

    pad_pull->cd();
    pad_pull->SetTopMargin(0.0);
    pad_pull->SetBottomMargin(0.35);

    // This is to make sure the tick length (which is normally proportional to the pad size)
    // is the same length in the lower (residual) plot as in the upper one
    frame_pull->GetXaxis()->SetTickLength(0.03 * pad_de->GetAbsHNDC() / pad_pull->GetAbsHNDC());
    frame_pull->GetXaxis()->SetTitle("#DeltaE [GeV]");
    frame_pull->GetXaxis()->SetTitleOffset(4.3);
    frame_pull->SetLabelOffset(0.01 * pad_de->GetAbsHNDC() / pad_pull->GetAbsHNDC());
    frame_pull->GetYaxis()->SetRangeUser(-5, 5);
    frame_pull->GetYaxis()->SetNdivisions(505);
    frame_pull->Draw();

    canvas->Write();
    canvas->SaveAs(".pdf");
    canvas->SaveAs(".png");
}

TPaveText* Fitter::CreateStatBox(double chi2) {
    const RooArgList results = fit_result_->floatParsFinal();
    const double x1 = 0.9 - results.getSize() * 0.06;
    TPaveText* stat_box = new TPaveText(0.28, x1, 0.33, 0.9, "NDC");
    stat_box->SetShadowColor(kWhite);
    stat_box->SetBorderSize(0);
    stat_box->SetFillColor(kWhite);
    stat_box->SetTextFont(43);
    stat_box->SetTextSize(14);
    stat_box->SetY1NDC(0.1);
    char line[1000];
    for (int i = 0; i < results.getSize(); i++) {
        snprintf(line, 1000, "%s = %.3f +- %.3f", results[i].GetTitle(),
                 static_cast<RooRealVar&>(results[i]).getVal(),
                 static_cast<RooRealVar&>(results[i]).getError());
        stat_box->AddText(line);
    }
    snprintf(line, 1000, "#chi^{2} = %.2f\n", chi2);
    stat_box->AddText(line);
    return stat_box;
}

void Fitter::SPlotFull(TChain* chain) {
    Log::print(Log::debug, "********** SPLOT FULL *********\n");
    RooStats::SPlot* sDataX = new RooStats::SPlot("sData", "An SPlot", *data_set_, &model_,
                                                  RooArgList(n_signal_plus_cross_model_, n_bkg_));

    RooDataSet* dataw_sig = new RooDataSet(data_set_->GetName(), data_set_->GetTitle(), data_set_,
                                           *data_set_->get(), 0, "n_signal_plus_cross_model_sw");
    RooDataSet* dataw_bg = new RooDataSet(data_set_->GetName(), data_set_->GetTitle(), data_set_,
                                          *data_set_->get(), 0, "n_bkg_sw");

    TCanvas* can = new TCanvas("splot", "splot", 500, 1000);

    RooPlot* frame_sig_thetab = thetab_.frame(100);
    dataw_sig->plotOn(frame_sig_thetab, RooFit::DataError(RooAbsData::SumW2));
    data_set_->plotOn(
        frame_sig_thetab, RooFit::MarkerColor(kRed),
        RooFit::Cut(
            "(evmcflag==1||evmcflag==7||evmcflag==8)||(candsel!=0&&(evmcflag!=1&&evmcflag!=7))"));

    RooPlot* frame_bkg_thetab = thetab_.frame(100);
    dataw_bg->plotOn(frame_bkg_thetab, RooFit::DataError(RooAbsData::SumW2));
    data_set_->plotOn(frame_bkg_thetab, RooFit::MarkerColor(kRed),
                      RooFit::Cut("!((evmcflag==1||evmcflag==7||evmcflag==8)||(candsel!=0&&("
                                  "evmcflag!=1&&evmcflag!=7)))"));

    can->Divide(1, 2);
    can->cd(1);
    frame_sig_thetab->Draw();
    can->cd(2);
    frame_bkg_thetab->Draw();
    can->Print(".pdf");

    TH2F* h2f = data_set_->createHistogram(de_, thetab_, 100, 100);
    TCanvas* can2 = new TCanvas("can2", "can2", 500, 500);
    h2f->Draw();
    Log::print(Log::debug, "**************Corr full = %f\n", h2f->GetCorrelationFactor());
    can2->Print("corr.pdf");

    delete sDataX;
}

void Fitter::SPlotSC(TChain* chain) {
    Log::print(Log::debug, "********** SPLOT SC *********\n");

    this->Setup(Components::signal_plus_crossfeed);
    data_set_ = new RooDataSet("data_set", "data_set", chain,
                               RooArgSet(de_, evmcflag_, candsel_, thetab_), data_cut_);

    RooRealVar n_signal_model_("n_signal_model", "n_{sig}", 20000, 1000, 100000);
    RooRealVar n_cross_model_("n_cross_model", "n_{cf}", 20000, 1000, 100000);
    RooAddPdf sc_model_("sc_model", "", RooArgList(signal_model_, cross_model_),
                        RooArgList(n_signal_model_, n_cross_model_));

    RooStats::SPlot* sDataX = new RooStats::SPlot("sData", "An SPlot", *data_set_, &sc_model_,
                                                  RooArgList(n_signal_model_, n_cross_model_));

    RooDataSet* dataw_sig = new RooDataSet(data_set_->GetName(), data_set_->GetTitle(), data_set_,
                                           *data_set_->get(), 0, "n_signal_model_sw");
    RooDataSet* dataw_bg = new RooDataSet(data_set_->GetName(), data_set_->GetTitle(), data_set_,
                                          *data_set_->get(), 0, "n_cross_model_sw");

    TCanvas* can = new TCanvas("splot_sc", "splot_sc", 500, 1000);

    RooPlot* frame_sig_thetab = thetab_.frame(100);
    dataw_sig->plotOn(frame_sig_thetab, RooFit::DataError(RooAbsData::SumW2));
    data_set_->plotOn(frame_sig_thetab, RooFit::MarkerColor(kRed),
                      RooFit::Cut("(evmcflag==1||evmcflag==7||evmcflag==8)"));

    RooPlot* frame_bkg_thetab = thetab_.frame(100);
    dataw_bg->plotOn(frame_bkg_thetab, RooFit::DataError(RooAbsData::SumW2));
    data_set_->plotOn(frame_bkg_thetab, RooFit::MarkerColor(kRed),
                      RooFit::Cut("(candsel!=0&&(evmcflag!=1&&evmcflag!=7))"));

    can->Divide(1, 2);
    can->cd(1);
    frame_sig_thetab->Draw();
    can->cd(2);
    frame_bkg_thetab->Draw();
    can->SaveAs(".pdf");

    delete sDataX;
}

void Fitter::SPlotSB(TChain* chain) {
    Log::print(Log::debug, "********** SPLOT SB *********\n");

    data_cut_ =
        "(evmcflag==1||evmcflag==7||evmcflag==8)||!((evmcflag==1||evmcflag==7||evmcflag==8)||("
        "candsel!=0&&(evmcflag!=1&&evmcflag!=7)))";
    data_set_ = new RooDataSet("data_set", "data_set", chain,
                               RooArgSet(de_, evmcflag_, candsel_, thetab_), data_cut_);

    RooRealVar n_signal_model_("n_signal_model", "n_{sig}", 20000, 1000, 100000);
    RooRealVar n_bkg_model_("n_bkg_model", "n_{bkg}", 20000, 1000, 100000);
    RooAddPdf sb_model_("sb_model", "", RooArgList(signal_model_, bkg_model_),
                        RooArgList(n_signal_model_, n_bkg_model_));

    RooStats::SPlot* sDataX = new RooStats::SPlot("sData", "An SPlot", *data_set_, &sb_model_,
                                                  RooArgList(n_signal_model_, n_bkg_model_));

    RooDataSet* dataw_sig = new RooDataSet(data_set_->GetName(), data_set_->GetTitle(), data_set_,
                                           *data_set_->get(), 0, "n_signal_model_sw");
    RooDataSet* dataw_bg = new RooDataSet(data_set_->GetName(), data_set_->GetTitle(), data_set_,
                                          *data_set_->get(), 0, "n_bkg_model_sw");

    TCanvas* can = new TCanvas("splot_sb", "splot_sb", 500, 1000);

    RooPlot* frame_sig_thetab = thetab_.frame(100);
    dataw_sig->plotOn(frame_sig_thetab, RooFit::DataError(RooAbsData::SumW2));
    data_set_->plotOn(frame_sig_thetab, RooFit::MarkerColor(kRed),
                      RooFit::Cut("(evmcflag==1||evmcflag==7||evmcflag==8)"));

    RooPlot* frame_bkg_thetab = thetab_.frame(100);
    dataw_bg->plotOn(frame_bkg_thetab, RooFit::DataError(RooAbsData::SumW2));
    data_set_->plotOn(frame_bkg_thetab, RooFit::MarkerColor(kRed),
                      RooFit::Cut("!((evmcflag==1||evmcflag==7||evmcflag==8)||(candsel!=0&&("
                                  "evmcflag!=1&&evmcflag!=7)))"));

    can->Divide(1, 2);
    can->cd(1);
    frame_sig_thetab->Draw();
    can->cd(2);
    frame_bkg_thetab->Draw();
    can->SaveAs(".pdf");

    delete sDataX;
}

void Fitter::SPlotCB(TChain* chain) {
    Log::print(Log::debug, "********** SPLOT CB *********\n");

    data_cut_ = "!(evmcflag==1||evmcflag==7||evmcflag==8)";
    data_set_ = new RooDataSet("data_set", "data_set", chain,
                               RooArgSet(de_, evmcflag_, candsel_, thetab_), data_cut_);

    RooRealVar n_cross_model_("n_cross_model", "n_{cf}", 20000, 1000, 100000);
    RooRealVar n_bkg_model_("n_bkg_model", "n_{bkg}", 20000, 1000, 100000);
    RooAddPdf cb_model_("cb_model", "", RooArgList(cross_model_, bkg_model_),
                        RooArgList(n_cross_model_, n_bkg_model_));

    RooStats::SPlot* sDataX = new RooStats::SPlot("sData", "An SPlot", *data_set_, &cb_model_,
                                                  RooArgList(n_cross_model_, n_bkg_model_));

    RooDataSet* dataw_sig = new RooDataSet(data_set_->GetName(), data_set_->GetTitle(), data_set_,
                                           *data_set_->get(), 0, "n_cross_model_sw");
    RooDataSet* dataw_bg = new RooDataSet(data_set_->GetName(), data_set_->GetTitle(), data_set_,
                                          *data_set_->get(), 0, "n_bkg_model_sw");

    TCanvas* can = new TCanvas("splot_cb", "splot_cb", 500, 1000);

    RooPlot* frame_sig_thetab = thetab_.frame(100);
    dataw_sig->plotOn(frame_sig_thetab, RooFit::DataError(RooAbsData::SumW2));
    data_set_->plotOn(frame_sig_thetab, RooFit::MarkerColor(kRed),
                      RooFit::Cut("(candsel!=0&&(evmcflag!=1&&evmcflag!=7))"));

    RooPlot* frame_bkg_thetab = thetab_.frame(100);
    dataw_bg->plotOn(frame_bkg_thetab, RooFit::DataError(RooAbsData::SumW2));
    data_set_->plotOn(frame_bkg_thetab, RooFit::MarkerColor(kRed),
                      RooFit::Cut("!(candsel!=0&&(evmcflag!=1&&evmcflag!=7))"));

    can->Divide(1, 2);
    can->cd(1);
    frame_sig_thetab->Draw();
    can->cd(2);
    frame_bkg_thetab->Draw();
    can->SaveAs(".pdf");

    delete sDataX;
}

double Fitter::GetCorrelation(TChain* chain, const RooRealVar& var1, const RooRealVar& var2,
                              bool save_plot) {
    if (data_set_) delete data_set_;
    data_set_ = new RooDataSet("data_set", "data_set", chain,
                               RooArgSet(de_, evmcflag_, candsel_, thetab_), data_cut_);
    data_set_->Print();
    TH2F* correlation_histo = data_set_->createHistogram(var1, var2, 100, 100, data_cut_);
    double correlation = correlation_histo->GetCorrelationFactor();

    if (save_plot) {
        TString title("corr_");
        title += active_model_name_;
        title += "_";
        title += var1.GetName();
        title += "_";
        title += var2.GetName();
        TCanvas c_correlation(title, title, 500, 500);
        c_correlation.SetTicks(1, 1);

        correlation_histo->SetStats(false);
        correlation_histo->GetXaxis()->SetTitle(var1.GetTitle());
        correlation_histo->SetTitleOffset(1.2);
        correlation_histo->GetYaxis()->SetTitle(var2.GetTitle());

        TPaveText* stat_box = new TPaveText(0.7, 0.91, 0.9, 0.96, "NDC");
        stat_box->SetShadowColor(kWhite);
        stat_box->SetBorderSize(0);
        stat_box->SetFillColor(kWhite);
        stat_box->SetTextFont(43);
        stat_box->SetTextSize(14);
        stat_box->SetY1NDC(0.1);
        char line[1000];
        snprintf(line, 1000, "f_{corr} = %.3f", correlation);
        stat_box->AddText(line);

        correlation_histo->SetTitle("");
        correlation_histo->Draw();
        stat_box->Draw();

        c_correlation.SaveAs(".pdf");
        c_correlation.SaveAs(".png");
    }

    delete correlation_histo;
    return correlation;
}

/**
 * Print the resulting number of events and fractions
 */
void Fitter::PrintResults() {
    printf("\nResults:\n");

    if (!MC_) {
        printf("F_w   = %f +- %f\n", width_factor_.getVal(),
               width_factor_.getPropagatedError(*fit_result_));
    }
    printf("f_cr  = %f +- %f\n", f_cr_.getVal(), f_cr_.getPropagatedError(*fit_result_));
    printf("f_scf = %f +- %f\n", f_scf_.getVal(), f_scf_.getPropagatedError(*fit_result_));
    printf("f_bkg = %f +- %f\n", f_bkg_.getVal(), f_bkg_.getPropagatedError(*fit_result_));

    printf("n_cr  = %.0f +- %.0f\n", n_cr_.getVal(), n_cr_.getPropagatedError(*fit_result_));
    printf("n_scf = %.0f +- %.0f\n", n_scf_.getVal(), n_scf_.getPropagatedError(*fit_result_));
    printf("n_bkg = %.0f +- %.0f\n", n_bkg_.getVal(), n_bkg_.getError());
}
