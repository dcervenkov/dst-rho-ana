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
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooProdPdf.h"
#include "RooRealVar.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TEnv.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TTree.h"

// Local includes
#include "colors.h"
#include "constants.h"

Fitter::Fitter(const char* evtgen_filepath, const char* gsim_filepath, const char* output_dir) {
    dec_type_.defineType("a", 1);
    dec_type_.defineType("ab", 2);
    dec_type_.defineType("b", 3);
    dec_type_.defineType("bb", 4);

    if (strstr(evtgen_filepath, ".root")) {
        TFile* evtgen_file = new TFile(evtgen_filepath);
        evtgen_dataset_ = dynamic_cast<RooDataSet*>(evtgen_file->Get("dataset"));
    } else {
        evtgen_dataset_ =
            RooDataSet::read(evtgen_filepath, RooArgList(thetat_, thetab_, phit_, dt_, dec_type_));
    }

    TFile* gsim_file = new TFile(gsim_filepath);
    TTree* gsim_tree = dynamic_cast<TTree*>(gsim_file->Get("h2000"));
    gsim_dataset_ = new RooDataSet("gsim_dataset", "gsim dataset", gsim_tree,
                                   RooArgSet(thetat_, thetab_, phit_, vrvtxz_, vtvtxz_, evmcflag_),
                                   "evmcflag==1");
    // Add calculated dt to the dataset
    gsim_dataset_->addColumn(dt_formula_);
    delete gsim_tree;
    gsim_file->Close();

    gEnv->SetValue("Canvas.PrintDirectory", output_dir);
    printf("print dir: %s\n", gEnv->GetValue("Canvas.PrintDirectory", "not found"));
    output_file_ = new TFile(TString(output_dir) + "/plots.root", "RECREATE");
}

Fitter::~Fitter() { output_file_->Close(); }

void Fitter::PlotVar(RooRealVar& var) {
    if (canvas_var_) delete canvas_var_;
    canvas_var_ = new TCanvas(TString(var.GetName()) + "_canvas",
                              TString(var.GetTitle()) + " canvas", 500, 500);
    if (plot_var_) delete plot_var_;
    plot_var_ = var.frame();

    evtgen_dataset_->plotOn(plot_var_);
    gsim_dataset_->plotOn(plot_var_, RooFit::MarkerColor(2), RooFit::LineColor(2));
    plot_var_->Draw();

    plot_var_->SetTitle("");
    plot_var_->GetXaxis()->SetTitle(TString(var.GetTitle()));
    plot_var_->GetYaxis()->SetTitle("");

    canvas_var_->Write();
    canvas_var_->SaveAs(format);
}

void Fitter::PlotVars2D(RooRealVar& var1, RooRealVar& var2) {
    if (canvas_var_) delete canvas_var_;
    canvas_var_ = new TCanvas(
        TString(var1.GetName()) + "_" + TString(var2.GetName()) + "_evtgen_canvas",
        TString(var1.GetTitle()) + "_" + TString(var2.GetName()) + " evtgen canvas", 500, 500);

    TH2D* evtgen_histo = static_cast<TH2D*>(
        evtgen_dataset_->createHistogram("evtgen_histo", var1, RooFit::YVar(var2)));

    canvas_var_->SetRightMargin(0.14);

    evtgen_histo->SetMinimum(0);
    evtgen_histo->SetTitle("");
    evtgen_histo->Draw("colz");
    evtgen_histo->GetZaxis()->SetTitle("");

    canvas_var_->Write();
    canvas_var_->SaveAs(format);

    TH2D* gsim_histo =
        static_cast<TH2D*>(gsim_dataset_->createHistogram("gsim_histo", var1, RooFit::YVar(var2)));

    gsim_histo->SetMinimum(0);
    gsim_histo->SetTitle("");
    gsim_histo->Draw("colz");
    gsim_histo->GetZaxis()->SetTitle("");

    canvas_var_->SetName(TString(var1.GetName()) + "_" + TString(var2.GetName()) + "_gsim_canvas");
    canvas_var_->SetTitle(TString(var1.GetName()) + "_" + TString(var2.GetName()) + " gsim canvas");

    canvas_var_->Write();
    canvas_var_->SaveAs(format);

    delete evtgen_histo;
    delete gsim_histo;
}

void Fitter::PlotEfficiency(RooRealVar& var, bool plot_model, bool legend_position_top,
                            bool legend_position_left) {
    if (canvas_eff_) delete canvas_eff_;
    canvas_eff_ = new TCanvas(TString(var.GetName()) + "_eff_canvas",
                              TString(var.GetTitle()) + " eff canvas", 500, 500);
    if (plot_eff_) delete plot_eff_;
    plot_eff_ = var.frame();

    TPad* pad_eff;
    TPad* pad_pull;
    if (plot_model) {
        pad_eff = new TPad("pad_eff", "pad_eff", 0, 0.25, 1, 1);
        pad_pull = new TPad("pad_pull", "pad_pull", 0, 0, 1, 0.25);
        pad_eff->Draw();
        pad_pull->Draw();

        pad_eff->cd();
        pad_eff->SetBottomMargin(0.0);
        pad_eff->SetLeftMargin(0.12);
    }

    if (var == thetat_) {
        model_ = thetat_model_e_;
        efficiency_ = &efficiency_thetat_;
    } else if (var == thetab_) {
        model_ = thetab_model_e_;
        efficiency_ = &efficiency_thetab_;
    } else if (var == phit_) {
        model_ = phit_model_e_;
        efficiency_ = &efficiency_phit_;
    }

    (*efficiency_)->plotOn(plot_eff_);

    if (plot_model) {
        model_->plotOn(plot_eff_);
        plot_eff_->GetXaxis()->SetTitle("");
        plot_eff_->GetXaxis()->SetLabelSize(0);
    }

    // This line makes sure the 0 is not drawn as it would overlap with the lower pad
    plot_eff_->GetYaxis()->SetRangeUser(0.001, plot_eff_->GetMaximum());
    plot_eff_->SetTitle("");
    plot_eff_->GetYaxis()->SetTitle("Efficiency");
    plot_eff_->GetYaxis()->SetTitleOffset(1.60);
    plot_eff_->Draw();

    if (plot_model) {
        const double chi2 = plot_eff_->chiSquare();
        TPaveText* stat_box = CreateStatBox(chi2, legend_position_top, legend_position_left);
        stat_box->Draw();

        pad_pull->cd();
        pad_pull->SetTopMargin(0.0);
        pad_pull->SetBottomMargin(0.35);
        pad_pull->SetLeftMargin(0.12);

        // Create a new frame to draw the pull distribution and add the distribution to the frame
        RooPlot* plot_eff_pull_ = var.frame(RooFit::Title("Pull Distribution"));
        plot_eff_pull_->SetTitle("");
        RooHist* hpull = plot_eff_->pullHist();
        hpull->SetFillColor(kGray);
        // The only working way to get rid of error bars; HIST draw option doesn't work with RooPlot
        for (int i = 0; i < hpull->GetN(); i++) {
            hpull->SetPointError(i, 0, 0, 0, 0);
        }
        plot_eff_pull_->addPlotable(hpull, "B");
        // We plot again without bars, so the points are not half covered by
        // bars as in case of "BP" draw option. We need to create and plot a
        // clone, because the ROOT object ownership is transfered to the RooPlot
        // by addPlotable().
        // If we just added the hpull twice, we would get a segfault.
		RooHist* hpull_clone = dynamic_cast<RooHist*>(hpull->Clone());
		// We plot again without bars, so the points are not half covered by bars
        plot_eff_pull_->addPlotable(hpull_clone, "P"); 

        plot_eff_pull_->GetXaxis()->SetTickLength(0.03 * pad_eff->GetAbsHNDC() /
                                                  pad_pull->GetAbsHNDC());
        plot_eff_pull_->GetXaxis()->SetTitle(TString(var.GetTitle()));
        plot_eff_pull_->GetXaxis()->SetTitleOffset(4.0);
        plot_eff_pull_->GetXaxis()->SetLabelOffset(0.01 * pad_eff->GetAbsHNDC() /
                                                   pad_pull->GetAbsHNDC());
        plot_eff_pull_->GetYaxis()->SetRangeUser(-5, 5);
        plot_eff_pull_->GetYaxis()->SetNdivisions(505);
        plot_eff_pull_->Draw();
    }

    canvas_eff_->Write();
    canvas_eff_->SaveAs(format);
}

void Fitter::PlotEfficiency2D(RooRealVar& var1, RooRealVar& var2) {
    if (canvas_eff_) delete canvas_eff_;
    canvas_eff_ = new TCanvas(
        TString(var1.GetName()) + "_" + TString(var2.GetName()) + "_eff_canvas",
        TString(var1.GetTitle()) + "_" + TString(var2.GetName()) + " eff canvas", 500, 500);

    TH2D* gsim_histo = static_cast<TH2D*>(
		gsim_dataset_->createHistogram("gsim_histo", var1, RooFit::YVar(var2)));
    TH2D* evtgen_histo = static_cast<TH2D*>(
        evtgen_dataset_->createHistogram("evtgen_histo", var1, RooFit::YVar(var2)));

    TH2* eff_histo = dynamic_cast<TH2*>(gsim_histo->Clone("eff_histo"));
    eff_histo->Divide(eff_histo, evtgen_histo, 1.0, 1.0, "B");

    // Empty (whiten) bins with too few entries
    for (int i = 1; i <= eff_histo->GetNbinsX(); i++) {
        for (int j = 1; j <= eff_histo->GetNbinsY(); j++) {
            if (gsim_histo->GetBinContent(i, j) == 0) {
                eff_histo->SetBinContent(i, j, -1);
            }
        }
    }

    printf("Correlation %s : %s is %f\n", var1.GetName(), var2.GetName(),
           eff_histo->GetCorrelationFactor());
    
    // This clone must happen before Draw() of eff_histo, otherwise problems
    // with eff_pull_histo color legend range occur.
    TH2* eff_pull_histo = dynamic_cast<TH2*>(eff_histo->Clone("eff_pull_histo"));

    TString name = "eff_";
    name += var1.GetName();
    name += "_";
    name += var2.GetName();

    canvas_eff_->SetRightMargin(0.14);

    eff_histo->SetTitle("");
    eff_histo->SetMinimum(0);
    eff_histo->SetMaximum(0.35);
    eff_histo->Draw("colz");
    eff_histo->GetZaxis()->SetTitle("");

    TPaveText* stat_box = new TPaveText(0.2, 0.94, 0.2, 1, "NDC");
    stat_box->SetShadowColor(kWhite);
    stat_box->SetBorderSize(0);
    stat_box->SetFillColor(kWhite);
    stat_box->SetTextFont(43);
    stat_box->SetTextSize(14);
    stat_box->SetY1NDC(0.1);

    char line[1000];
    snprintf(line, 1000, "c_{f} = %.3f\n", eff_histo->GetCorrelationFactor());
    stat_box->AddText(line);
    stat_box->Draw();

    canvas_eff_->Write();
    canvas_eff_->SaveAs(format);

    RooProdPdf efficiency_3D("efficiency_3D", "efficiency_3D",
                             RooArgList(*thetat_model_e_, *thetab_model_e_, *phit_model_e_));
    RooDataHist* efficiency_pdf_binned = efficiency_3D.generateBinned(
        RooArgSet(thetat_, thetab_, phit_), eff_histo->Integral(), kTRUE);
    TH2D* pdf_eff_histo = static_cast<TH2D*>(
        efficiency_pdf_binned->createHistogram("pdf_eff_histo", var1, RooFit::YVar(var2)));

    name = "pdf_eff_";
    name += var1.GetName();
    name += "_";
    name += var2.GetName();

    pdf_eff_histo->SetTitle("");
    pdf_eff_histo->SetMinimum(0);
    pdf_eff_histo->SetMaximum(0.35);
    pdf_eff_histo->Draw("colz");
    pdf_eff_histo->GetZaxis()->SetTitle("");

    canvas_eff_->SetName(TString(var1.GetName()) + "_" + TString(var2.GetName()) +
                         "_pdf_eff_canvas");
    canvas_eff_->SetTitle(TString(var1.GetTitle()) + "_" + TString(var2.GetName()) +
                          " PDF eff canvas");

    canvas_eff_->Write();
    canvas_eff_->SaveAs(format);

    eff_pull_histo->Add(pdf_eff_histo, -1);

    TH2* eff_errors_histo = dynamic_cast<TH2*>(eff_histo->Clone("eff_errors_histo"));
    eff_errors_histo->Reset();
    for (int i = 1; i <= eff_histo->GetNbinsX(); i++) {
        for (int j = 1; j <= eff_histo->GetNbinsY(); j++) {
            eff_errors_histo->SetBinContent(i, j, gsim_histo->GetBinError(i, j));
        }
    }

    eff_pull_histo->Divide(eff_errors_histo);

    // Whiten bins with too few data entries
    for (int i = 1; i <= eff_histo->GetNbinsX(); i++) {
        for (int j = 1; j <= eff_histo->GetNbinsY(); j++) {
            if (eff_histo->GetBinContent(i, j) == 0) {
                eff_pull_histo->SetBinContent(i, j, -1000);
            }
        }
    }

    gStyle->SetPalette(kLightTemperature);
    eff_pull_histo->SetTitle("");
    eff_pull_histo->GetZaxis()->SetTitle();
    eff_pull_histo->SetMinimum(-0.04);
    eff_pull_histo->SetMaximum(0.04);
    eff_pull_histo->SetOption("colz");
    eff_pull_histo->Draw();
    eff_pull_histo->Write();

    canvas_eff_->SetName(TString(var1.GetName()) + "_" + TString(var2.GetName()) +
                         "_eff_pull_canvas");
    canvas_eff_->SetTitle(TString(var1.GetTitle()) + "_" + TString(var2.GetName()) +
                          " eff pull canvas");

    canvas_eff_->Write();
    canvas_eff_->SaveAs(format);
    colors::setColors();

    delete eff_histo;
    delete gsim_histo;
    delete evtgen_histo;
    delete stat_box;
    delete eff_pull_histo;
    delete eff_errors_histo;
}

void Fitter::FitEfficiency(RooRealVar& var) {
    TH1* gsim_histo = gsim_dataset_->createHistogram("gsim_histo", var);
    TH1* evtgen_histo = evtgen_dataset_->createHistogram("evtgen_histo", var);

    TH1* eff_histo = dynamic_cast<TH1*>(gsim_histo->Clone("eff_histo"));
    eff_histo->Divide(eff_histo, evtgen_histo, 1.0, 1.0, "B");

    TString name = "eff_";
    name += var.GetName();

    if (var == thetat_) {
        efficiency_ = &efficiency_thetat_;
        model_ = thetat_model_e_;
    } else if (var == thetab_) {
        efficiency_ = &efficiency_thetab_;
        model_ = thetab_model_e_;
    } else if (var == phit_) {
        efficiency_ = &efficiency_phit_;
        model_ = phit_model_e_;
    }

    RooDataHist* eff_hist = new RooDataHist(name, name, var, eff_histo);
    *efficiency_ = eff_hist;

    result_ = model_->fitTo(*eff_hist, RooFit::Save(), RooFit::Minimizer("Minuit2"),
                            RooFit::SumW2Error(kTRUE), RooFit::NumCPU(4));

    delete eff_histo;
    delete gsim_histo;
    delete evtgen_histo;
}

TPaveText* Fitter::CreateStatBox(double chi2, bool position_top, bool position_left) {
    const RooArgList results = result_->floatParsFinal();
    double x_left, x_right, y_bottom, y_top;
    const double line_height = 0.06;

    if (position_top) {
        y_top = 0.9;
        y_bottom = y_top - results.getSize() * line_height;
    } else {
        y_bottom = 0.023;
        y_top = y_bottom + results.getSize() * line_height;
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
    snprintf(line, 1000, "#chi^{2} = %.2f\n", chi2);
    stat_box->AddText(line);
    return stat_box;
}

void Fitter::SetEfficiencyModel(const int model_num) {
    switch (model_num) {
        case 1:
            thetat_model_e_ = &thetat_model1_e_;
            thetab_model_e_ = &thetab_model1_e_;
            phit_model_e_ = &phit_model1_e_;
            break;
        case 2:
            thetat_model_e_ = &thetat_model2_e_;
            thetab_model_e_ = &thetab_model2_e_;
            phit_model_e_ = &phit_model2_e_;
            break;
        case 3:
            thetat_model_e_ = &thetat_model3_e_;
            thetab_model_e_ = &thetab_model3_e_;
            phit_model_e_ = &phit_model3_e_;
            break;
        case 4:
            thetat_model_e_ = &thetat_model4_e_;
            thetab_model_e_ = &thetab_model3_e_;
            phit_model_e_ = &phit_model2_e_;
            break;
    }
}
