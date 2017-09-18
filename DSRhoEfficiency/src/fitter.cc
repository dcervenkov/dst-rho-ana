/**
 *  @file    fitter.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2015-07-31
 *
 *  @brief This class performs the yield fitting itself as well as plotting
 *
 */

#include "fitter.h"

// Standard includes
#include <algorithm>
#include <cmath>

// ROOT includes
#include "RVersion.h"
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
#include "TEfficiency.h"
#include "TEnv.h"
#include "TFile.h"
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

void Fitter::EnlargeVarRanges(const double margin) {
    for (auto&& var : vars_) {
        const double min = var->getMin();
        const double max = var->getMax();
        const double range = max - min;
        var->setMin(min - range * margin);
        var->setMax(max + range * margin);
    }
}

void Fitter::MirrorDataAtEdges(RooDataSet* data) {
    int num_added = 0;
    const int num_entries = data->numEntries();
    for (int evt = 0; evt < num_entries; evt++) {
        // for (int evt = 0; evt < 1000; evt++) {
        // if (evt % 10000 == 0) {
        //     printf("Mirror processed %i events of %i\n", evt, num_entries);
        // }
        const RooArgSet* row = data->get(evt);
        double vals[3] = {row->getRealValue("thetat"), row->getRealValue("thetab"),
                          row->getRealValue("phit")};
        // delete row;

        for (int var_num1 = 0; var_num1 < 3; var_num1++) {
            if (CloseToEdge(vals, var_num1)) {
                double* mirror_vals = GetMirror(vals, var_num1);
                *vars_[0] = mirror_vals[0];
                *vars_[1] = mirror_vals[1];
                *vars_[2] = mirror_vals[2];
                // printf("1 org: %f %f %f, new: %f, %f, %f\n", vals[0], vals[1], vals[2],
                //        mirror_vals[0], mirror_vals[1], mirror_vals[2]);
                delete[] mirror_vals;
                data->add(RooArgSet(*vars_[0], *vars_[1], *vars_[2]));
                num_added++;
                for (int var_num2 = var_num1 + 1; var_num2 < 3; var_num2++) {
                    if (CloseToEdge(vals, var_num2)) {
                        double* mirror_vals = GetMirror(vals, var_num1, var_num2);
                        *vars_[0] = mirror_vals[0];
                        *vars_[1] = mirror_vals[1];
                        *vars_[2] = mirror_vals[2];
                        // printf("2 org: %f %f %f, new: %f, %f, %f\n", vals[0], vals[1], vals[2],
                        //        mirror_vals[0], mirror_vals[1], mirror_vals[2]);
                        delete[] mirror_vals;
                        data->add(RooArgSet(*vars_[0], *vars_[1], *vars_[2]));
                        num_added++;
                        for (int var_num3 = var_num2 + 1; var_num3 < 3; var_num3++) {
                            if (CloseToEdge(vals, var_num3)) {
                                double* mirror_vals = GetMirror(vals, var_num1, var_num2, var_num3);
                                *vars_[0] = mirror_vals[0];
                                *vars_[1] = mirror_vals[1];
                                *vars_[2] = mirror_vals[2];
                                // printf("3 org: %f %f %f, new: %f, %f, %f\n", vals[0], vals[1],
                                //        vals[2], mirror_vals[0], mirror_vals[1], mirror_vals[2]);
                                delete[] mirror_vals;
                                data->add(RooArgSet(*vars_[0], *vars_[1], *vars_[2]));
                                num_added++;
                            }
                        }
                    }
                }
            }
        }
    }
    printf("DBG: Added %i mirrored datapoints.\n", num_added);
}

double* Fitter::GetMirror(double vals[], const int var_num) {
    const double min = orig_vars_[var_num].getMin();
    const double max = orig_vars_[var_num].getMax();
    double* new_vals = new double[3];
    new_vals[0] = vals[0];
    new_vals[1] = vals[1];
    new_vals[2] = vals[2];
    switch (CloseToEdge(vals, var_num)) {
        case 1:
            new_vals[var_num] = min - (vals[var_num] - min);
            break;
        case 2:
            new_vals[var_num] = max + (max - vals[var_num]);
            break;
    }
    return new_vals;
}

double* Fitter::GetMirror(double vals[], const int var_num1, const int var_num2) {
    double* new_vals1 = GetMirror(vals, var_num1);
    double* new_vals2 = GetMirror(new_vals1, var_num2);
    delete[] new_vals1;
    return new_vals2;
}

double* Fitter::GetMirror(double vals[], const int var_num1, const int var_num2,
                          const int var_num3) {
    double* new_vals1 = GetMirror(vals, var_num1);
    double* new_vals2 = GetMirror(new_vals1, var_num2);
    double* new_vals3 = GetMirror(new_vals2, var_num3);
    delete[] new_vals1;
    delete[] new_vals2;
    return new_vals3;
}

int Fitter::CloseToEdge(double vals[], const int var_num) {
    const double margin = 0.1;
    const double min = orig_vars_[var_num].getMin();
    const double max = orig_vars_[var_num].getMax();
    const double range = max - min;
    if (vals[var_num] < min + range * margin) {
        return 1;
    } else if (vals[var_num] > max - range * margin) {
        return 2;
    }
    return 0;
}

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

void Fitter::PlotVar(RooRealVar& var, RooDataHist& data1, RooDataHist& data2) {
    if (canvas_var_) delete canvas_var_;
    canvas_var_ = new TCanvas(
        TString(var.GetName()) + "_" + TString(data1.GetName()) + "_" + TString(data2.GetName()),
        TString(var.GetTitle()) + " canvas", 500, 500);
    if (plot_var_) delete plot_var_;
    plot_var_ = var.frame();

    TPad* pad_var;
    TPad* pad_pull;

    pad_var = new TPad("pad_var", "pad_var", 0, 0.25, 1, 1);
    pad_pull = new TPad("pad_pull", "pad_pull", 0, 0, 1, 0.25);
    pad_var->Draw();
    pad_pull->Draw();

    pad_var->cd();
    pad_var->SetBottomMargin(0.0);
    pad_var->SetLeftMargin(0.12);

    data1.plotOn(plot_var_, RooFit::DataError(RooAbsData::None));
    data2.plotOn(plot_var_, RooFit::MarkerColor(2), RooFit::LineColor(2), RooFit::DataError(RooAbsData::None));

    plot_var_->SetTitle("");
    plot_var_->GetYaxis()->SetTitle("");
    plot_var_->GetXaxis()->SetTitle("");
    plot_var_->GetXaxis()->SetLabelSize(0);

    // This line makes sure the 0 is not drawn as it would overlap with the lower
    // pad
    plot_var_->GetYaxis()->SetRangeUser(0.001, plot_var_->GetMaximum());
    plot_var_->GetYaxis()->SetTitleOffset(1.60);
    plot_var_->Draw();

    pad_pull->cd();
    pad_pull->SetTopMargin(0.0);
    pad_pull->SetBottomMargin(0.35);
    pad_pull->SetLeftMargin(0.12);

    // Create a new frame to draw the pull distribution and add the distribution
    // to the frame
    RooPlot* plot_var_pull_ = var.frame(RooFit::Title("Pull Distribution"));
    plot_var_pull_->SetTitle("");

    TH1* data1_histo = data1.createHistogram("data1", var);
    TH1* data2_histo = data2.createHistogram("data2", var);
    TH1* pull_histo = dynamic_cast<TH1*>(data1_histo->Clone("pull_histo"));

    double p1;
    double p2;
    for (int i = 1; i <= pull_histo->GetNbinsX(); i++) {
        p1 = data1_histo->GetBinContent(i);
        p2 = data2_histo->GetBinContent(i);
        // pull_histo->SetBinContent(i, p1 / p2);
        pull_histo->SetBinContent(i, p1 - p2);
        // pull_histo->SetBinContent(i, (p1 - p2) / std::sqrt((p1 + p2) / 2));
    }
    RooHist* hpull = new RooHist(*pull_histo);

    hpull->SetFillColor(kGray);
    // The only working way to get rid of error bars; HIST draw option doesn't
    // work with RooPlot
    for (int i = 0; i < hpull->GetN(); i++) {
        hpull->SetPointError(i, 0, 0, 0, 0);
    }
    plot_var_pull_->addPlotable(hpull, "B");
    // We plot again without bars, so the points are not half covered by
    // bars as in case of "BP" draw option. We need to create and plot a
    // clone, because the ROOT object ownership is transfered to the RooPlot
    // by addPlotable().
    // If we just added the hpull twice, we would get a segfault.
    RooHist* hpull_clone = dynamic_cast<RooHist*>(hpull->Clone());
    // We plot again without bars, so the points are not half covered by bars
    plot_var_pull_->addPlotable(hpull_clone, "P");

    plot_var_pull_->GetXaxis()->SetTickLength(0.03 * pad_var->GetAbsHNDC() /
                                              pad_pull->GetAbsHNDC());
    plot_var_pull_->GetXaxis()->SetTitle(TString(var.GetTitle()));
    plot_var_pull_->GetXaxis()->SetTitleOffset(4.0);
    plot_var_pull_->GetXaxis()->SetLabelOffset(0.01 * pad_var->GetAbsHNDC() /
                                               pad_pull->GetAbsHNDC());
    plot_var_pull_->GetYaxis()->SetTitle("");
    // plot_var_pull_->GetYaxis()->SetRangeUser(0.8, 1.19);
    double max = std::max(- pull_histo->GetMinimum(), pull_histo->GetMaximum());
    plot_var_pull_->GetYaxis()->SetRangeUser(-max, +max);
    // plot_var_pull_->GetYaxis()->SetRangeUser(-5, 5);
    plot_var_pull_->GetYaxis()->SetNdivisions(505);
    plot_var_pull_->Draw();

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

void Fitter::PlotVars2D(RooRealVar& var1, RooRealVar& var2, RooDataHist& data1,
                        RooDataHist& data2) {
    if (canvas_var_) delete canvas_var_;
    canvas_var_ = new TCanvas(
        TString(var1.GetName()) + "_" + TString(var2.GetName()) + "_" + TString(data1.GetName()),
        TString(var1.GetTitle()) + "_" + TString(var2.GetName()) + "_" + TString(data1.GetTitle()),
        500, 500);

    TH2D* histo1 = static_cast<TH2D*>(data1.createHistogram("histo1", var1, RooFit::YVar(var2)));

    TH2D* histo2 = static_cast<TH2D*>(data2.createHistogram("histo2", var1, RooFit::YVar(var2)));

    // For some reason histo1->Clone("newname") doesn't return anything!
    TH2D* pull_histo =
        static_cast<TH2D*>(data2.createHistogram("pull_histo", var1, RooFit::YVar(var2)));

    const double max = std::max(histo1->GetMaximum(), histo2->GetMaximum());

    canvas_var_->SetRightMargin(0.14);

    histo1->SetMinimum(0);
    histo1->SetMaximum(max);
    histo1->SetTitle("");
    histo1->Draw("colz");
    histo1->GetZaxis()->SetTitle("");

    canvas_var_->Write();
    canvas_var_->SaveAs(format);

    histo2->SetMinimum(0);
    histo2->SetMaximum(max);
    histo2->SetTitle("");
    histo2->Draw("colz");
    histo2->GetZaxis()->SetTitle("");

    canvas_var_->SetName(TString(var1.GetName()) + "_" + TString(var2.GetName()) + "_" +
                         TString(data2.GetName()));
    canvas_var_->SetTitle(TString(var1.GetName()) + "_" + TString(var2.GetName()) + "_" +
                          TString(data2.GetTitle()));

    canvas_var_->Write();
    canvas_var_->SaveAs(format);

    double p1;
    double p2;
    for (int i = 1; i <= pull_histo->GetNbinsX(); i++) {
        for (int j = 1; j <= pull_histo->GetNbinsY(); j++) {
            p1 = histo1->GetBinContent(i, j);
            p2 = histo2->GetBinContent(i, j);
            // pull_histo->SetBinContent(i, j, p1/p2);
            pull_histo->SetBinContent(i, j, p1 - p2);
            // pull_histo->SetBinContent(i, j, (p1 - p2) / std::sqrt((p1 + p2) / 2));
        }
    }

// Check if we are running a newer version of ROOT which has new palettes
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 4, 0)
    gStyle->SetPalette(kLightTemperature);
#else
    gStyle->SetPalette(87);
#endif
    pull_histo->SetTitle("");
    // pull_histo->SetMinimum(0.8);
    // pull_histo->SetMaximum(1.2);
    const double pull_max = std::max(- pull_histo->GetMinimum(), pull_histo->GetMaximum());
    pull_histo->SetMinimum(-pull_max);
    pull_histo->SetMaximum(+pull_max);
    // pull_histo->SetMinimum(-5);
    // pull_histo->SetMaximum(5);
    pull_histo->GetZaxis()->SetTitle("");
    pull_histo->Draw("colz");
    pull_histo->Write();

    canvas_var_->SetName(TString(var1.GetName()) + "_" + TString(var2.GetName()) + "_" + TString(data1.GetName()) + "_" + TString(data2.GetName()));
    canvas_var_->SetTitle(TString(var1.GetTitle()) + "_" + TString(var2.GetTitle()) + " " + TString(data1.GetName()) + " " + TString(data2.GetName()));

    canvas_var_->Update();
    canvas_var_->Write();
    canvas_var_->SaveAs(format);
    colors::setColors();

    delete histo1;
    delete histo2;
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

    // This line makes sure the 0 is not drawn as it would overlap with the lower
    // pad
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

        // Create a new frame to draw the pull distribution and add the distribution
        // to the frame
        RooPlot* plot_eff_pull_ = var.frame(RooFit::Title("Pull Distribution"));
        plot_eff_pull_->SetTitle("");
        RooHist* hpull = plot_eff_->pullHist();
        hpull->SetFillColor(kGray);
        // The only working way to get rid of error bars; HIST draw option doesn't
        // work with RooPlot
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

    TH2D* gsim_histo =
        static_cast<TH2D*>(gsim_dataset_->createHistogram("gsim_histo", var1, RooFit::YVar(var2)));
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
    TH2* eff_pull_histo = dynamic_cast<TH2*>(gsim_histo->Clone("eff_pull_histo"));

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

    // We are not using eff_pull_histo->Add(pdf_eff_histo, -1) because we want
    // to skip bins which are empty in gsim_histo.
    for (int i = 1; i <= eff_pull_histo->GetNbinsX(); i++) {
        for (int j = 1; j <= eff_pull_histo->GetNbinsY(); j++) {
            if (eff_pull_histo->GetBinContent(i, j) != 0) {
                const double eff = eff_histo->GetBinContent(i, j);
                const double pdf = pdf_eff_histo->GetBinContent(i, j);
                eff_pull_histo->SetBinContent(i, j, eff - pdf);
            }
        }
    }

// Check if we are running a newer version of ROOT which has new palettes
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 4, 0)
    gStyle->SetPalette(kLightTemperature);
#else
    gStyle->SetPalette(87);
#endif
    eff_pull_histo->SetTitle("");
    eff_pull_histo->GetZaxis()->SetTitle();
    const double max = std::max(-eff_pull_histo->GetMinimum(), eff_pull_histo->GetMaximum());
    eff_pull_histo->SetMinimum(-max);
    eff_pull_histo->SetMaximum(max);
    eff_pull_histo->SetOption("colz 1");
    eff_pull_histo->Draw();
    eff_pull_histo->Write();

    canvas_eff_->SetName(TString(var1.GetName()) + "_" + TString(var2.GetName()) +
                         "_eff_residual_canvas");
    canvas_eff_->SetTitle(TString(var1.GetTitle()) + "_" + TString(var2.GetName()) +
                          " eff residual canvas");

    canvas_eff_->Update();
    canvas_eff_->Write();
    canvas_eff_->SaveAs(format);
    colors::setColors();

    delete eff_histo;
    delete gsim_histo;
    delete evtgen_histo;
    delete stat_box;
    delete eff_pull_histo;
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

void Fitter::ProcessBinnedEfficiency() {
    TEfficiency* binned_efficiency;

    TH3F* evtgen_histo =
        new TH3F("evtgen_histo", "evtgen_histo", 10, 0, kPi, 10, 0.5, 2.95, 10, -kPi, kPi);
    for (int i = 0; i < evtgen_dataset_->numEntries(); i++) {
        const RooArgSet* row = evtgen_dataset_->get(i);
        // printf("thetab = %f\n", row->getRealValue("thetab"));
        double thetat = row->getRealValue("thetat");
        double thetab = row->getRealValue("thetab");
        double phit = row->getRealValue("phit");
        evtgen_histo->Fill(thetat, thetab, phit);
    }

    TH3F* gsim_histo =
        new TH3F("gsim_histo", "gsim_histo", 10, 0, kPi, 10, 0.5, 2.95, 10, -kPi, kPi);
    for (int i = 0; i < gsim_dataset_->numEntries(); i++) {
        const RooArgSet* row = gsim_dataset_->get(i);
        // printf("thetab = %f\n", row->getRealValue("thetab"));
        double thetat = row->getRealValue("thetat");
        double thetab = row->getRealValue("thetab");
        double phit = row->getRealValue("phit");
        gsim_histo->Fill(thetat, thetab, phit);
    }

    // Check how many bins are empty
    printf("***** Bins where evtgen < gsim\n");
    for (int i = 0; i < 12 * 12 * 12; i++) {
        if (evtgen_histo->GetBinContent(i) < gsim_histo->GetBinContent(i)) {
            printf("Before: bin = %i, evtgen = %f, gsim = %f\n", i, evtgen_histo->GetBinContent(i),
                   gsim_histo->GetBinContent(i));
            gsim_histo->SetBinContent(i, evtgen_histo->GetBinContent(i));
            printf("After:  bin = %i, evtgen = %f, gsim = %f\n", i, evtgen_histo->GetBinContent(i),
                   gsim_histo->GetBinContent(i));
        }
    }
    printf("***** That's all\n");

    if (TEfficiency::CheckConsistency(*gsim_histo, *evtgen_histo)) {
        binned_efficiency = new TEfficiency(*gsim_histo, *evtgen_histo);
    }

    int global_bin = binned_efficiency->GetGlobalBin(kPi / 2, kPi / 2, kPi / 2);
    double particular_eff = binned_efficiency->GetEfficiency(global_bin);
    printf("bin = %i, Eff = %f\n", global_bin, particular_eff);
    for (int i = 0; i < 1000; i++) {
        printf("bin = %i, Eff = %f\n", i, binned_efficiency->GetEfficiency(i));
    }

    double thetat_bin_width = kPi / 10;
    double thetab_bin_width = (2.95 - 0.5) / 10;
    double phit_bin_width = 2 * kPi / 10;
    double thetat_start = 0 + thetat_bin_width / 2;
    double thetab_start = 0.5 + thetab_bin_width / 2;
    double phit_start = -kPi + phit_bin_width / 2;

    for (int t = 1; t <= 10; t++) {
        for (int b = 1; b <= 10; b++) {
            for (int p = 1; p <= 10; p++) {
                double tht = thetat_start + thetat_bin_width * (t - 1);
                double thb = thetab_start + thetab_bin_width * (b - 1);
                double pht = phit_start + phit_bin_width * (p - 1);
                int bin = evtgen_histo->GetBin(t, b, p);
                double evt = evtgen_histo->GetBinContent(bin);
                double gsim = gsim_histo->GetBinContent(bin);
                double eff = binned_efficiency->GetEfficiency(bin);
                printf(
                    "t = %i, b = %i, p = %i, tht = %f, thb = %f, pht = %f, bin = "
                    "%i, evt = %f, "
                    "gsim = %f, eff = %f\n",
                    t, b, p, tht, thb, pht, bin, evt, gsim, eff);
            }
        }
    }

    TCanvas* eff_canvas = new TCanvas("eff_evtgen_canvas", "eff evtgen canvas", 500, 500);
    TH2D* proj2D = dynamic_cast<TH2D*>(evtgen_histo->Project3D("xy"));
    proj2D->Draw();
    eff_canvas->SaveAs("test.png");
}

void Fitter::ProcessKDEEfficiency() {
    EnlargeVarRanges(0.1);
    MirrorDataAtEdges(evtgen_dataset_);
    MirrorDataAtEdges(gsim_dataset_);

    OneDimPhaseSpace phasespace_thetat("phasespace_thetat", thetat_.getMin(), thetat_.getMax());
    OneDimPhaseSpace phasespace_thetab("phasespace_thetab", thetab_.getMin(), thetab_.getMax());
    OneDimPhaseSpace phasespace_phit("phasespace_phit", phit_.getMin(), phit_.getMax());
    CombinedPhaseSpace phasespace("phasespace", &phasespace_thetat, &phasespace_thetab,
                                  &phasespace_phit);
    UniformDensity uniform_density("Uniform Density", &phasespace);
    FormulaDensity formula_density("Formula Density", &phasespace,
                                   "(x - 1.57)^2 * (y - 1.57)^2 * (z)^2");

    TH3F* eff_histo = GetBinned3DEfficiency();
    TTree* eff_tree = Histogram2TTree(eff_histo);

    BinnedKernelDensity bin_kde("BinKernelPDF", &phasespace, eff_tree, "thetat", "thetab", "phit",
                                "weight", 50, 50, 50, 0.20, 0.20, 0.40, 0);

    AdaptiveKernelDensity kde("KernelPDF", &phasespace, eff_tree, "thetat", "thetab", "phit",
                              "weight", 50, 50, 50, 0.10, 0.10, 0.20, &bin_kde);

    kde.writeToTextFile("efficiency.root");

    thetat_.setMin(0);
    thetat_.setMax(kPi);
    thetab_.setMin(0.5);
    thetab_.setMax(2.95);
    phit_.setMin(-kPi);
    phit_.setMax(kPi);
    vars_[0]->setMin(0);
    vars_[0]->setMax(kPi);
    vars_[1]->setMin(0.5);
    vars_[1]->setMax(2.95);
    vars_[2]->setMin(-kPi);
    vars_[2]->setMax(+kPi);

    eff_histo = GetBinned3DEfficiency();

    OneDimPhaseSpace phasespace_thetat_new("phasespace_thetat", thetat_.getMin(), thetat_.getMax());
    OneDimPhaseSpace phasespace_thetab_new("phasespace_thetab", thetab_.getMin(), thetab_.getMax());
    OneDimPhaseSpace phasespace_phit_new("phasespace_phit", phit_.getMin(), phit_.getMax());
    CombinedPhaseSpace phasespace_new("phasespace", &phasespace_thetat_new, &phasespace_thetab_new,
                                  &phasespace_phit_new);

    // printf("eff_tree->GetEntriesFast() = %i\n", eff_tree->GetEntriesFast());
    // KernelDensity kde("KernelPDF", &phasespace, eff_tree->GetEntriesFast(), 0.2,
    // 0.2, 0.4);
    // kde.readTuple(eff_tree, "thetat", "thetab", "phit");
    // kde.generateApproximation(eff_tree->GetEntriesFast());

    TH3F* binned_pdf = DensityToHisto(kde, phasespace_new);
    binned_pdf->Print();

    TH3F* simulated_histo = SimulateEfficiency(kde, evtgen_dataset_);
    TH3F* simulated_histo_dummy = SimulateEfficiencyDummy(eff_histo, evtgen_dataset_);
    TH3F* gsim_histo = Create3DHisto(gsim_dataset_);
    TH3F* evtgen_histo = Create3DHisto(evtgen_dataset_);
    TCanvas* eff_canvas = new TCanvas("eff_evtgen_canvas", "eff evtgen canvas", 500, 500);
    TH1D* simulated_projection;
    TH1D* gsim_projection;
    TH1D* eff_projection;
    TString name;

    printf("*** Numeric check ***\n");
    for (int i = 5000; i < 5010; i++) {
        printf("i = %i, eff = %f, gsim = %f, evtgen = %f, sim_dummy = %f, sim = %f\n", i,
               eff_histo->GetBinContent(i), gsim_histo->GetBinContent(i),
               evtgen_histo->GetBinContent(i), simulated_histo_dummy->GetBinContent(i),
               simulated_histo->GetBinContent(i));
    }
    for (int i = 50000; i < 50010; i++) {
        printf("i = %i, eff = %f, gsim = %f, evtgen = %f, sim_dummy = %f, sim = %f\n", i,
               eff_histo->GetBinContent(i), gsim_histo->GetBinContent(i),
               evtgen_histo->GetBinContent(i), simulated_histo_dummy->GetBinContent(i),
               simulated_histo->GetBinContent(i));
    }
    // printf("*** Non-matching bins ***\n");
    // for (int i = 0; i < 52 * 52 * 52; i++) {
    //     if (((int)gsim_histo->GetBinContent(i)) !=
    //     ((int)simulated_histo_dummy->GetBinContent(i))) {
    //         printf("i = %i, eff = %f, gsim = %f, evtgen = %f, sim_dummy = %f, sim
    //         = %f\n", i,
    //                eff_histo->GetBinContent(i), gsim_histo->GetBinContent(i),
    //                evtgen_histo->GetBinContent(i),
    //                simulated_histo_dummy->GetBinContent(i),
    //                simulated_histo->GetBinContent(i));
    //     }
    // }

    // printf("*** < 1 gsim bins ***\n");
    // for (int i = 0; i < 52 * 52 * 52; i++) {
    //     if (gsim_histo->GetBinContent(i) < 1.0) {
    //         printf("i = %i, eff = %f, gsim = %f, evtgen = %f, sim_dummy = %f, sim
    //         = %f\n", i,
    //                eff_histo->GetBinContent(i), gsim_histo->GetBinContent(i),
    //                evtgen_histo->GetBinContent(i),
    //                simulated_histo_dummy->GetBinContent(i),
    //                simulated_histo->GetBinContent(i));
    //     }
    // }

    printf("total: gsim = %f, sim dummy = %f, sim = %f\n", gsim_histo->GetSumOfWeights(),
           simulated_histo_dummy->GetSumOfWeights(), simulated_histo->GetSumOfWeights());

    // TString var[] = {'x', 'y', 'z'};
    // for (int var_num = 0; var_num < 3; var_num++) {
    //     name = "proj_1D_sim_";
    //     name += vars_[var_num]->GetName();
    //     eff_canvas->SetName(name);
    //     simulated_projection = dynamic_cast<TH1D*>(simulated_histo->Project3D(var[var_num]));
    //     simulated_projection->Draw();
    //     eff_canvas->SaveAs(format);
    // }

    // for (int var_num = 0; var_num < 3; var_num++) {
    //     name = "proj_1D_sim_dummy_";
    //     name += vars_[var_num]->GetName();
    //     eff_canvas->SetName(name);
    //     simulated_projection = dynamic_cast<TH1D*>(simulated_histo_dummy->Project3D(var[var_num]));
    //     simulated_projection->Draw();
    //     eff_canvas->SaveAs(format);
    // }

    // for (int var_num = 0; var_num < 3; var_num++) {
    //     name = "proj_1D_gsim_";
    //     name += vars_[var_num]->GetName();
    //     eff_canvas->SetName(name);
    //     gsim_projection = dynamic_cast<TH1D*>(gsim_histo->Project3D(var[var_num]));
    //     gsim_projection->Draw();
    //     eff_canvas->SaveAs(format);
    // }

    // for (int var_num = 0; var_num < 3; var_num++) {
    //     name = "proj_1D_gsim_and_sim_";
    //     name += vars_[var_num]->GetName();
    //     eff_canvas->SetName(name);
    //     gsim_projection = dynamic_cast<TH1D*>(gsim_histo->Project3D(var[var_num]));
    //     gsim_projection->Draw();
    //     simulated_projection = dynamic_cast<TH1D*>(simulated_histo->Project3D(var[var_num]));
    //     simulated_projection->SetLineColor(kRed);
    //     simulated_projection->Draw("same");
    //     eff_canvas->SaveAs(format);
    // }

    // for (int var_num = 0; var_num < 3; var_num++) {
    //     name = "dummy_vs_gsim_1D_";
    //     name += vars_[var_num]->GetName();
    //     eff_canvas->SetName(name);
    //     simulated_projection = dynamic_cast<TH1D*>(simulated_histo_dummy->Project3D(var[var_num]));
    //     gsim_projection = dynamic_cast<TH1D*>(gsim_histo->Project3D(var[var_num]));
    //     eff_projection = dynamic_cast<TH1D*>(simulated_projection->Clone("eff projection"));
    //     eff_projection->Divide(simulated_projection, gsim_projection);
    //     eff_projection->Draw();
    //     eff_canvas->SaveAs(format);
    // }

    // for (int var_num = 0; var_num < 3; var_num++) {
    //     name = "dummy_minus_gsim_1D_";
    //     name += vars_[var_num]->GetName();
    //     eff_canvas->SetName(name);
    //     simulated_projection = dynamic_cast<TH1D*>(simulated_histo_dummy->Project3D(var[var_num]));
    //     gsim_projection = dynamic_cast<TH1D*>(gsim_histo->Project3D(var[var_num]));
    //     eff_projection = dynamic_cast<TH1D*>(simulated_projection->Clone("eff projection"));
    //     eff_projection->Add(gsim_projection, -1);
    //     eff_projection->Draw();
    //     eff_canvas->SaveAs(format);
    // }

    // for (int var_num = 0; var_num < 3; var_num++) {
    //     name = "sim_vs_gsim_1D_";
    //     name += vars_[var_num]->GetName();
    //     eff_canvas->SetName(name);
    //     simulated_projection = dynamic_cast<TH1D*>(simulated_histo->Project3D(var[var_num]));
    //     gsim_projection = dynamic_cast<TH1D*>(gsim_histo->Project3D(var[var_num]));
    //     eff_projection = dynamic_cast<TH1D*>(simulated_projection->Clone("eff projection"));
    //     eff_projection->Divide(simulated_projection, gsim_projection);
    //     eff_projection->Draw();
    //     eff_canvas->SaveAs(format);
    // }

    // for (int var_num = 0; var_num < 3; var_num++) {
    //     name = "sim_minus_gsim_1D_";
    //     name += vars_[var_num]->GetName();
    //     eff_canvas->SetName(name);
    //     simulated_projection = dynamic_cast<TH1D*>(simulated_histo->Project3D(var[var_num]));
    //     gsim_projection = dynamic_cast<TH1D*>(gsim_histo->Project3D(var[var_num]));
    //     eff_projection = dynamic_cast<TH1D*>(simulated_projection->Clone("eff projection"));
    //     eff_projection->Add(gsim_projection, -1);
    //     eff_projection->Draw();
    //     eff_canvas->SaveAs(format);
    // }

    // for (int var_num = 0; var_num < 3; var_num++) {
    //     name = "eff_1D_";
    //     name += vars_[var_num]->GetName();
    //     eff_canvas->SetName(name);
    //     eff_projection = dynamic_cast<TH1D*>(eff_histo->Project3D(var[var_num]));
    //     eff_projection->Draw();
    //     eff_canvas->SaveAs(format);
    // }

    // for (int var_num = 0; var_num < 3; var_num++) {
    //     name = "eff_pdf_1D_";
    //     name += vars_[var_num]->GetName();
    //     eff_canvas->SetName(name);
    //     eff_projection = dynamic_cast<TH1D*>(binned_pdf->Project3D(var[var_num]));
    //     eff_projection->Draw();
    //     eff_canvas->SaveAs(format);
    // }

    // TH2D* pdf_projection_2D;
    // for (int var_num = 0; var_num < 3; var_num++) {
    //     name = "proj_2D";
    //     for (int i = 0; i < 3; i++) {
    //         if (i != var_num) {
    //             name += "_";
    //             name += vars_[i].GetName();
    //         }
    //     }
    //     name += ".png";
    //     pdf_projection_2D = ProjectDensityTo2D(kde, phasespace, var_num);
    //     pdf_projection_2D->Draw("colz");
    //     eff_canvas->SaveAs(name);
    // }

    RooDataHist roo_simulated_histo("roo_simulated_histo", "roo_simulated_histo",
                                    RooArgList(thetat_, thetab_, phit_), simulated_histo);
    RooDataHist roo_gsim_histo("roo_gsim_histo", "roo_gsim_histo",
                               RooArgList(thetat_, thetab_, phit_), gsim_histo);

    for (auto&& var : vars_) {
        PlotVar(*var, roo_simulated_histo, roo_gsim_histo);
    }

    for (int i = 0; i < 3; i++) {
        for (int j = i + 1; j < 3; j++) {
            PlotVars2D(*vars_[i], *vars_[j], roo_simulated_histo, roo_gsim_histo);
        }
    }

    // Scale the histos so that they have efficiency averages on y (or z) axes
    eff_histo->Scale(1./50.);
    RooDataHist roo_eff_histo_2D("roo_eff_histo", "roo_eff_histo", RooArgList(thetat_, thetab_, phit_),
                              eff_histo);
    TH3F* scaled_binned_pdf = dynamic_cast<TH3F*>(binned_pdf->Clone("scaled_binned_pdf"));

    double scale = eff_histo->GetSumOfWeights()/binned_pdf->GetSumOfWeights();
    scaled_binned_pdf->Scale(scale);

    RooDataHist roo_eff_pdf_histo_2D("roo_eff_pdf_histo", "roo_eff_pdf_histo",
                                  RooArgList(thetat_, thetab_, phit_), scaled_binned_pdf);

    eff_histo->Scale(1./50.);
    RooDataHist roo_eff_histo_1D("roo_eff_histo", "roo_eff_histo", RooArgList(thetat_, thetab_, phit_),
                              eff_histo);
    scaled_binned_pdf->Scale(1./50.);
    RooDataHist roo_eff_pdf_histo_1D("roo_eff_pdf_histo", "roo_eff_pdf_histo",
                                  RooArgList(thetat_, thetab_, phit_), scaled_binned_pdf);

    for (auto&& var : vars_) {
        PlotVar(*var, roo_eff_histo_1D, roo_eff_pdf_histo_1D);
    }

    for (int i = 0; i < 3; i++) {
        for (int j = i + 1; j < 3; j++) {
            PlotVars2D(*vars_[i], *vars_[j], roo_eff_histo_2D, roo_eff_pdf_histo_2D);
        }
    }
}

TH3F* Fitter::GetBinned3DEfficiency() {
    // printf("Trying to get binned efficiency.\n");
    const int num_bins = 50;
    TH3F* evtgen_histo = new TH3F("evtgen_histo", "evtgen_histo", num_bins, thetat_.getMin(),
                                  thetat_.getMax(), num_bins, thetab_.getMin(), thetab_.getMax(),
                                  num_bins, phit_.getMin(), phit_.getMax());
    for (int i = 0; i < evtgen_dataset_->numEntries(); i++) {
        const RooArgSet* row = evtgen_dataset_->get(i);
        // printf("thetab = %f\n", row->getRealValue("thetab"));
        double thetat = row->getRealValue("thetat");
        double thetab = row->getRealValue("thetab");
        double phit = row->getRealValue("phit");
        evtgen_histo->Fill(thetat, thetab, phit);
    }

    TH3F* gsim_histo =
        new TH3F("gsim_histo", "gsim_histo", num_bins, thetat_.getMin(), thetat_.getMax(), num_bins,
                 thetab_.getMin(), thetab_.getMax(), num_bins, phit_.getMin(), phit_.getMax());
    for (int i = 0; i < gsim_dataset_->numEntries(); i++) {
        const RooArgSet* row = gsim_dataset_->get(i);
        // printf("thetab = %f\n", row->getRealValue("thetab"));
        double thetat = row->getRealValue("thetat");
        double thetab = row->getRealValue("thetab");
        double phit = row->getRealValue("phit");
        gsim_histo->Fill(thetat, thetab, phit);
    }

    // TH1D* evtgen_histo_1D = evtgen_histo->ProjectionZ();
    // TH1D* gsim_histo_1D = gsim_histo->ProjectionZ();
    // TH1D* eff_histo_1D =
    // dynamic_cast<TH1D*>(gsim_histo_1D->Clone("eff_histo_1D"));
    // eff_histo_1D->Divide(gsim_histo_1D, evtgen_histo_1D, 1.0, 1.0, "B");

    // TCanvas* canvas = new TCanvas("canvas", "canvas", 500, 500);
    // eff_histo_1D->Draw();
    // canvas->SaveAs("eff_1D.png");

    // printf("Trying to clone and divide binned histos.\n");
    TH3F* eff_histo = dynamic_cast<TH3F*>(gsim_histo->Clone("eff_histo"));
    eff_histo->Divide(gsim_histo, evtgen_histo);  //, 1.0, 1.0, "B");

    // TH1D* histo_projection;
    // histo_projection = evtgen_histo->ProjectionZ();
    // histo_projection->Draw();
    // canvas->SaveAs("evtgen.png");

    // histo_projection = gsim_histo->ProjectionZ();
    // histo_projection->Draw();
    // canvas->SaveAs("gsim.png");

    // histo_projection = eff_histo->ProjectionZ();
    // histo_projection->Draw();
    // canvas->SaveAs("eff.png");

    // evtgen_histo->Print();
    // gsim_histo->Print();
    // eff_histo->Print();
    // printf("Trying to delete binned histos.\n");
    delete gsim_histo;
    delete evtgen_histo;

    return eff_histo;
}

TTree* Fitter::Histogram2TTree(TH3F* histo) {
    // Creates a TTree and fills it with the coordinates of all
    // filled bins. The tree will have one branch for each dimension,
    // and one for the bin content.
    TString name(histo->GetName());
    name += "_tree";
    TString title(histo->GetTitle());
    title += " tree";
    TTree* tree = new TTree(name, title);
    float thetat;
    float thetab;
    float phit;
    float weight;
    tree->Branch("thetat", &thetat, "thetat/F");
    tree->Branch("thetab", &thetab, "thetab/F");
    tree->Branch("phit", &phit, "phit/F");
    tree->Branch("weight", &weight, "weight/F");

    for (int x = 1; x <= histo->GetXaxis()->GetNbins(); x++) {
        for (int y = 1; y <= histo->GetYaxis()->GetNbins(); y++) {
            for (int z = 1; z <= histo->GetZaxis()->GetNbins(); z++) {
                thetat = histo->GetXaxis()->GetBinCenter(x);
                thetab = histo->GetYaxis()->GetBinCenter(y);
                phit = histo->GetZaxis()->GetBinCenter(z);
                weight = histo->GetBinContent(histo->GetBin(x, y, z));
                tree->Fill();
            }
        }
    }
    tree->Print("v");

    return tree;
}

TH1D* Fitter::ProjectDensityTo1D(BinnedKernelDensity pdf, CombinedPhaseSpace phasespace,
                                 const int var_num) {
    TH3F* binned_pdf =
        new TH3F("binned_pdf", "Binned PDF", 100, phasespace.lowerLimit(0),
                 phasespace.upperLimit(0), 100, phasespace.lowerLimit(1), phasespace.upperLimit(1),
                 100, phasespace.lowerLimit(2), phasespace.upperLimit(2));

    double thetat;
    double thetab;
    double phit;
    for (int x = 1; x <= binned_pdf->GetXaxis()->GetNbins(); x++) {
        for (int y = 1; y <= binned_pdf->GetYaxis()->GetNbins(); y++) {
            for (int z = 1; z <= binned_pdf->GetZaxis()->GetNbins(); z++) {
                thetat = binned_pdf->GetXaxis()->GetBinCenter(x);
                thetab = binned_pdf->GetYaxis()->GetBinCenter(y);
                phit = binned_pdf->GetZaxis()->GetBinCenter(z);
                std::vector<double> coords;
                coords.push_back(thetat);
                coords.push_back(thetab);
                coords.push_back(phit);
                binned_pdf->SetBinContent(binned_pdf->GetBin(x, y, z), pdf.density(coords));
            }
        }
    }

    TH1D* pdf_projection;
    switch (var_num) {
        case 0:
            pdf_projection = binned_pdf->ProjectionX();
            break;
        case 1:
            pdf_projection = binned_pdf->ProjectionY();
            break;
        case 2:
            pdf_projection = binned_pdf->ProjectionZ();
            break;
    }
    return pdf_projection;
}

TH2D* Fitter::ProjectDensityTo2D(BinnedKernelDensity pdf, CombinedPhaseSpace phasespace,
                                 const int var_num) {
    TH3F* binned_pdf =
        new TH3F("binned_pdf", "Binned PDF", 100, phasespace.lowerLimit(0),
                 phasespace.upperLimit(0), 100, phasespace.lowerLimit(1), phasespace.upperLimit(1),
                 100, phasespace.lowerLimit(2), phasespace.upperLimit(2));

    double thetat;
    double thetab;
    double phit;
    for (int x = 1; x <= binned_pdf->GetXaxis()->GetNbins(); x++) {
        for (int y = 1; y <= binned_pdf->GetYaxis()->GetNbins(); y++) {
            for (int z = 1; z <= binned_pdf->GetZaxis()->GetNbins(); z++) {
                thetat = binned_pdf->GetXaxis()->GetBinCenter(x);
                thetab = binned_pdf->GetYaxis()->GetBinCenter(y);
                phit = binned_pdf->GetZaxis()->GetBinCenter(z);
                std::vector<double> coords;
                coords.push_back(thetat);
                coords.push_back(thetab);
                coords.push_back(phit);
                binned_pdf->SetBinContent(binned_pdf->GetBin(x, y, z), pdf.density(coords));
            }
        }
    }

    TH2D* pdf_projection;
    switch (var_num) {
        case 0:
            pdf_projection = dynamic_cast<TH2D*>(binned_pdf->Project3D("yz"));
            break;
        case 1:
            pdf_projection = dynamic_cast<TH2D*>(binned_pdf->Project3D("xz"));
            break;
        case 2:
            pdf_projection = dynamic_cast<TH2D*>(binned_pdf->Project3D("xy"));
            break;
    }
    return pdf_projection;
}

void Fitter::Process1DKDEEfficiency() {
    OneDimPhaseSpace phasespace("phasespace", phit_.getMin(), phit_.getMax());

    TH3F* eff_histo = GetBinned3DEfficiency();
    TTree* eff_tree_3D = Histogram2TTree(eff_histo);
    eff_tree_3D->SetBranchStatus("*", 0);
    eff_tree_3D->SetBranchStatus("phit", 1);
    eff_tree_3D->SetBranchStatus("weight", 1);
    TTree* eff_tree = eff_tree_3D->CloneTree();

    BinnedKernelDensity kde("KernelPDF", &phasespace, eff_tree, "phit", "weight", 100, 0.05, 0, 0);

    TCanvas* eff_canvas = new TCanvas("eff_evtgen_canvas", "eff evtgen canvas", 500, 500);
    TH1F* projection = new TH1F("projection", "projection", 100, -kPi, kPi);
    TH1D* histo_projection;
    kde.project(projection);
    // projection->Draw();
    histo_projection = eff_histo->ProjectionZ();
    histo_projection->Draw();
    // eff_tree->Draw("weight");
    eff_canvas->SaveAs("phit.png");
}

TH3F* Fitter::SimulateEfficiency(AdaptiveKernelDensity pdf, const RooDataSet* dataset) {
    TH3F* evtgen_histo = Create3DHisto(evtgen_dataset_);
    const int num_bins = 50;
    TH3F* simulated_histo = new TH3F("simulated_histo", "simulated_histo", num_bins,
                                     thetat_.getMin(), thetat_.getMax(), num_bins, thetab_.getMin(),
                                     thetab_.getMax(), num_bins, phit_.getMin(), phit_.getMax());

    const double total_efficiency = 0.105;
    double thetat;
    double thetab;
    double phit;
    int bin;
    double evtgen_value;
    for (int x = 1; x <= simulated_histo->GetXaxis()->GetNbins(); x++) {
        for (int y = 1; y <= simulated_histo->GetYaxis()->GetNbins(); y++) {
            for (int z = 1; z <= simulated_histo->GetZaxis()->GetNbins(); z++) {
                thetat = simulated_histo->GetXaxis()->GetBinCenter(x);
                thetab = simulated_histo->GetYaxis()->GetBinCenter(y);
                phit = simulated_histo->GetZaxis()->GetBinCenter(z);
                std::vector<double> coords;
                coords.push_back(thetat);
                coords.push_back(thetab);
                coords.push_back(phit);
                bin = simulated_histo->GetBin(x, y, z);
                evtgen_value = evtgen_histo->GetBinContent(bin);
                simulated_histo->SetBinContent(
                    bin, pdf.density(coords) * evtgen_value * total_efficiency);
            }
        }
    }

    return simulated_histo;
}

TH3F* Fitter::SimulateEfficiencyDummy(TH3F* eff_histo, const RooDataSet* dataset) {
    TH3F* evtgen_histo = Create3DHisto(evtgen_dataset_);
    const int num_bins = 50;
    TH3F* simulated_histo = new TH3F("simulated_histo", "simulated_histo", num_bins,
                                     thetat_.getMin(), thetat_.getMax(), num_bins, thetab_.getMin(),
                                     thetab_.getMax(), num_bins, phit_.getMin(), phit_.getMax());

    const double total_efficiency = 1;
    double thetat;
    double thetab;
    double phit;
    int bin;
    double eff_value;
    double evtgen_value;
    for (int x = 1; x <= simulated_histo->GetXaxis()->GetNbins(); x++) {
        for (int y = 1; y <= simulated_histo->GetYaxis()->GetNbins(); y++) {
            for (int z = 1; z <= simulated_histo->GetZaxis()->GetNbins(); z++) {
                thetat = simulated_histo->GetXaxis()->GetBinCenter(x);
                thetab = simulated_histo->GetYaxis()->GetBinCenter(y);
                phit = simulated_histo->GetZaxis()->GetBinCenter(z);
                std::vector<double> coords;
                coords.push_back(thetat);
                coords.push_back(thetab);
                coords.push_back(phit);
                bin = simulated_histo->GetBin(x, y, z);
                evtgen_value = evtgen_histo->GetBinContent(bin);
                eff_value = eff_histo->GetBinContent(bin);
                simulated_histo->SetBinContent(bin, eff_value * evtgen_value * total_efficiency);
            }
        }
    }

    return simulated_histo;
}

TH3F* Fitter::Create3DHisto(const RooDataSet* dataset) {
    const int num_bins = 50;
    TH3F* histo =
        new TH3F("histo", "histo", num_bins, thetat_.getMin(), thetat_.getMax(), num_bins,
                 thetab_.getMin(), thetab_.getMax(), num_bins, phit_.getMin(), phit_.getMax());
    for (int i = 0; i < dataset->numEntries(); i++) {
        const RooArgSet* row = dataset->get(i);
        // printf("thetab = %f\n", row->getRealValue("thetab"));
        double thetat = row->getRealValue("thetat");
        double thetab = row->getRealValue("thetab");
        double phit = row->getRealValue("phit");
        histo->Fill(thetat, thetab, phit);
    }

    return histo;
}

TH3F* Fitter::DensityToHisto(AdaptiveKernelDensity pdf, CombinedPhaseSpace phasespace) {
    TH3F* binned_pdf =
        new TH3F("binned_pdf", "Binned PDF", 50, phasespace.lowerLimit(0),
                 phasespace.upperLimit(0), 50, phasespace.lowerLimit(1), phasespace.upperLimit(1),
                 50, phasespace.lowerLimit(2), phasespace.upperLimit(2));

    double thetat;
    double thetab;
    double phit;
    for (int x = 1; x <= binned_pdf->GetXaxis()->GetNbins(); x++) {
        for (int y = 1; y <= binned_pdf->GetYaxis()->GetNbins(); y++) {
            for (int z = 1; z <= binned_pdf->GetZaxis()->GetNbins(); z++) {
                thetat = binned_pdf->GetXaxis()->GetBinCenter(x);
                thetab = binned_pdf->GetYaxis()->GetBinCenter(y);
                phit = binned_pdf->GetZaxis()->GetBinCenter(z);
                std::vector<double> coords;
                coords.push_back(thetat);
                coords.push_back(thetab);
                coords.push_back(phit);
                binned_pdf->SetBinContent(binned_pdf->GetBin(x, y, z), pdf.density(coords));
            }
        }
    }

    return binned_pdf;
}
