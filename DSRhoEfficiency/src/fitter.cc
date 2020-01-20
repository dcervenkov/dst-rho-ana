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
#include <cstring>

// Boost includes
#include <boost/algorithm/string/replace.hpp>

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
#include "RooTreeDataStore.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TEnv.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3F.h"
#include "TLegend.h"
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
#include "log.h"
#include "tools.h"

Fitter::Fitter(const char* evtgen_filepath, const char* gsim_filepath, const char* output_dir) {
    dec_type_.defineType("a", 1);
    dec_type_.defineType("ab", 2);
    dec_type_.defineType("b", 3);
    dec_type_.defineType("bb", 4);

    if (strstr(evtgen_filepath, ".root")) {
        TFile* evtgen_file = new TFile(evtgen_filepath);
        TTree* evtgen_tree = dynamic_cast<TTree*>(evtgen_file->Get("h2000"));
        evtgen_dataset_ = new RooDataSet(
            "evtgen_dataset", "evtgen dataset", evtgen_tree,
            RooArgSet(thetat_, thetab_, phit_, vrvtxz_, vtvtxz_, evmcflag_), "evmcflag==1");
        delete evtgen_tree;
        evtgen_file->Close();
    } else {
        evtgen_dataset_ =
            RooDataSet::read(evtgen_filepath, RooArgList(thetat_, thetab_, phit_, dt_, dec_type_), "Q");
    }
    Log::print(Log::info, "Read %i Evtgen events\n", evtgen_dataset_->numEntries());

    TString cuts = tools::GetCommonCutsString();
    cuts += "&&evmcflag==1";

    RooArgSet dataset_vars;
    dataset_vars.add(csbdtg_);
    dataset_vars.add(de_);
    dataset_vars.add(dt_);
    dataset_vars.add(evmcflag_);
    dataset_vars.add(thetat_);
    dataset_vars.add(thetab_);
    dataset_vars.add(phit_);

    dataset_vars.add(vrusable_);
    dataset_vars.add(vrvtxz_);
    dataset_vars.add(vrerr6_);
    dataset_vars.add(vrchi2_);
    dataset_vars.add(vrndf_);
    dataset_vars.add(vrntrk_);

    dataset_vars.add(vtusable_);
    dataset_vars.add(vtvtxz_);
    dataset_vars.add(vterr6_);
    dataset_vars.add(vtchi2_);
    dataset_vars.add(vtndf_);
    dataset_vars.add(vtntrk_);

    TFile* gsim_file = new TFile(gsim_filepath);
    TTree* gsim_tree = dynamic_cast<TTree*>(gsim_file->Get("h2000"));
    gsim_dataset_ = new RooDataSet("gsim_dataset", "gsim dataset", gsim_tree, dataset_vars, cuts);

    // Add calculated dt to the dataset
    // gsim_dataset_->addColumn(dt_formula_);
    delete gsim_tree;
    gsim_file->Close();
    Log::print(Log::info, "Read %i Gsim events\n", gsim_dataset_->numEntries());

    tools::SetPlotDir(output_dir);

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
        const RooArgSet* row = data->get(evt);
        double vals[3] = {row->getRealValue("thetat"), row->getRealValue("thetab"),
                          row->getRealValue("phit")};
        for (int var_num1 = 0; var_num1 < 3; var_num1++) {
            if (CloseToEdge(vals, var_num1)) {
                double* mirror_vals = GetMirrorVals(vals, var_num1);
                *vars_[0] = mirror_vals[0];
                *vars_[1] = mirror_vals[1];
                *vars_[2] = mirror_vals[2];
                delete[] mirror_vals;
                data->add(RooArgSet(*vars_[0], *vars_[1], *vars_[2]));
                num_added++;
                for (int var_num2 = var_num1 + 1; var_num2 < 3; var_num2++) {
                    if (CloseToEdge(vals, var_num2)) {
                        double* mirror_vals = GetMirrorVals(vals, var_num1, var_num2);
                        *vars_[0] = mirror_vals[0];
                        *vars_[1] = mirror_vals[1];
                        *vars_[2] = mirror_vals[2];
                        delete[] mirror_vals;
                        data->add(RooArgSet(*vars_[0], *vars_[1], *vars_[2]));
                        num_added++;
                        for (int var_num3 = var_num2 + 1; var_num3 < 3; var_num3++) {
                            if (CloseToEdge(vals, var_num3)) {
                                double* mirror_vals =
                                    GetMirrorVals(vals, var_num1, var_num2, var_num3);
                                *vars_[0] = mirror_vals[0];
                                *vars_[1] = mirror_vals[1];
                                *vars_[2] = mirror_vals[2];
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
    Log::print(Log::debug, "Added %i mirrored datapoints.\n", num_added);
}

double* Fitter::GetMirrorVals(double vals[], const int var_num) {
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

double* Fitter::GetMirrorVals(double vals[], const int var_num1, const int var_num2) {
    double* new_vals1 = GetMirrorVals(vals, var_num1);
    double* new_vals2 = GetMirrorVals(new_vals1, var_num2);
    delete[] new_vals1;
    return new_vals2;
}

double* Fitter::GetMirrorVals(double vals[], const int var_num1, const int var_num2,
                              const int var_num3) {
    double* new_vals1 = GetMirrorVals(vals, var_num1);
    double* new_vals2 = GetMirrorVals(new_vals1, var_num2);
    double* new_vals3 = GetMirrorVals(new_vals2, var_num3);
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

void Fitter::PlotVar(const RooRealVar& var) const {
    RooDataHist evtgen_datahist("evtgen_datahist", "evtgen_datahist", RooArgSet(var),
                                *evtgen_dataset_);
    RooDataHist gsim_datahist("gsim_datahist", "gsim_datahist", RooArgSet(var), *gsim_dataset_);
    tools::PlotVar(var, evtgen_datahist, gsim_datahist, false, false);
}

void Fitter::PlotVars2D(const RooRealVar& var1, const RooRealVar& var2) const {
    tools::PlotVars2D(var1, var2, *evtgen_dataset_, "", constants::format);
    tools::PlotVars2D(var1, var2, *gsim_dataset_, "", constants::format);
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

    // Set the current directory back to the one for plots (ugly ROOT stuff)
    if (output_file_) {
        output_file_->cd();
    }
    canvas_eff_->Write();
    canvas_eff_->SaveAs(constants::format);
}

void Fitter::PlotEfficiency2D(RooRealVar& var1, RooRealVar& var2) {
    if (canvas_eff_) delete canvas_eff_;
    canvas_eff_ = new TCanvas(
        TString(var1.GetName()) + "_" + TString(var2.GetName()) + "_eff_canvas",
        TString(var1.GetTitle()) + "_" + TString(var2.GetName()) + " eff canvas", 500, 500);

    TH3F* eff_histo_3D = GetBinned3DEfficiency();

    std::string projection_string;
    projection_string += var2.GetName();
    projection_string += var1.GetName();
    boost::replace_all(projection_string, "thetat", "x");
    boost::replace_all(projection_string, "thetab", "y");
    boost::replace_all(projection_string, "phit", "z");

    TH2D* eff_histo = dynamic_cast<TH2D*>(eff_histo_3D->Project3D(projection_string.c_str()));

    double normalization = 1;
    if (projection_string.find("x") == std::string::npos) {
        normalization /= eff_histo_3D->GetNbinsX();
    } else if (projection_string.find("y") == std::string::npos) {
        normalization /= eff_histo_3D->GetNbinsY();
    } else if (projection_string.find("z") == std::string::npos) {
        normalization /= eff_histo_3D->GetNbinsZ();
    }

    eff_histo->Scale(normalization);


    Log::print(Log::info, "Correlation %s : %s is %f\n", var1.GetName(), var2.GetName(),
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
    canvas_eff_->SaveAs(constants::format);

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
    pdf_eff_histo->SetMaximum(eff_histo->GetMaximum());
    pdf_eff_histo->Draw("colz");
    pdf_eff_histo->GetZaxis()->SetTitle("");

    canvas_eff_->SetName(TString(var1.GetName()) + "_" + TString(var2.GetName()) +
                         "_pdf_eff_canvas");
    canvas_eff_->SetTitle(TString(var1.GetTitle()) + "_" + TString(var2.GetName()) +
                          " PDF eff canvas");

    // Set the current directory back to the one for plots (ugly ROOT stuff)
    if (output_file_) {
        output_file_->cd();
    }
    canvas_eff_->Write();
    canvas_eff_->SaveAs(constants::format);

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
    // eff_pull_histo->SetMinimum(-0.1);
    // eff_pull_histo->SetMaximum(0.1);
    eff_pull_histo->SetOption("colz 1");
    eff_pull_histo->Draw();
    eff_pull_histo->Write();

    canvas_eff_->SetName(TString(var1.GetName()) + "_" + TString(var2.GetName()) +
                         "_eff_residual_canvas");
    canvas_eff_->SetTitle(TString(var1.GetTitle()) + "_" + TString(var2.GetName()) +
                          " eff residual canvas");

    canvas_eff_->Update();
    canvas_eff_->Write();
    canvas_eff_->SaveAs(constants::format);
    colors::setColors();

    delete eff_histo;
    delete eff_histo_3D;
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
            thetat_model_e_ = &thetat_model1_;
            thetab_model_e_ = &thetab_model1_;
            phit_model_e_ = &phit_model1_;
            break;
        case 2:
            thetat_model_e_ = &thetat_model2_;
            thetab_model_e_ = &thetab_model2_;
            phit_model_e_ = &phit_model2_;
            break;
        case 3:
            thetat_model_e_ = &thetat_model3_;
            thetab_model_e_ = &thetab_model3_;
            phit_model_e_ = &phit_model3_;
            break;
        case 4:
            thetat_model_e_ = &thetat_model4_;
            thetab_model_e_ = &thetab_model3_;
            phit_model_e_ = &phit_model2_;
            break;
    }
}

void Fitter::ProcessKDEEfficiency(const char* efficiency_file,
                                  const std::array<double, 6> bin_kde_pars,
                                  const std::array<double, 6> ada_kde_pars,
                                  const double mirror_margin) {
    if (mirror_margin) {
        EnlargeVarRanges(mirror_margin);
        MirrorDataAtEdges(evtgen_dataset_);
        MirrorDataAtEdges(gsim_dataset_);
    }

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
                                "weight", bin_kde_pars[0], bin_kde_pars[1], bin_kde_pars[2],
                                bin_kde_pars[3], bin_kde_pars[4], bin_kde_pars[5], 0);

    AdaptiveKernelDensity kde("KernelPDF", &phasespace, eff_tree, "thetat", "thetab", "phit",
                              "weight", ada_kde_pars[0], ada_kde_pars[1], ada_kde_pars[2],
                              ada_kde_pars[3], ada_kde_pars[4], ada_kde_pars[5], &bin_kde);
    kde.writeToTextFile(efficiency_file);

    if (mirror_margin) {
        thetat_.setMin(constants::cuts::thetat_low);
        thetat_.setMax(constants::cuts::thetat_high);
        thetab_.setMin(constants::cuts::thetab_low);
        thetab_.setMax(constants::cuts::thetab_high);
        phit_.setMin(constants::cuts::phit_low);
        phit_.setMax(constants::cuts::phit_high);
        vars_[0]->setMin(constants::cuts::thetat_low);
        vars_[0]->setMax(constants::cuts::thetat_high);
        vars_[1]->setMin(constants::cuts::thetab_low);
        vars_[1]->setMax(constants::cuts::thetab_high);
        vars_[2]->setMin(constants::cuts::phit_low);
        vars_[2]->setMax(constants::cuts::phit_high);
        delete eff_histo;
        eff_histo = GetBinned3DEfficiency();
    }

    // TFile new_efficiency_file("efficiency.root", "read");
    // TH3F* trans_histo = (TH3F*)new_efficiency_file.Get("eff_histo");
    // trans_histo->SetDirectory(0);
    // new_efficiency_file.Close();
    // TH3F* binned_pdf = ConvertTransHisto(trans_histo);


    TH3F* gsim_histo = Create3DHisto(gsim_dataset_);
    TH3F* evtgen_histo = Create3DHisto(evtgen_dataset_);
    TH3F* binned_pdf = ConvertDensityToHisto(kde);
    TH3F* simulated_histo = Create3DHisto(evtgen_dataset_);
    simulated_histo->Multiply(binned_pdf);
    simulated_histo->Scale(gsim_histo->GetSumOfWeights() / simulated_histo->GetSumOfWeights());

    RooDataHist roo_simulated_histo("simulated", "simulated",
                                    RooArgList(thetat_, thetab_, phit_), simulated_histo);
    RooDataHist roo_gsim_histo("gsim", "gsim",
                               RooArgList(thetat_, thetab_, phit_), gsim_histo);
    RooDataHist roo_evtgen_histo("evtgen", "evtgen",
                               RooArgList(thetat_, thetab_, phit_), evtgen_histo);

    for (auto&& var : vars_) {
        tools::PlotVar(*var, roo_simulated_histo, roo_gsim_histo, true, false);
        tools::PlotVar(*var, roo_evtgen_histo, roo_gsim_histo, true, false);
    }

    for (int i = 0; i < 3; i++) {
        for (int j = i + 1; j < 3; j++) {
            tools::PlotVars2D(*vars_[i], *vars_[j], roo_simulated_histo, roo_gsim_histo, "", constants::format);
            tools::PlotPull2D(*vars_[i], *vars_[j], roo_simulated_histo, roo_gsim_histo, "", constants::format, true);
        }
    }

    // Scale the histos so that they have efficiency averages on y (or z) axes.
    // We would have to be more clever if there was a different number of bins
    // on each axis. For now, we assume they are equal.
    assert(eff_histo->GetNbinsX() == eff_histo->GetNbinsY());
    assert(eff_histo->GetNbinsX() == eff_histo->GetNbinsZ());
    const double normalization = 1. / eff_histo->GetNbinsX();

    eff_histo->Scale(normalization);
    RooDataHist roo_eff_histo_2D("eff", "eff",
                                 RooArgList(thetat_, thetab_, phit_), eff_histo);

    TH3F* scaled_binned_pdf = dynamic_cast<TH3F*>(binned_pdf->Clone("scaled_binned_pdf"));
    double scale = eff_histo->GetSumOfWeights() / binned_pdf->GetSumOfWeights();
    scaled_binned_pdf->Scale(scale);
    RooDataHist roo_eff_pdf_histo_2D("eff_pdf", "eff_pdf",
                                     RooArgList(thetat_, thetab_, phit_), scaled_binned_pdf);

    eff_histo->Scale(normalization);
    RooDataHist roo_eff_histo_1D("eff1D", "eff1D",
                                 RooArgList(thetat_, thetab_, phit_), eff_histo);
    scaled_binned_pdf->Scale(normalization);
    RooDataHist roo_eff_pdf_histo_1D("eff_pdf1D", "eff_pdf1D",
                                     RooArgList(thetat_, thetab_, phit_), scaled_binned_pdf);

    for (auto&& var : vars_) {
        tools::PlotVar(*var, roo_eff_histo_1D, roo_eff_pdf_histo_1D, true, false);
    }

    for (int i = 0; i < 3; i++) {
        for (int j = i + 1; j < 3; j++) {
            tools::PlotVars2D(*vars_[i], *vars_[j], roo_eff_histo_2D, roo_eff_pdf_histo_2D, "", constants::format);
            tools::PlotPull2D(*vars_[i], *vars_[j], roo_eff_histo_2D, roo_eff_pdf_histo_2D, "", constants::format, true);
        }
    }
}

void Fitter::ProcessKDEEfficiency2(const char* efficiency_file,
                                  const std::array<double, 6> bin_kde_pars,
                                  const std::array<double, 6> ada_kde_pars,
                                  const double mirror_margin) {
    TH3F* gsim_histo = Create3DHisto(gsim_dataset_);
    TH3F* gsim_kde = GetKDEHisto(gsim_dataset_, bin_kde_pars);
    TH3F* evtgen_histo = Create3DHisto(evtgen_dataset_);
    TH3F* evtgen_kde = GetKDEHisto(evtgen_dataset_, bin_kde_pars);

    TH3F* eff_histo = dynamic_cast<TH3F*>(gsim_histo->Clone("eff_histo"));
    eff_histo->Divide(gsim_histo, evtgen_histo);  //, 1.0, 1.0, "B");

    TH3F* eff_kde = dynamic_cast<TH3F*>(gsim_kde->Clone("eff_kde"));
    eff_kde->Divide(gsim_kde, evtgen_kde);  //, 1.0, 1.0, "B");

    TFile f(efficiency_file, "RECREATE");
    eff_kde->Write();
    f.Close();

    RooDataHist roo_gsim_kde_histo("gsim_kde", "gsim_kde",
                                    RooArgList(thetat_, thetab_, phit_), gsim_kde);
    RooDataHist roo_gsim_histo("gsim", "gsim",
                               RooArgList(thetat_, thetab_, phit_), gsim_histo);

    RooDataHist roo_evtgen_kde_histo("evtgen_kde", "evtgen_kde",
                                     RooArgList(thetat_, thetab_, phit_), evtgen_kde);
    RooDataHist roo_evtgen_histo("evtgen", "evtgen",
                                 RooArgList(thetat_, thetab_, phit_), evtgen_histo);

    for (auto&& var : vars_) {
        tools::PlotVar(*var, roo_gsim_kde_histo, roo_gsim_histo, true, false);
        tools::PlotVar(*var, roo_evtgen_kde_histo, roo_evtgen_histo, true, false);
    }

    for (int i = 0; i < 3; i++) {
        for (int j = i + 1; j < 3; j++) {
            tools::PlotVars2D(*vars_[i], *vars_[j], roo_gsim_kde_histo, roo_gsim_histo, "", constants::format);
            tools::PlotPull2D(*vars_[i], *vars_[j], roo_gsim_kde_histo, roo_gsim_histo, "", constants::format, true);
            tools::PlotVars2D(*vars_[i], *vars_[j], roo_evtgen_kde_histo, roo_evtgen_histo, "", constants::format);
            tools::PlotPull2D(*vars_[i], *vars_[j], roo_evtgen_kde_histo, roo_evtgen_histo, "", constants::format, true);
        }
    }

    // Scale the histos so that they have efficiency averages on y (or z) axes.
    // We would have to be more clever if there was a different number of bins
    // on each axis. For now, we assume they are equal.
    assert(eff_histo->GetNbinsX() == eff_histo->GetNbinsY());
    assert(eff_histo->GetNbinsX() == eff_histo->GetNbinsZ());
    const double normalization = 1. / eff_histo->GetNbinsX();

    eff_histo->Scale(normalization);
    RooDataHist roo_eff_histo_2D("eff", "eff",
                                 RooArgList(thetat_, thetab_, phit_), eff_histo);

    eff_kde->Scale(normalization);
    RooDataHist roo_eff_kde_histo_2D("eff_kde", "eff_kde",
                                     RooArgList(thetat_, thetab_, phit_), eff_kde);

    eff_histo->Scale(normalization);
    RooDataHist roo_eff_histo_1D("eff1D", "eff1D",
                                 RooArgList(thetat_, thetab_, phit_), eff_histo);
    eff_kde->Scale(normalization);
    RooDataHist roo_eff_kde_histo_1D("eff_pdf1D", "eff_pdf1D",
                                     RooArgList(thetat_, thetab_, phit_), eff_kde);

    for (auto&& var : vars_) {
        tools::PlotVar(*var, roo_eff_histo_1D, roo_eff_kde_histo_1D, true, false);
    }

    for (int i = 0; i < 3; i++) {
        for (int j = i + 1; j < 3; j++) {
            tools::PlotVars2D(*vars_[i], *vars_[j], roo_eff_histo_2D, roo_eff_kde_histo_2D, "", constants::format);
            tools::PlotPull2D(*vars_[i], *vars_[j], roo_eff_histo_2D, roo_eff_kde_histo_2D, "", constants::format, true);
        }
    }
}

void Fitter::ProcessNormalizedEfficiency(const char* efficiency_file) {


    TH3F* eff_histo = GetBinned3DEfficiency();

    TH3F* gsim_histo = Create3DHisto(gsim_dataset_);
    TH3F* evtgen_histo = Create3DHisto(evtgen_dataset_);
    TH3F* binned_pdf = NormalizePDF(eff_histo, 0, 1);

    TFile f(efficiency_file, "RECREATE");
    binned_pdf->Write();
    f.Close();
    Log::print(Log::info, "Efficiency written to '%s'\n", efficiency_file);

    TH3F* simulated_histo = Create3DHisto(evtgen_dataset_, "simulated_histo");
    simulated_histo->Multiply(binned_pdf);
    simulated_histo->Scale(gsim_histo->GetSumOfWeights() / simulated_histo->GetSumOfWeights());

    RooDataHist roo_simulated_histo("simulated", "simulated",
                                    RooArgList(thetat_, thetab_, phit_), simulated_histo);
    RooDataHist roo_gsim_histo("gsim", "gsim",
                               RooArgList(thetat_, thetab_, phit_), gsim_histo);
    RooDataHist roo_evtgen_histo("evtgen", "evtgen",
                               RooArgList(thetat_, thetab_, phit_), evtgen_histo);

    for (auto&& var : vars_) {
        tools::PlotVar(*var, roo_simulated_histo, roo_gsim_histo, true, false);
        tools::PlotVar(*var, roo_evtgen_histo, roo_gsim_histo, true, false);
    }

    for (int i = 0; i < 3; i++) {
        for (int j = i + 1; j < 3; j++) {
            tools::PlotVars2D(*vars_[i], *vars_[j], roo_simulated_histo, roo_gsim_histo, "", constants::format);
            tools::PlotPull2D(*vars_[i], *vars_[j], roo_simulated_histo, roo_gsim_histo, "", constants::format, true);
        }
    }

    // Scale the histos so that they have efficiency averages on y (or z) axes.
    // We would have to be more clever if there was a different number of bins
    // on each axis. For now, we assume they are equal.
    assert(eff_histo->GetNbinsX() == eff_histo->GetNbinsY());
    assert(eff_histo->GetNbinsX() == eff_histo->GetNbinsZ());
    const double normalization = 1. / eff_histo->GetNbinsX();

    eff_histo->Scale(normalization);
    RooDataHist roo_eff_histo_2D("eff", "eff",
                                 RooArgList(thetat_, thetab_, phit_), eff_histo);

    TH3F* scaled_binned_pdf = dynamic_cast<TH3F*>(binned_pdf->Clone("scaled_binned_pdf"));
    double scale = eff_histo->GetSumOfWeights() / binned_pdf->GetSumOfWeights();
    scaled_binned_pdf->Scale(scale);
    RooDataHist roo_eff_pdf_histo_2D("eff_pdf", "eff_pdf",
                                     RooArgList(thetat_, thetab_, phit_), scaled_binned_pdf);

    eff_histo->Scale(normalization);
    RooDataHist roo_eff_histo_1D("eff1D", "eff1D",
                                 RooArgList(thetat_, thetab_, phit_), eff_histo);
    scaled_binned_pdf->Scale(normalization);
    RooDataHist roo_eff_pdf_histo_1D("eff_pdf1D", "eff_pdf1D",
                                     RooArgList(thetat_, thetab_, phit_), scaled_binned_pdf);

    for (auto&& var : vars_) {
        tools::PlotVar(*var, roo_eff_histo_1D, roo_eff_pdf_histo_1D, false, true);
    }

    for (int i = 0; i < 3; i++) {
        for (int j = i + 1; j < 3; j++) {
            tools::PlotVars2D(*vars_[i], *vars_[j], roo_eff_histo_2D, roo_eff_pdf_histo_2D, "", constants::format);
            tools::PlotPull2D(*vars_[i], *vars_[j], roo_eff_histo_2D, roo_eff_pdf_histo_2D, "", constants::format, true);
        }
    }

}



TH3F* Fitter::GetBinned3DEfficiency() {
    TH3F* evtgen_histo = Create3DHisto(evtgen_dataset_);
    TH3F* gsim_histo = Create3DHisto(gsim_dataset_);

    TH3F* eff_histo = dynamic_cast<TH3F*>(gsim_histo->Clone("eff_histo"));
    eff_histo->Divide(gsim_histo, evtgen_histo);  //, 1.0, 1.0, "B");

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

TH3F* Fitter::Create3DHisto(const RooDataSet* dataset, const char* name) const {
    const int num_bins = 50;
    TH3F* histo = new TH3F(name ? name : dataset->GetName(), name ? name : dataset->GetTitle(),
                           num_bins, thetat_.getMin(), thetat_.getMax(), num_bins, thetab_.getMin(),
                           thetab_.getMax(), num_bins, phit_.getMin(), phit_.getMax());
    for (int i = 0; i < dataset->numEntries(); i++) {
        const RooArgSet* row = dataset->get(i);
        double thetat = row->getRealValue("thetat");
        double thetab = row->getRealValue("thetab");
        double phit = row->getRealValue("phit");
        histo->Fill(thetat, thetab, phit);
    }

    return histo;
}

TH3F* Fitter::ConvertDensityToHisto(AdaptiveKernelDensity pdf) const {
    const int num_bins = 50;
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

TH3F* Fitter::ConvertDensityToHisto(BinnedKernelDensity pdf) const {
    const int num_bins = 50;
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

TH3F* Fitter::ConvertTransHisto(TH3F* trans_histo) {
    const int num_bins = 50;
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

                double eff = 0;
                double transtht = thetat / TMath::Pi() * 2 - 1;
                double transthb = thetab / TMath::Pi() * 2 - 1;
                if (CanUseInterpolation(trans_histo, phit, transtht, transthb)) {
                    eff = trans_histo->Interpolate(phit, transtht, transthb);
                } else {
                    int binx = trans_histo->GetXaxis()->FindBin(phit);
                    int biny = trans_histo->GetYaxis()->FindBin(transtht);
                    int binz = trans_histo->GetZaxis()->FindBin(transthb);
                    int bin = trans_histo->GetBin(binx, biny, binz);
                    eff = trans_histo->GetBinContent(bin);
                }

                pdf_histo->SetBinContent(pdf_histo->GetBin(x, y, z), eff);
            }
        }
    }

    return pdf_histo;
}

/**
 * Check whether any of the vars is too close to the histogram edge to use
 * interpolation.
 * 
 * From the ROOT documentation: The given values (x,y,z) must be between first
 * bin center and last bin center for each coordinate.
 */
bool Fitter::CanUseInterpolation(const TH3F* histo, const double& phit, const double& transtht,
                               const double& transthb) const {
    double vars[3] = {phit, transtht, transthb};
    const TAxis* axes[3] = {histo->GetXaxis(), histo->GetYaxis(), histo->GetZaxis()};

    for (int i = 0; i < 3; i++) {
        int last_bin = axes[i]->GetNbins();
        double low_center = axes[i]->GetBinCenter(1);
        double high_center = axes[i]->GetBinCenter(last_bin);
        if (vars[i] < low_center || vars[i] > high_center) {
            return false;
        }
    }
    return true;
}

TTree* Fitter::DoubleTreeToFloatTree(TTree* double_tree) const {
    double dthetat, dthetab, dphit;
    double_tree->SetBranchAddress("thetat", &dthetat);
    double_tree->SetBranchAddress("thetab", &dthetab);
    double_tree->SetBranchAddress("phit", &dphit);

    TTree* tree = new TTree("floattree", "floattree");
    float thetat;
    float thetab;
    float phit;
    tree->Branch("thetat", &thetat, "thetat/F");
    tree->Branch("thetab", &thetab, "thetab/F");
    tree->Branch("phit", &phit, "phit/F");

    Long64_t num_entries = double_tree->GetEntries();
    for (Long64_t i = 0; i < num_entries; i++) {
        double_tree->GetEntry(i);
        thetat = (float)dthetat;
        thetab = (float)dthetab;
        phit = (float)dphit;
        tree->Fill();
    }

    return tree;
}

TH3F* Fitter::GetKDEHisto(RooDataSet* dataset, const std::array<double, 6> bin_kde_pars) const {
    OneDimPhaseSpace phasespace_thetat("phasespace_thetat", thetat_.getMin(), thetat_.getMax());
    OneDimPhaseSpace phasespace_thetab("phasespace_thetab", thetab_.getMin(), thetab_.getMax());
    OneDimPhaseSpace phasespace_phit("phasespace_phit", phit_.getMin(), phit_.getMax());
    CombinedPhaseSpace phasespace("phasespace", &phasespace_thetat, &phasespace_thetab,
                                  &phasespace_phit);
    
    // This line needs to be here to make the new dataset have TTree as a storage backend
    RooAbsData::setDefaultStorageType(RooAbsData::Tree);
    RooDataSet* dataset_with_tree = new RooDataSet("dataset_with_tree", "dataset_with_tree", dataset, *dataset->get()); 
    // TTree* double_tree = (TTree*)((RooTreeDataStore*)dataset_with_tree->store())->tree();
    RooAbsDataStore* ds = dataset_with_tree->store();
    TTree* double_tree = const_cast<TTree*>(ds->tree());

    // KernelDensity takes float TTrees not double TTrees
    TTree* tree = DoubleTreeToFloatTree(double_tree);

    BinnedKernelDensity kde("gsim_kde", &phasespace, tree, "thetat", "thetab", "phit",
                                "weight", bin_kde_pars[0], bin_kde_pars[1], bin_kde_pars[2],
                                bin_kde_pars[3], bin_kde_pars[4], bin_kde_pars[5], 0);

    // AdaptiveKernelDensity kde("KernelPDF", &phasespace, eff_tree, "thetat", "thetab", "phit",
    //                           "weight", ada_kde_pars[0], ada_kde_pars[1], ada_kde_pars[2],
    //                           ada_kde_pars[3], ada_kde_pars[4], ada_kde_pars[5], &bin_kde);

    TH3F* pdf = ConvertDensityToHisto(kde);
    TH3F* histo = Create3DHisto(dataset);
    pdf->Scale(histo->Integral()/pdf->Integral());

    return pdf;
}

TH3F* Fitter::NormalizePDF(const TH3F* pdf, const double low, const double high) {
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
                    // Log::print(Log::debug, "bin %i: value = %f, interpolation = %f\n", bin, pdf->GetBinContent(bin), interpolation);
                    normalized_pdf->SetBinContent(bin, interpolation);
                }
            }
        }
    }

    Log::print(Log::info, "Fixed %i/%i (%.2f%%) bins whose contents were not in [%f, %f].\n",
               fixed_bins, total_bins, (double)fixed_bins / total_bins * 100, low, high);
    return normalized_pdf;
}

double Fitter::Interpolate(const TH3F* histo, int x_org, int y_org, int z_org, int size) {
    double new_value = 0;
    int points = 0;
    int num_bins_x = histo->GetXaxis()->GetNbins();
    int num_bins_y = histo->GetYaxis()->GetNbins();
    int num_bins_z = histo->GetZaxis()->GetNbins();

    int skipped_bins = 0;
    int all_bins = 0;
    for (int x = x_org - size; x <= x_org + size; x++) {
        for (int y = y_org - size; y <= y_org + size; y++) {
            for (int z = z_org - size; z <= z_org + size; z++) {
                all_bins++;
                if (x == x_org && y == y_org && z == z_org) {
                    skipped_bins++;
                    continue;
                }
                if (x < 1 || y < 1 || z < 1 || x > num_bins_x || y > num_bins_y || z > num_bins_z) {
                    skipped_bins++;
                    continue;
                }
                int bin = histo->GetBin(x, y, z);
                new_value += histo->GetBinContent(bin);
                // Log::print(Log::debug, "delta = %f\n", histo->GetBinContent(bin));
                points++;
            }
        }
    }
    // if (skipped_bins) {
    //     Log::print(Log::debug, "Skipped %i/%i (%f %%) bins in Interpolate()\n", skipped_bins,
    //                all_bins, (double)skipped_bins / all_bins * 100);
    // }
    if (points == 0) return 0;
    return new_value/points;
}
