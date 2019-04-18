/**
 *  @file    fittercpv.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2016-11-03
 *
 *  @brief Class that performs the CP fit itself as well as plotting
 *
 */

#include "fittercpv.h"

// Standard includes
#include <unistd.h>
#include <array>
#include <boost/filesystem.hpp>
#include <ctime>
#include <fstream>
#include <regex>
#include <sstream>
#include <string>

// Belle includes
#include "tatami/tatami.h"

// ROOT includes
#include "RVersion.h"
#include "RooArgSet.h"
#include "RooBifurGauss.h"
#include "RooCategory.h"
#include "RooChi2Var.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooHist.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooProdPdf.h"
#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooTFnBinding.h"
#include "RooTreeDataStore.h"
#include "RooVoigtian.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TEnv.h"
#include "TF3.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TTree.h"

// Meerkat includes
#include "RooMeerkatPdf.hh"

// Local includes
#include "angularpdf.h"
#include "cksum.h"
#include "constants.h"
#include "dtcppdf.h"
#include "dtscfpdf.h"
#include "gitversion.h"
#include "log.h"
#include "rapidjson/document.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/stringbuffer.h"
#include "RooHistPdfFast.h"
#include "tools.h"

FitterCPV::FitterCPV() {

    num_CPUs_ = 1;
    do_lifetime_fit_ = false;
    do_mixing_fit_ = false;
    make_plots_ = false;
    perfect_tagging_ = false;
    generator_level_ = false;

    vrzerr_ = nullptr;
    vtzerr_ = nullptr;
}

FitterCPV::~FitterCPV() {
    if (output_file_) {
        if (output_file_->IsOpen()) {
            output_file_->Close();
        }
    }
}

void FitterCPV::InitVars(std::array<double, 16> par_input) {
    ap_ = new RooRealVar("ap", "|a_{#parallel}|", par_input[0], 0, 0.5);
    apa_ = new RooRealVar("apa", "arg(a_{#parallel})", par_input[1], 0, 1);
    a0_ = new RooRealVar("a0", "|a_{0}|", par_input[2], 0.8, 1);
    a0a_ = new RooRealVar("a0a", "arg(a_{0})", 0);
    at_ = new RooFormulaVar("at", "sqrt(1-ap*ap-a0*a0)", RooArgSet(*ap_, *a0_));
    ata_ = new RooRealVar("ata", "arg(a_{#perp})", par_input[3], 2, 4);

    xp_ = new RooRealVar("xp", "x_{#parallel}", par_input[4], -0.4, 0.4);
    x0_ = new RooRealVar("x0", "x_{0}", par_input[5], -0.4, 0.4);
    xt_ = new RooRealVar("xt", "x_{#perp}", par_input[6], -0.4, 0.4);

    yp_ = new RooRealVar("yp", "y_{#parallel}", par_input[7], -0.4, 0.4);
    y0_ = new RooRealVar("y0", "y_{0}", par_input[8], -0.4, 0.4);
    yt_ = new RooRealVar("yt", "y_{#perp}", par_input[9], -0.4, 0.4);

    xpb_ = new RooRealVar("xpb", "#bar x_{#parallel}", par_input[10], -0.4, 0.4);
    x0b_ = new RooRealVar("x0b", "#bar x_{0}", par_input[11], -0.4, 0.4);
    xtb_ = new RooRealVar("xtb", "#bar x_{#perp}", par_input[12], -0.4, 0.4);

    ypb_ = new RooRealVar("ypb", "#bar y_{#parallel}", par_input[13], -0.4, 0.4);
    y0b_ = new RooRealVar("y0b", "#bar y_{0}", par_input[14], -0.4, 0.4);
    ytb_ = new RooRealVar("ytb", "#bar y_{#perp}", par_input[15], -0.4, 0.4);

    dt_ = new RooRealVar("dt", "#Deltat [ps]", constants::cuts::dt_low, constants::cuts::dt_high);

    if (generator_level_) {
        thetat_ = new RooRealVar("thetatg", "#theta_{t}^{g} [rad]", 0, TMath::Pi());
        thetab_ = new RooRealVar("thetabg", "#theta_{b}^{g} [rad]", constants::cuts::thetab_low, constants::cuts::thetab_high);
        phit_ = new RooRealVar("phitg", "#phi_{t}^{g} [rad]", -TMath::Pi(), TMath::Pi());
    } else {
        thetat_ = new RooRealVar("thetat", "#theta_{t} [rad]", 0, TMath::Pi());
        thetab_ = new RooRealVar("thetab", "#theta_{b} [rad]", constants::cuts::thetab_low, constants::cuts::thetab_high);
        phit_ = new RooRealVar("phit", "#phi_{t} [rad]", -TMath::Pi(), TMath::Pi());
    }

    vrusable_ = new RooRealVar("vrusable", "vrusable", 0, 1);
    vrvtxz_ = new RooRealVar("vrvtxz", "vrvtxz", -10, 10);
    vrerr6_ = new RooRealVar("vrerr6", "vrerr6", -1, 1);
    vrchi2_ = new RooRealVar("vrchi2", "vrchi2", 0, 10000000);
    //  vreffxi_ = new RooRealVar( "vreffxi", "vreffxi", 0, 10000000);
    vrndf_ = new RooRealVar("vrndf", "vrndf", 0, 100);
    //  vreffndf_ = new RooRealVar( "vreffndf", "vreffndf", 0, 100);
    vrntrk_ = new RooRealVar("vrntrk", "vrntrk", 0, 100);

    vtusable_ = new RooRealVar("vtusable", "vtusable", 0, 1);
    vtvtxz_ = new RooRealVar("vtvtxz", "vtvtxz", -10, 10);
    vtndf_ = new RooRealVar("vtndf", "vtndf", 4);
    vterr6_ = new RooRealVar("vterr6", "vterr6", -1, 1);
    vtchi2_ = new RooRealVar("vtchi2", "vtchi2", 0, 10000000);
    vtndf_ = new RooRealVar("vtndf", "vtndf", 0, 100);
    vtntrk_ = new RooRealVar("vtntrk", "vtntrk", 0, 100);
    vtistagl_ = new RooRealVar("vtistagl", "vtistagl", 0, 100);

    expno_ = new RooRealVar("expno", "expno", 0, 100);
    expmc_ = new RooRealVar("expmc", "expmc", 0, 10);

    evmcflag_ = new RooRealVar("evmcflag", "evmcflag", 0, 10);
    brecflav_ = new RooRealVar("brecflav", "brecflav", -1, 1);
    btagmcli_ = new RooRealVar("btagmcli", "btagmcli", -1000, 1000);
    tagqr_ = new RooRealVar("tagqr", "tagqr", -1, 1);
    tagwtag_ = new RooRealVar("tagwtag", "tagwtag", 0, 0.5);

    benergy_ = new RooRealVar("benergy", "benergy", 0, 100);
    mbc_ = new RooRealVar("mbc", "mbc", 5, 6);

    shcosthb_ = new RooRealVar("shcosthb", "shcosthb", -1, 1);

    vrzerr_formula_ =
        new RooFormulaVar("vrzerr", "#sigma z_{rec} [cm]", "sqrt(vrerr6)", RooArgSet(*vrerr6_));
    vtzerr_formula_ =
        new RooFormulaVar("vtzerr", "#sigma z_{tag} [cm]", "sqrt(vterr6)", RooArgSet(*vterr6_));

    tau_ = new RooRealVar("tau", "#tau", constants::tau - 1, constants::tau + 1);
    dm_ = new RooRealVar("dm", "#Deltam", constants::dm - 1, constants::dm + 1);

    decaytype_ = new RooCategory("decaytype", "decaytype");
    decaytype_->defineType("a", 1);
    decaytype_->defineType("ab", 2);
    decaytype_->defineType("b", 3);
    decaytype_->defineType("bb", 4);

    // Make a copy of the input parameters for saving results, etc.
    par_input_ = par_input;

    PrepareVarArgsets();

    scf_angular_pdf_ = CreateAngularSCFPDF();
    bkg_angular_pdf_ = CreateAngularBKGPDF();
}
/**
 * Populate various vectors and argsets with all parameters
 * This populates conditional_vars_, dataset_vars_, parameters_ and their pair
 * argsets which are used elsewhere.
 */
void FitterCPV::PrepareVarArgsets() {
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
    dataset_vars_.push_back(&dt_);
    dataset_vars_.push_back(&thetat_);
    dataset_vars_.push_back(&thetab_);
    dataset_vars_.push_back(&phit_);

    parameters_.push_back(&ap_);
    parameters_.push_back(&apa_);
    parameters_.push_back(&a0_);
    parameters_.push_back(&ata_);
    parameters_.push_back(&xp_);
    parameters_.push_back(&x0_);
    parameters_.push_back(&xt_);
    parameters_.push_back(&yp_);
    parameters_.push_back(&y0_);
    parameters_.push_back(&yt_);
    parameters_.push_back(&xpb_);
    parameters_.push_back(&x0b_);
    parameters_.push_back(&xtb_);
    parameters_.push_back(&ypb_);
    parameters_.push_back(&y0b_);
    parameters_.push_back(&ytb_);

    // The variables present in the ntuple are added to an RooArgSet that will be needed
    // when we create a RooDataSet from the input_tree
    for (auto var : conditional_vars_) {
        conditional_vars_argset_.add(**var);
    }

    // TODO: Comments
    for (auto var : dataset_vars_) {
        dataset_vars_argset_.add(**var);
    }

    for (auto par : parameters_) {
        parameters_argset_.add(**par);
    }
}

/**
 * Make a simple plot of a variable and save it to a file
 * NOTE: @var can't be const (without const_cast) because of dumb RooFit design
 */
void FitterCPV::PlotVar(RooRealVar& var, const RooAbsData& data) const {
    TCanvas* canvas = new TCanvas(var.GetName(), var.GetTitle(), 500, 500);

    double range_min;
    double range_max;
    data.getRange(var, range_min, range_max, 0.05);

    RooPlot* plot = var.frame(range_min, range_max);

    data.plotOn(plot);
    plot->Draw();

    plot->SetTitle("");
    plot->GetXaxis()->SetTitle(TString(var.GetTitle()));
    plot->GetYaxis()->SetTitle("");

    canvas->Write();
    canvas->SaveAs(constants::format);
}

/**
 * Make a simple plot of a variable and save it to a file
 * NOTE: @var can't be const (without const_cast) because of dumb RooFit design
 */
void FitterCPV::PlotVar(RooRealVar& var, const RooAbsPdf& pdf) const {
    TString name = pdf.GetName();
    name += "_";
    name += var.GetName();
    TCanvas canvas(name, name, 500, 500);

    RooPlot* plot = var.frame();

    // From RooFit manual: "No variables are projected by default when PDF is
    // plotted on an empty frame" One has to be careful to explicitly project
    // over intended variables in this case, otherwise they are only slices.
    plot->updateNormVars(dataset_vars_argset_);

    pdf.plotOn(plot);
    plot->Draw();

    plot->SetTitle("");
    plot->GetXaxis()->SetTitle(TString(var.GetTitle()));
    plot->GetYaxis()->SetTitle("");

    canvas.Write();
    canvas.SaveAs(constants::format);
}

// TODO: Remove/refactor
void FitterCPV::FitCR() {
    RooDataSet* temp_dataset = static_cast<RooDataSet*>(dataset_->reduce("evmcflag==1"));
    dataset_ = temp_dataset;

    DtCPPDF cr_pdf_a("cr_pdf_a", "cr_pdf_a", false, true, perfect_tagging_, efficiency_model_, efficiency_files_,
                     *thetat_, *thetab_, *phit_, *ap_, *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_,
                     *y0_, *yt_, *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_,
                     *benergy_, *mbc_, *vrntrk_, *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_,
                     *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF cr_pdf_ab("cr_pdf_ab", "cr_pdf_ab", true, true, perfect_tagging_, efficiency_model_, efficiency_files_,
                      *thetat_, *thetab_, *phit_, *ap_, *apa_, *a0_, *ata_, *xpb_, *x0b_, *xtb_,
                      *ypb_, *y0b_, *ytb_, *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_,
                      *shcosthb_, *benergy_, *mbc_, *vrntrk_, *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_,
                      *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF cr_pdf_b("cr_pdf_b", "cr_pdf_b", false, false, perfect_tagging_, efficiency_model_, efficiency_files_,
                     *thetat_, *thetab_, *phit_, *ap_, *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_,
                     *y0_, *yt_, *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_,
                     *benergy_, *mbc_, *vrntrk_, *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_,
                     *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF cr_pdf_bb("cr_pdf_bb", "cr_pdf_bb", true, false, perfect_tagging_, efficiency_model_, efficiency_files_,
                      *thetat_, *thetab_, *phit_, *ap_, *apa_, *a0_, *ata_, *xpb_, *x0b_, *xtb_,
                      *ypb_, *y0b_, *ytb_, *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_,
                      *shcosthb_, *benergy_, *mbc_, *vrntrk_, *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_,
                      *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    RooSimultaneous sim_pdf("sim_pdf", "sim_pdf", *decaytype_);
    sim_pdf.addPdf(cr_pdf_a, "a");
    sim_pdf.addPdf(cr_pdf_ab, "ab");
    sim_pdf.addPdf(cr_pdf_b, "b");
    sim_pdf.addPdf(cr_pdf_bb, "bb");

    tau_->setConstant(true);
    dm_->setConstant(true);

    if (do_mixing_fit_) {
        result_ = sim_pdf.fitTo(*dataset_, RooFit::ConditionalObservables(conditional_vars_argset_),
                                RooFit::Minimizer("Minuit2"),
                                RooFit::Hesse(false), RooFit::Minos(false), RooFit::Save(true),
                                RooFit::NumCPU(num_CPUs_));
        if (result_) {
            result_->Print();
        }

        if (make_plots_) {
            RooDataSet* dataset_a =
                static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::a"));
            PlotWithPull(*dt_, *dataset_a, cr_pdf_a);

            RooDataSet* dataset_b =
                static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::b"));
            PlotWithPull(*dt_, *dataset_b, cr_pdf_b);

            RooDataSet* dataset_ab =
                static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::ab"));
            PlotWithPull(*dt_, *dataset_ab, cr_pdf_ab);

            RooDataSet* dataset_bb =
                static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::bb"));
            PlotWithPull(*dt_, *dataset_bb, cr_pdf_bb);

            RooDataHist* cr_hist = cr_pdf_a.generateBinned(
                RooArgSet(*thetat_, *thetab_, *phit_), 1000, RooFit::ExpectedData(true));

            RooHistPdf cr_histpdf("cr_histpdf", "cr_histpdf", RooArgSet(*thetat_, *thetab_, *phit_),
                                  *cr_hist);

            PlotWithPull(*thetat_, *dataset_, cr_histpdf);
            PlotWithPull(*thetab_, *dataset_, cr_histpdf);
            PlotWithPull(*phit_, *dataset_, cr_histpdf);
        }
    }
}

// TODO: Remove/refactor
void FitterCPV::FitCRSCF() {
    // RooDataSet* temp_dataset = static_cast<RooDataSet*>(dataset_->reduce("evmcflag!=1"));
    // dataset_ = temp_dataset;

    DtCPPDF cr_pdf_a("cr_pdf_a", "cr_pdf_a", false, true, perfect_tagging_, efficiency_model_, efficiency_files_,
                     *thetat_, *thetab_, *phit_, *ap_, *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_,
                     *y0_, *yt_, *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_,
                     *benergy_, *mbc_, *vrntrk_, *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_,
                     *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF cr_pdf_ab("cr_pdf_ab", "cr_pdf_ab", true, true, perfect_tagging_, efficiency_model_, efficiency_files_,
                      *thetat_, *thetab_, *phit_, *ap_, *apa_, *a0_, *ata_, *xpb_, *x0b_, *xtb_,
                      *ypb_, *y0b_, *ytb_, *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_,
                      *shcosthb_, *benergy_, *mbc_, *vrntrk_, *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_,
                      *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF cr_pdf_b("cr_pdf_b", "cr_pdf_b", false, false, perfect_tagging_, efficiency_model_, efficiency_files_,
                     *thetat_, *thetab_, *phit_, *ap_, *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_,
                     *y0_, *yt_, *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_,
                     *benergy_, *mbc_, *vrntrk_, *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_,
                     *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF cr_pdf_bb("cr_pdf_bb", "cr_pdf_bb", true, false, perfect_tagging_, efficiency_model_, efficiency_files_,
                      *thetat_, *thetab_, *phit_, *ap_, *apa_, *a0_, *ata_, *xpb_, *x0b_, *xtb_,
                      *ypb_, *y0b_, *ytb_, *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_,
                      *shcosthb_, *benergy_, *mbc_, *vrntrk_, *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_,
                      *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtSCFPDF scf_dt_pdf_a("cr_pdf_a", "cr_pdf_a", false, true, perfect_tagging_, *ap_, *apa_, *a0_,
                          *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_, *tagwtag_, *dt_, *tau_, *dm_,
                          *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_, *vrzerr_,
                          *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtSCFPDF scf_dt_pdf_ab("cr_pdf_ab", "cr_pdf_ab", true, true, perfect_tagging_, *ap_, *apa_,
                           *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_, *tagwtag_, *dt_, *tau_,
                           *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_, *vrzerr_,
                           *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtSCFPDF scf_dt_pdf_b("cr_pdf_b", "cr_pdf_b", false, false, perfect_tagging_, *ap_, *apa_, *a0_,
                          *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_, *tagwtag_, *dt_, *tau_, *dm_,
                          *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_, *vrzerr_,
                          *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtSCFPDF scf_dt_pdf_bb("cr_pdf_bb", "cr_pdf_bb", true, false, perfect_tagging_, *ap_, *apa_,
                           *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_, *tagwtag_, *dt_, *tau_,
                           *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_, *vrzerr_,
                           *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    // Self-cross-feed dt CF model
    RooRealVar scf_dt_cf_voigt_mu("scf_dt_cf_voigt_mu", "v_{#mu}", -0.249);
    RooRealVar scf_dt_cf_voigt_sigma("scf_dt_cf_voigt_sigma_", "v_{#sigma}", 1.910);
    RooRealVar scf_dt_cf_voigt_width("scf_dt_cf_voigt_width_", "v_{w}", 0.684);
    RooVoigtian scf_dt_cf_voigt("scf_dt_cf_voigt", "scf_dt_cf_voigt", *dt_, scf_dt_cf_voigt_mu,
                                scf_dt_cf_voigt_width, scf_dt_cf_voigt_sigma);

    RooRealVar scf_dt_cf_gaus_mu("scf_dt_cf_gaus_mu", "g_{#mu}", -0.105);
    RooRealVar scf_dt_cf_gaus_sigma("scf_dt_cf_gaus_sigma_", "g_{#sigma}", 0.881);
    RooGaussian scf_dt_cf_gaus("scf_dt_cf_gaus", "scf_dt_cf_gaus", *dt_, scf_dt_cf_gaus_mu,
                               scf_dt_cf_gaus_sigma);

    RooRealVar scf_dt_cf_f("scf_dt_cf_f", "f_{v/g}", 0.720);
    RooAddPdf scf_dt_cf_model("scf_dt_cf_model", "scf_dt_cf_model",
                              RooArgList(scf_dt_cf_voigt, scf_dt_cf_gaus), RooArgList(scf_dt_cf_f));

    // Self-cross-feed dt DCS model
    RooRealVar scf_dt_dcs_voigt_mu("scf_dt_dcs_voigt_mu", "v_{#mu}", -0.303);
    RooRealVar scf_dt_dcs_voigt_sigma("scf_dt_dcs_voigt_sigma_", "v_{#sigma}", 2.739);
    RooRealVar scf_dt_dcs_voigt_width("scf_dt_dcs_voigt_width_", "v_{w}", 0.712);
    RooVoigtian scf_dt_dcs_voigt("scf_dt_dcs_voigt", "scf_dt_dcs_voigt", *dt_, scf_dt_dcs_voigt_mu,
                                 scf_dt_dcs_voigt_width, scf_dt_dcs_voigt_sigma);

    RooRealVar scf_dt_dcs_gaus_mu("scf_dt_dcs_gaus_mu", "g_{#mu}", -0.177);
    RooRealVar scf_dt_dcs_gaus_sigma("scf_dt_dcs_gaus_sigma_", "g_{#sigma}", 1.156);
    RooGaussian scf_dt_dcs_gaus("scf_dt_dcs_gaus", "scf_dt_dcs_gaus", *dt_, scf_dt_dcs_gaus_mu,
                                scf_dt_dcs_gaus_sigma);

    RooRealVar scf_dt_dcs_f("scf_dt_dcs_f", "f_{v/g}", 0.712);
    RooAddPdf scf_dt_dcs_model("scf_dt_dcs_model", "scf_dt_dcs_model",
                               RooArgList(scf_dt_dcs_voigt, scf_dt_dcs_gaus),
                               RooArgList(scf_dt_dcs_f));

    RooProdPdf scf_pdf_a("scf_pdf_a", "scf_pdf_a", RooArgList(scf_dt_cf_model, *scf_angular_pdf_));
    RooProdPdf scf_pdf_ab("scf_pdf_ab", "scf_pdf_ab", RooArgList(scf_dt_cf_model, *scf_angular_pdf_));
    RooProdPdf scf_pdf_b("scf_pdf_b", "scf_pdf_b", RooArgList(scf_dt_dcs_model, *scf_angular_pdf_));
    RooProdPdf scf_pdf_bb("scf_pdf_bb", "scf_pdf_bb", RooArgList(scf_dt_dcs_model, *scf_angular_pdf_));

    // RooProdPdf scf_pdf_a("scf_pdf_a", "scf_pdf_a", RooArgList(scf_dt_pdf_a, *scf_angular_pdf));
    // RooProdPdf scf_pdf_ab("scf_pdf_ab", "scf_pdf_ab", RooArgList(scf_dt_pdf_ab, *scf_angular_pdf));
    // RooProdPdf scf_pdf_b("scf_pdf_b", "scf_pdf_b", RooArgList(scf_dt_pdf_b, *scf_angular_pdf));
    // RooProdPdf scf_pdf_bb("scf_pdf_bb", "scf_pdf_bb", RooArgList(scf_dt_pdf_bb, *scf_angular_pdf));

    RooRealVar cr_scf_f("cr_scf_f", "f_{cr}", constants::fraction_cr_of_crscf, 0.80, 0.99);
    cr_scf_f.setConstant();

    RooAddPdf pdf_a("pdf_a", "pdf_a", RooArgList(cr_pdf_a, scf_pdf_a), RooArgList(cr_scf_f));
    RooAddPdf pdf_ab("pdf_ab", "pdf_ab", RooArgList(cr_pdf_ab, scf_pdf_ab), RooArgList(cr_scf_f));
    RooAddPdf pdf_b("pdf_b", "pdf_b", RooArgList(cr_pdf_b, scf_pdf_b), RooArgList(cr_scf_f));
    RooAddPdf pdf_bb("pdf_bb", "pdf_bb", RooArgList(cr_pdf_bb, scf_pdf_bb), RooArgList(cr_scf_f));

    RooSimultaneous sim_pdf("sim_pdf", "sim_pdf", *decaytype_);
    sim_pdf.addPdf(pdf_a, "a");
    sim_pdf.addPdf(pdf_ab, "ab");
    sim_pdf.addPdf(pdf_b, "b");
    sim_pdf.addPdf(pdf_bb, "bb");

    tau_->setConstant(true);
    dm_->setConstant(true);

    if (do_mixing_fit_) {
        result_ = sim_pdf.fitTo(*dataset_,
        RooFit::ConditionalObservables(conditional_vars_argset_),
                                RooFit::Minimizer("Minuit2"),
                                RooFit::Hesse(false), RooFit::Minos(false),
                                RooFit::Save(true), RooFit::NumCPU(num_CPUs_));
        if (result_) {
            result_->Print();
        }

        if (make_plots_) {
            // // PDF for B_bar differs, so we have to generate them separately. (One
            // // could also use *extended* RooSimultaneous or 'ProtoData' for
            // // RooSimultaneous.)
            // RooDataHist* cr_hist_a = cr_pdf_a.generateBinned(RooArgSet(*dt_, *thetat_, *thetab_, *phit_), 1000,
            //                                             RooFit::ExpectedData(true));
            // RooDataHist* cr_hist_B_bar = cr_pdf_B_bar.generateBinned(RooArgSet(*thetat_, *thetab_, *phit_),
            //                                                     1000, RooFit::ExpectedData(true));
            // // Add histos from both particle and anti-particle PDFs to create the
            // // final RooHistPdf.
            // cr_hist->add(*cr_hist_B_bar);
            // RooHistPdf cr_histpdf("cr_histpdf", "cr_histpdf", RooArgSet(*thetat_, *thetab_, *phit_),
            //                     *cr_hist);

            // RooDataHist* scf_hist = scf_pdf->generateBinned(RooArgSet(*thetat_, *thetab_, *phit_), 1000,
            //                                             RooFit::ExpectedData(true));
            // RooHistPdf scf_histpdf("scf_histpdf", "scf_histpdf", RooArgSet(*thetat_, *thetab_, *phit_),
            //                     *scf_hist);

            // RooRealVar cr_scf_f_gen("cr_scf_f_gen", "cr_scf_f_gen", (1 + cr_scf_f.getVal())/2);
            // RooAddPdf all_histpdf("all_histpdf", "all_histpdf", RooArgList(cr_histpdf, scf_histpdf),
            //                     cr_scf_f);

            // Set the current directory back to the one for plots (ugly ROOT stuff)
            if (output_file_) {
                output_file_->cd();
            }

            RooDataSet* dataset_a =
                static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::a"));
            PlotWithPull(*dt_, *dataset_a, cr_pdf_a);

            // RooDataSet* dataset_b =
            //     static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::b"));
            // PlotWithPull(*dt_, *dataset_b, cr_pdf_b);

            RooDataSet* dataset_ab =
                static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::ab"));
            PlotWithPull(*dt_, *dataset_ab, cr_pdf_ab);

            // RooDataSet* dataset_bb =
            //     static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::bb"));
            // PlotWithPull(*dt_, *dataset_bb, cr_pdf_bb);

            PlotWithPull(*thetat_, *dataset_, cr_pdf_a);
            PlotWithPull(*thetab_, *dataset_, cr_pdf_a);
            PlotWithPull(*phit_, *dataset_, cr_pdf_a);
        }
    }
}

void FitterCPV::FitAll() {
    DtCPPDF cr_pdf_a("cr_pdf_a", "cr_pdf_a", false, true, perfect_tagging_, efficiency_model_, efficiency_files_,
                     *thetat_, *thetab_, *phit_, *ap_, *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_,
                     *y0_, *yt_, *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_,
                     *benergy_, *mbc_, *vrntrk_, *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_,
                     *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF cr_pdf_ab("cr_pdf_ab", "cr_pdf_ab", true, true, perfect_tagging_, efficiency_model_, efficiency_files_,
                      *thetat_, *thetab_, *phit_, *ap_, *apa_, *a0_, *ata_, *xpb_, *x0b_, *xtb_,
                      *ypb_, *y0b_, *ytb_, *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_,
                      *shcosthb_, *benergy_, *mbc_, *vrntrk_, *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_,
                      *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF cr_pdf_b("cr_pdf_b", "cr_pdf_b", false, false, perfect_tagging_, efficiency_model_, efficiency_files_,
                     *thetat_, *thetab_, *phit_, *ap_, *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_,
                     *y0_, *yt_, *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_,
                     *benergy_, *mbc_, *vrntrk_, *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_,
                     *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF cr_pdf_bb("cr_pdf_bb", "cr_pdf_bb", true, false, perfect_tagging_, efficiency_model_, efficiency_files_,
                      *thetat_, *thetab_, *phit_, *ap_, *apa_, *a0_, *ata_, *xpb_, *x0b_, *xtb_,
                      *ypb_, *y0b_, *ytb_, *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_,
                      *shcosthb_, *benergy_, *mbc_, *vrntrk_, *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_,
                      *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtSCFPDF scf_dt_pdf_a("scf_dt_pdf_a", "scf_dt_pdf_a", false, true, perfect_tagging_, *ap_,
                          *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_, *tagwtag_, *dt_,
                          *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
                          *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_,
                          *vtistagl_);

    DtSCFPDF scf_dt_pdf_ab("scf_dt_pdf_ab", "scf_dt_pdf_ab", true, true, perfect_tagging_, *ap_,
                           *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_, *tagwtag_, *dt_,
                           *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
                           *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_,
                           *vtistagl_);

    DtSCFPDF scf_dt_pdf_b("scf_dt_pdf_b", "scf_dt_pdf_b", false, false, perfect_tagging_, *ap_,
                          *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_, *tagwtag_, *dt_,
                          *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
                          *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_,
                          *vtistagl_);

    DtSCFPDF scf_dt_pdf_bb("scf_dt_pdf_bb", "scf_dt_pdf_bb", true, false, perfect_tagging_, *ap_,
                           *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_, *tagwtag_, *dt_,
                           *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
                           *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_,
                           *vtistagl_);

    // Self-cross-feed dt CF model
    RooRealVar scf_dt_cf_voigt_mu("scf_dt_cf_voigt_mu", "v_{#mu}", -0.249);
    RooRealVar scf_dt_cf_voigt_sigma("scf_dt_cf_voigt_sigma_", "v_{#sigma}", 1.910);
    RooRealVar scf_dt_cf_voigt_width("scf_dt_cf_voigt_width_", "v_{w}", 0.684);
    RooVoigtian scf_dt_cf_voigt("scf_dt_cf_voigt", "scf_dt_cf_voigt", *dt_, scf_dt_cf_voigt_mu,
                                scf_dt_cf_voigt_width, scf_dt_cf_voigt_sigma);

    RooRealVar scf_dt_cf_gaus_mu("scf_dt_cf_gaus_mu", "g_{#mu}", -0.105);
    RooRealVar scf_dt_cf_gaus_sigma("scf_dt_cf_gaus_sigma_", "g_{#sigma}", 0.881);
    RooGaussian scf_dt_cf_gaus("scf_dt_cf_gaus", "scf_dt_cf_gaus", *dt_, scf_dt_cf_gaus_mu,
                               scf_dt_cf_gaus_sigma);

    RooRealVar scf_dt_cf_f("scf_dt_cf_f", "f_{v/g}", 0.720);
    RooAddPdf scf_dt_cf_model("scf_dt_cf_model", "scf_dt_cf_model",
                              RooArgList(scf_dt_cf_voigt, scf_dt_cf_gaus), RooArgList(scf_dt_cf_f));

    // Self-cross-feed dt DCS model
    RooRealVar scf_dt_dcs_voigt_mu("scf_dt_dcs_voigt_mu", "v_{#mu}", -0.303);
    RooRealVar scf_dt_dcs_voigt_sigma("scf_dt_dcs_voigt_sigma_", "v_{#sigma}", 2.739);
    RooRealVar scf_dt_dcs_voigt_width("scf_dt_dcs_voigt_width_", "v_{w}", 0.712);
    RooVoigtian scf_dt_dcs_voigt("scf_dt_dcs_voigt", "scf_dt_dcs_voigt", *dt_, scf_dt_dcs_voigt_mu,
                                 scf_dt_dcs_voigt_width, scf_dt_dcs_voigt_sigma);

    RooRealVar scf_dt_dcs_gaus_mu("scf_dt_dcs_gaus_mu", "g_{#mu}", -0.177);
    RooRealVar scf_dt_dcs_gaus_sigma("scf_dt_dcs_gaus_sigma_", "g_{#sigma}", 1.156);
    RooGaussian scf_dt_dcs_gaus("scf_dt_dcs_gaus", "scf_dt_dcs_gaus", *dt_, scf_dt_dcs_gaus_mu,
                                scf_dt_dcs_gaus_sigma);

    RooRealVar scf_dt_dcs_f("scf_dt_dcs_f", "f_{v/g}", 0.712);
    RooAddPdf scf_dt_dcs_model("scf_dt_dcs_model", "scf_dt_dcs_model",
                               RooArgList(scf_dt_dcs_voigt, scf_dt_dcs_gaus),
                               RooArgList(scf_dt_dcs_f));

    RooProdPdf scf_pdf_a("scf_pdf_a", "scf_pdf_a", RooArgList(scf_dt_cf_model, *scf_angular_pdf_));
    RooProdPdf scf_pdf_ab("scf_pdf_ab", "scf_pdf_ab", RooArgList(scf_dt_cf_model, *scf_angular_pdf_));
    RooProdPdf scf_pdf_b("scf_pdf_b", "scf_pdf_b", RooArgList(scf_dt_dcs_model, *scf_angular_pdf_));
    RooProdPdf scf_pdf_bb("scf_pdf_bb", "scf_pdf_bb", RooArgList(scf_dt_dcs_model, *scf_angular_pdf_));

    // RooProdPdf scf_pdf_a("scf_pdf_a", "scf_pdf_a", RooArgList(scf_dt_pdf_a, *scf_angular_pdf));
    // RooProdPdf scf_pdf_ab("scf_pdf_ab", "scf_pdf_ab", RooArgList(scf_dt_pdf_ab, *scf_angular_pdf));
    // RooProdPdf scf_pdf_b("scf_pdf_b", "scf_pdf_b", RooArgList(scf_dt_pdf_b, *scf_angular_pdf));
    // RooProdPdf scf_pdf_bb("scf_pdf_bb", "scf_pdf_bb", RooArgList(scf_dt_pdf_bb, *scf_angular_pdf));


    // Background dt CF model
    RooRealVar bkg_dt_cf_voigt_mu("bkg_dt_cf_voigt_mu", "v_{#mu}", -0.159);
    RooRealVar bkg_dt_cf_voigt_sigma("bkg_dt_cf_voigt_sigma_", "v_{#sigma}", 1.644);
    RooRealVar bkg_dt_cf_voigt_width("bkg_dt_cf_voigt_width_", "v_{w}", 0.874);
    RooVoigtian bkg_dt_cf_voigt("bkg_dt_cf_voigt", "bkg_dt_cf_voigt", *dt_, bkg_dt_cf_voigt_mu,
                                bkg_dt_cf_voigt_width, bkg_dt_cf_voigt_sigma);

    RooRealVar bkg_dt_cf_gaus_mu("bkg_dt_cf_gaus_mu", "g_{#mu}", -0.026);
    RooRealVar bkg_dt_cf_gaus_sigma("bkg_dt_cf_gaus_sigma_", "g_{#sigma}", 0.613);
    RooGaussian bkg_dt_cf_gaus("bkg_dt_cf_gaus", "bkg_dt_cf_gaus", *dt_, bkg_dt_cf_gaus_mu,
                               bkg_dt_cf_gaus_sigma);

    RooRealVar bkg_dt_cf_f("bkg_dt_cf_f", "f_{v/g}", 0.632);
    RooAddPdf bkg_dt_cf_model("bkg_dt_cf_model", "bkg_dt_cf_model",
                              RooArgList(bkg_dt_cf_voigt, bkg_dt_cf_gaus), RooArgList(bkg_dt_cf_f));

    // Background dt DCS model
    RooRealVar bkg_dt_dcs_voigt_mu("bkg_dt_dcs_voigt_mu", "v_{#mu}", -0.234);
    RooRealVar bkg_dt_dcs_voigt_sigma("bkg_dt_dcs_voigt_sigma_", "v_{#sigma}", 2.210);
    RooRealVar bkg_dt_dcs_voigt_width("bkg_dt_dcs_voigt_width_", "v_{w}", 1.143);
    RooVoigtian bkg_dt_dcs_voigt("bkg_dt_dcs_voigt", "bkg_dt_dcs_voigt", *dt_, bkg_dt_dcs_voigt_mu,
                                 bkg_dt_dcs_voigt_width, bkg_dt_dcs_voigt_sigma);

    RooRealVar bkg_dt_dcs_gaus_mu("bkg_dt_dcs_gaus_mu", "g_{#mu}", -0.020);
    RooRealVar bkg_dt_dcs_gaus_sigma("bkg_dt_dcs_gaus_sigma_", "g_{#sigma}", 0.639);
    RooGaussian bkg_dt_dcs_gaus("bkg_dt_dcs_gaus", "bkg_dt_dcs_gaus", *dt_, bkg_dt_dcs_gaus_mu,
                                bkg_dt_dcs_gaus_sigma);

    RooRealVar bkg_dt_dcs_f("bkg_dt_dcs_f", "f_{v/g}", 0.632);
    RooAddPdf bkg_dt_dcs_model("bkg_dt_dcs_model", "bkg_dt_dcs_model",
                               RooArgList(bkg_dt_dcs_voigt, bkg_dt_dcs_gaus),
                               RooArgList(bkg_dt_dcs_f));

    RooProdPdf bkg_pdf_a("bkg_pdf_a", "bkg_pdf_a", RooArgList(bkg_dt_cf_model, *bkg_angular_pdf_));
    RooProdPdf bkg_pdf_ab("bkg_pdf_ab", "bkg_pdf_ab", RooArgList(bkg_dt_cf_model, *bkg_angular_pdf_));
    RooProdPdf bkg_pdf_b("bkg_pdf_b", "bkg_pdf_b", RooArgList(bkg_dt_dcs_model, *bkg_angular_pdf_));
    RooProdPdf bkg_pdf_bb("bkg_pdf_bb", "bkg_pdf_bb", RooArgList(bkg_dt_dcs_model, *bkg_angular_pdf_));


    RooRealVar cr_f("cr_f", "f_{cr}", constants::fraction_cr_of_crscfbkg, 0.10, 0.99);
    RooRealVar scf_f("scf_f", "f_{scf}", constants::fraction_scf_of_crscfbkg, 0.10, 0.99);

    cr_f.setConstant();
    scf_f.setConstant();

    RooAddPdf pdf_a("pdf_a", "pdf_a", RooArgList(cr_pdf_a, scf_pdf_a, bkg_pdf_a), RooArgList(cr_f, scf_f));
    RooAddPdf pdf_ab("pdf_ab", "pdf_ab", RooArgList(cr_pdf_ab, scf_pdf_ab, bkg_pdf_a), RooArgList(cr_f, scf_f));
    RooAddPdf pdf_b("pdf_b", "pdf_b", RooArgList(cr_pdf_b, scf_pdf_b, bkg_pdf_a), RooArgList(cr_f, scf_f));
    RooAddPdf pdf_bb("pdf_bb", "pdf_bb", RooArgList(cr_pdf_bb, scf_pdf_bb, bkg_pdf_a), RooArgList(cr_f, scf_f));

    RooSimultaneous sim_pdf("sim_pdf", "sim_pdf", *decaytype_);
    sim_pdf.addPdf(pdf_a, "a");
    sim_pdf.addPdf(pdf_ab, "ab");
    sim_pdf.addPdf(pdf_b, "b");
    sim_pdf.addPdf(pdf_bb, "bb");

    tau_->setConstant(true);
    dm_->setConstant(true);

    if (do_mixing_fit_) {
        result_ =
            sim_pdf.fitTo(*dataset_, RooFit::ConditionalObservables(conditional_vars_argset_),
                          RooFit::Hesse(false), RooFit::Minos(false), RooFit::Minimizer("Minuit2"),
                          RooFit::Save(true), RooFit::NumCPU(num_CPUs_));
        if (result_) {
            result_->Print();
        }

        if (make_plots_) {
            // TODO: Refactor. Now! Really!
            // This plotting code actually works! In finite time! Have to start from here.
            // RooDataSet* cr_dataset_ = static_cast<RooDataSet*>(dataset_->reduce("evmcflag==1"));
            // RooDataSet* scf_dataset_ = static_cast<RooDataSet*>(dataset_->reduce("evmcflag!=1"));

            // RooDataSet* cr_dataset_a =
            //     static_cast<RooDataSet*>(cr_dataset_->reduce("decaytype==decaytype::a"));
            // PlotWithPull(*dt_, *cr_dataset_a, cr_pdf_a);

            // RooDataSet* cr_dataset_b =
            //     static_cast<RooDataSet*>(cr_dataset_->reduce("decaytype==decaytype::b"));
            // PlotWithPull(*dt_, *cr_dataset_b, cr_pdf_b);

            // RooDataSet* cr_dataset_ab =
            //     static_cast<RooDataSet*>(cr_dataset_->reduce("decaytype==decaytype::ab"));
            // PlotWithPull(*dt_, *cr_dataset_ab, cr_pdf_ab);

            // RooDataSet* cr_dataset_bb =
            //     static_cast<RooDataSet*>(cr_dataset_->reduce("decaytype==decaytype::bb"));
            // PlotWithPull(*dt_, *cr_dataset_bb, cr_pdf_bb);

            // RooDataSet* scf_dataset_a =
            //     static_cast<RooDataSet*>(scf_dataset_->reduce("decaytype==decaytype::a"));
            // PlotWithPull(*dt_, *scf_dataset_a, scf_pdf_a);

            // RooDataSet* scf_dataset_b =
            //     static_cast<RooDataSet*>(scf_dataset_->reduce("decaytype==decaytype::b"));
            // PlotWithPull(*dt_, *scf_dataset_b, scf_pdf_b);

            // RooDataSet* scf_dataset_ab =
            //     static_cast<RooDataSet*>(scf_dataset_->reduce("decaytype==decaytype::ab"));
            // PlotWithPull(*dt_, *scf_dataset_ab, scf_pdf_ab);

            // RooDataSet* scf_dataset_bb =
            //     static_cast<RooDataSet*>(scf_dataset_->reduce("decaytype==decaytype::bb"));
            // PlotWithPull(*dt_, *scf_dataset_bb, scf_pdf_bb);

            // thetab_->setBins(300);
            RooDataHist* cr_hist =
                cr_pdf_a.generateBinned(RooArgSet(*dt_, *thetat_, *thetab_, *phit_), 1000 * cr_f.getVal(),
                                        RooFit::ExpectedData(true));

            RooDataHist* scf_hist = scf_pdf_a.generateBinned(
                RooArgSet(*dt_, *thetat_, *thetab_, *phit_), 1000 * scf_f.getVal(),
                RooFit::ExpectedData(true));

            RooDataHist* bkg_hist = bkg_pdf_a.generateBinned(
                RooArgSet(*dt_, *thetat_, *thetab_, *phit_), 1000 * (1 - cr_f.getVal() - scf_f.getVal()),
                RooFit::ExpectedData(true));

            // RooHistPdf all_histpdf("all_histpdf", "all_histpdf", RooArgSet(*thetat_, *thetab_,
            // *phit_), *all_hist);
            RooHistPdf cr_histpdf("cr_histpdf", "cr_histpdf", RooArgSet(*dt_, *thetat_, *thetab_, *phit_),
                                  *cr_hist);
            RooHistPdf scf_histpdf("scf_histpdf", "scf_histpdf",
                                   RooArgSet(*dt_, *thetat_, *thetab_, *phit_), *scf_hist);
            RooHistPdf bkg_histpdf("bkg_histpdf", "bkg_histpdf",
                                   RooArgSet(*dt_, *thetat_, *thetab_, *phit_), *bkg_hist);
            RooAddPdf all_histpdf("all_histpdf", "all_histpdf", RooArgList(cr_histpdf, scf_histpdf, bkg_histpdf),
                                  RooArgList(cr_f, scf_f));

            std::vector<RooAbsPdf*> components;
            components.push_back(&cr_histpdf);
            components.push_back(&scf_histpdf);
            components.push_back(&bkg_histpdf);

            // thetab_->setBins(100);

            // Set the current directory back to the one for plots (ugly ROOT stuff)
            if (output_file_) {
                output_file_->cd();
            }

            PlotWithPull(*dt_, *dataset_, all_histpdf, components);
            PlotWithPull(*thetat_, *dataset_, all_histpdf, components);
            PlotWithPull(*thetab_, *dataset_, all_histpdf, components);
            PlotWithPull(*phit_, *dataset_, all_histpdf, components);
        }
    }
}

void FitterCPV::FitAngularCR() {
    RooDataSet* temp_dataset = static_cast<RooDataSet*>(dataset_->reduce("evmcflag==1"));
    dataset_ = temp_dataset;

    Log::print(Log::debug, "Fitting %i events.\n", dataset_->numEntries());

    AngularPDF pdf_B("pdf_B", "pdf_B", false, efficiency_model_, efficiency_files_, *thetat_, *thetab_, *phit_, *ap_,
                     *apa_, *a0_, *ata_);
    AngularPDF pdf_B_bar("pdf_B_bar", "pdf_B_bar", true, efficiency_model_, efficiency_files_, *thetat_, *thetab_,
                         *phit_, *ap_, *apa_, *a0_, *ata_);

    RooSimultaneous sim_pdf("sim_pdf", "sim_pdf", *decaytype_);
    sim_pdf.addPdf(pdf_B, "a");
    sim_pdf.addPdf(pdf_B_bar, "ab");
    sim_pdf.addPdf(pdf_B, "b");
    sim_pdf.addPdf(pdf_B_bar, "bb");

    result_ = sim_pdf.fitTo(*dataset_, RooFit::Minimizer("Minuit2"), RooFit::Hesse(false),
                            RooFit::Minos(false), RooFit::Save(true), RooFit::NumCPU(num_CPUs_));

    if (result_) {
        result_->Print();
    }

    if (make_plots_) {
        // PDF for B_bar differs, so we have to generate them separately. (One
        // could also use *extended* RooSimultaneous or 'ProtoData' for
        // RooSimultaneous.)
        RooDataHist* cr_hist = pdf_B.generateBinned(RooArgSet(*thetat_, *thetab_, *phit_), 1000,
                                                    RooFit::ExpectedData(true));
        RooDataHist* cr_hist_B_bar = pdf_B_bar.generateBinned(RooArgSet(*thetat_, *thetab_, *phit_),
                                                              1000, RooFit::ExpectedData(true));
        // Add histos from both particle and anti-particle PDFs to create the
        // final RooHistPdf.
        cr_hist->add(*cr_hist_B_bar);
        RooHistPdf cr_histpdf("cr_histpdf", "cr_histpdf", RooArgSet(*thetat_, *thetab_, *phit_),
                              *cr_hist);

        // Set the current directory back to the one for plots (ugly ROOT stuff)
        if (output_file_) {
            output_file_->cd();
        }

        // const double margin_ap = 0.005;
        // const double margin_apa = 0.25;
        // const double margin_a0 = 0.0015;
        // const double margin_ata = 0.3;
        // SaveLikelihoodScan(sim_pdf, ap_, margin_ap);
        // SaveLikelihoodScan(sim_pdf, apa_, margin_apa);
        // SaveLikelihoodScan(sim_pdf, a0_, margin_a0);
        // SaveLikelihoodScan(sim_pdf, ata_, margin_ata);
        // SaveLikelihoodScan(sim_pdf, ap_, apa_, margin_ap, margin_apa);
        // SaveLikelihoodScan(sim_pdf, ap_, a0_, margin_ap, margin_a0);
        // SaveLikelihoodScan(sim_pdf, ap_, ata_, margin_ap, margin_ata);
        // SaveLikelihoodScan(sim_pdf, apa_, a0_, margin_apa, margin_a0);
        // SaveLikelihoodScan(sim_pdf, apa_, ata_, margin_ap, margin_ata);
        // SaveLikelihoodScan(sim_pdf, a0_, ata_, margin_a0, margin_ata);

        // Plot data with PDF averaged over B and Bbar (just like in a standard dataset)
        PlotWithPull(*thetat_, *dataset_, cr_histpdf);
        PlotWithPull(*thetab_, *dataset_, cr_histpdf);
        PlotWithPull(*phit_, *dataset_, cr_histpdf);

        // Plot 1D PDF (just for a B not Bbar) projections
        // PlotVar(*thetat_, pdf_B);
        // PlotVar(*thetab_, pdf_B);
        // PlotVar(*phit_, pdf_B);

        // We need this to get the exact number of events to be generated for
        // the PDF distribution. This is unnecessary for the above, because we
        // create a RooHistPdf from it and that is normalized to the data when
        // plotted.
        RooDataSet* dataset_B = dynamic_cast<RooDataSet*>(
            dataset_->reduce("decaytype==decaytype::a||decaytype==decaytype::b"));
        RooDataSet* dataset_B_bar = dynamic_cast<RooDataSet*>(
            dataset_->reduce("decaytype==decaytype::ab||decaytype==decaytype::bb"));

        thetat_->setBins(40);
        thetab_->setBins(40);
        phit_->setBins(40);
        RooDataHist hist("hist", "hist", RooArgSet(*thetat_, *thetab_, *phit_), *dataset_);
        tools::PlotVars2D(*thetat_, *thetab_, hist);
        tools::PlotVars2D(*thetat_, *phit_, hist);
        tools::PlotVars2D(*thetab_, *phit_, hist);

        // We change the binning for 2D plots, so we have to generate new binned
        // dataset, to avoid dealing with rebinning.
        delete cr_hist;
        delete cr_hist_B_bar;
        cr_hist = pdf_B.generateBinned(RooArgSet(*thetat_, *thetab_, *phit_),
                                       dataset_B->sumEntries(), RooFit::ExpectedData(true));
        cr_hist_B_bar =
            pdf_B_bar.generateBinned(RooArgSet(*thetat_, *thetab_, *phit_),
                                     dataset_B_bar->sumEntries(), RooFit::ExpectedData(true));
        cr_hist->add(*cr_hist_B_bar);

        // Set the current directory back to the one for plots (ugly ROOT stuff)
        if (output_file_) {
            output_file_->cd();
        }
        tools::PlotPull2D(*thetat_, *thetab_, hist, *cr_hist);
        tools::PlotPull2D(*thetat_, *phit_, hist, *cr_hist);
        tools::PlotPull2D(*thetab_, *phit_, hist, *cr_hist);

        tools::PlotVars2D(*thetat_, *thetab_, *cr_hist);
        tools::PlotVars2D(*thetat_, *phit_, *cr_hist);
        tools::PlotVars2D(*thetab_, *phit_, *cr_hist);

        const double chi2 = Calculate3DChi2(hist, *cr_hist);
        Log::print(Log::info, "Chi2 = %f\n", chi2);

        // SaveChi2Scan(sim_pdf, ap_);
        // SaveChi2Scan(sim_pdf, apa_);
        // SaveChi2Scan(sim_pdf, a0_);
        // SaveChi2Scan(sim_pdf, ata_);
        // SaveChi2Scan(sim_pdf, ap_, margin_ap);
        // SaveChi2Scan(sim_pdf, apa_, margin_apa);
        // SaveChi2Scan(sim_pdf, a0_, margin_a0);
        // SaveChi2Scan(sim_pdf, ata_, margin_ata);

        delete dataset_B;
        delete dataset_B_bar;
    }
}

void FitterCPV::FitAngularCRSCF() {
    Log::print(Log::info, "Fitting %i events.\n", dataset_->numEntries());

    AngularPDF cr_pdf_B("cr_pdf_B", "cr_pdf_B", false, efficiency_model_, efficiency_files_, *thetat_, *thetab_, *phit_, *ap_,
                     *apa_, *a0_, *ata_);
    AngularPDF cr_pdf_B_bar("cr_pdf_B_bar", "cr_pdf_B_bar", true, efficiency_model_, efficiency_files_, *thetat_, *thetab_,
                         *phit_, *ap_, *apa_, *a0_, *ata_);

    RooRealVar cr_scf_f("cr_scf_f", "f_{cr}", constants::fraction_cr_of_crscf, 0.80, 0.99);
    cr_scf_f.setConstant();

    RooAddPdf pdf_B("pdf_B", "pdf_B", RooArgList(cr_pdf_B, *scf_angular_pdf_), RooArgList(cr_scf_f));
    RooAddPdf pdf_B_bar("pdf_B_bar", "pdf_B_bar", RooArgList(cr_pdf_B_bar, *scf_angular_pdf_), RooArgList(cr_scf_f));

    RooSimultaneous sim_pdf("sim_pdf", "sim_pdf", *decaytype_);
    sim_pdf.addPdf(pdf_B, "a");
    sim_pdf.addPdf(pdf_B_bar, "ab");
    sim_pdf.addPdf(pdf_B, "b");
    sim_pdf.addPdf(pdf_B_bar, "bb");

    result_ = sim_pdf.fitTo(*dataset_, RooFit::Minimizer("Minuit2"), RooFit::Hesse(false),
                            RooFit::Minos(false), RooFit::Save(true), RooFit::NumCPU(num_CPUs_));

    if (result_) {
        result_->Print();
    }

    if (make_plots_) {
        // PDF for B_bar differs, so we have to generate them separately. (One
        // could also use *extended* RooSimultaneous or 'ProtoData' for
        // RooSimultaneous.)
        RooDataHist* cr_hist = cr_pdf_B.generateBinned(RooArgSet(*thetat_, *thetab_, *phit_), 1000,
                                                    RooFit::ExpectedData(true));
        RooDataHist* cr_hist_B_bar = cr_pdf_B_bar.generateBinned(RooArgSet(*thetat_, *thetab_, *phit_),
                                                              1000, RooFit::ExpectedData(true));
        // Add histos from both particle and anti-particle PDFs to create the
        // final RooHistPdf.
        cr_hist->add(*cr_hist_B_bar);
        RooHistPdf cr_histpdf("cr_histpdf", "cr_histpdf", RooArgSet(*thetat_, *thetab_, *phit_),
                              *cr_hist);

        RooDataHist* scf_hist = scf_angular_pdf_->generateBinned(RooArgSet(*thetat_, *thetab_, *phit_), 1000,
                                                       RooFit::ExpectedData(true));
        RooHistPdf scf_histpdf("scf_histpdf", "scf_histpdf", RooArgSet(*thetat_, *thetab_, *phit_),
                               *scf_hist);

        RooRealVar cr_scf_f_gen("cr_scf_f_gen", "cr_scf_f_gen", (1 + cr_scf_f.getVal())/2);
        RooAddPdf all_histpdf("all_histpdf", "all_histpdf", RooArgList(cr_histpdf, scf_histpdf),
                              cr_scf_f);

        // Set the current directory back to the one for plots (ugly ROOT stuff)
        if (output_file_) {
            output_file_->cd();
        }

        // const double margin_ap = 0.005;
        // const double margin_apa = 0.25;
        // const double margin_a0 = 0.0015;
        // const double margin_ata = 0.3;
        // SaveLikelihoodScan(sim_pdf, ap_, margin_ap);
        // SaveLikelihoodScan(sim_pdf, apa_, margin_apa);
        // SaveLikelihoodScan(sim_pdf, a0_, margin_a0);
        // SaveLikelihoodScan(sim_pdf, ata_, margin_ata);
        // SaveLikelihoodScan(sim_pdf, ap_, apa_, margin_ap, margin_apa);
        // SaveLikelihoodScan(sim_pdf, ap_, a0_, margin_ap, margin_a0);
        // SaveLikelihoodScan(sim_pdf, ap_, ata_, margin_ap, margin_ata);
        // SaveLikelihoodScan(sim_pdf, apa_, a0_, margin_apa, margin_a0);
        // SaveLikelihoodScan(sim_pdf, apa_, ata_, margin_ap, margin_ata);
        // SaveLikelihoodScan(sim_pdf, a0_, ata_, margin_a0, margin_ata);

        // PlotWithPull(*thetat_, *dataset_, cr_histpdf);
        // PlotWithPull(*thetab_, *dataset_, cr_histpdf);
        // PlotWithPull(*phit_, *dataset_, cr_histpdf);

        // RooHistPdf all_histpdf("all_histpdf", "all_histpdf", RooArgSet(*thetat_, *thetab_,
        // *phit_), *all_hist);
        std::vector<RooAbsPdf*> components;
        components.push_back(&cr_histpdf);
        components.push_back(&scf_histpdf);
        // thetab_->setBins(100);
        PlotWithPull(*thetat_, *dataset_, all_histpdf, components);
        PlotWithPull(*thetab_, *dataset_, all_histpdf, components);
        PlotWithPull(*phit_, *dataset_, all_histpdf, components);

        // We need this to get the exact number of events to be generated for
        // the PDF distribution. This is unnecessary for the above, because we
        // create a RooHistPdf from it and that is normalized to the data when
        // plotted.
        RooDataSet* dataset_B = dynamic_cast<RooDataSet*>(
            dataset_->reduce("decaytype==decaytype::a||decaytype==decaytype::b"));
        RooDataSet* dataset_B_bar = dynamic_cast<RooDataSet*>(
            dataset_->reduce("decaytype==decaytype::ab||decaytype==decaytype::bb"));

        thetat_->setBins(40);
        thetab_->setBins(40);
        phit_->setBins(40);
        RooDataHist hist("hist", "hist", RooArgSet(*thetat_, *thetab_, *phit_), *dataset_);
        tools::PlotVars2D(*thetat_, *thetab_, hist);
        tools::PlotVars2D(*thetat_, *phit_, hist);
        tools::PlotVars2D(*thetab_, *phit_, hist);

        // We change the binning for 2D plots, so we have to generate new binned
        // dataset, to avoid dealing with rebinning.
        delete cr_hist;
        delete cr_hist_B_bar;
        cr_hist = pdf_B.generateBinned(RooArgSet(*thetat_, *thetab_, *phit_),
                                       dataset_B->sumEntries(), RooFit::ExpectedData(true));
        cr_hist_B_bar =
            pdf_B_bar.generateBinned(RooArgSet(*thetat_, *thetab_, *phit_),
                                     dataset_B_bar->sumEntries(), RooFit::ExpectedData(true));
        cr_hist->add(*cr_hist_B_bar);

        // Set the current directory back to the one for plots (ugly ROOT stuff)
        if (output_file_) {
            output_file_->cd();
        }
        tools::PlotPull2D(*thetat_, *thetab_, hist, *cr_hist);
        tools::PlotPull2D(*thetat_, *phit_, hist, *cr_hist);
        tools::PlotPull2D(*thetab_, *phit_, hist, *cr_hist);

        const double chi2 = Calculate3DChi2(hist, *cr_hist);
        std::cout << "Chi2 = " << chi2 << std::endl;

        // SaveChi2Scan(sim_pdf, ap_, margin_apa);
        // SaveChi2Scan(sim_pdf, apa_, margin_apa);
        // SaveChi2Scan(sim_pdf, a0_, margin_a0);
        // SaveChi2Scan(sim_pdf, ata_, margin_ata);

        delete dataset_B;
        delete dataset_B_bar;
    }
}

void FitterCPV::FitAngularAll() {
    Log::print(Log::info, "Fitting %i events.\n", dataset_->numEntries());

    AngularPDF cr_pdf_B("cr_pdf_B", "cr_pdf_B", false, efficiency_model_, efficiency_files_, *thetat_, *thetab_, *phit_, *ap_,
                     *apa_, *a0_, *ata_);
    AngularPDF cr_pdf_B_bar("cr_pdf_B_bar", "cr_pdf_B_bar", true, efficiency_model_, efficiency_files_, *thetat_, *thetab_,
                         *phit_, *ap_, *apa_, *a0_, *ata_);

    RooRealVar cr_f("cr_f", "f_{cr}", constants::fraction_cr_of_crscfbkg, 0.10, 0.99);
    RooRealVar scf_f("scf_f", "f_{scf}", constants::fraction_scf_of_crscfbkg, 0.10, 0.99);

    cr_f.setConstant();
    scf_f.setConstant();

    RooAddPdf pdf_B("pdf_B", "pdf_B", RooArgList(cr_pdf_B, *scf_angular_pdf_, *bkg_angular_pdf_), RooArgList(cr_f, scf_f));
    RooAddPdf pdf_B_bar("pdf_B_bar", "pdf_B_bar", RooArgList(cr_pdf_B_bar, *scf_angular_pdf_, *bkg_angular_pdf_), RooArgList(cr_f, scf_f));

    RooSimultaneous sim_pdf("sim_pdf", "sim_pdf", *decaytype_);
    sim_pdf.addPdf(pdf_B, "a");
    sim_pdf.addPdf(pdf_B_bar, "ab");
    sim_pdf.addPdf(pdf_B, "b");
    sim_pdf.addPdf(pdf_B_bar, "bb");

    result_ = sim_pdf.fitTo(*dataset_, RooFit::Minimizer("Minuit2"), RooFit::Hesse(false),
                            RooFit::Minos(false), RooFit::Save(true), RooFit::NumCPU(num_CPUs_));

    if (result_) {
        result_->Print();
    }

    if (make_plots_) {
        // PDF for B_bar differs, so we have to generate them separately. (One
        // could also use *extended* RooSimultaneous or 'ProtoData' for
        // RooSimultaneous.)
        RooDataHist* cr_hist = cr_pdf_B.generateBinned(RooArgSet(*thetat_, *thetab_, *phit_), 1000,
                                                    RooFit::ExpectedData(true));
        RooDataHist* cr_hist_B_bar = cr_pdf_B_bar.generateBinned(RooArgSet(*thetat_, *thetab_, *phit_),
                                                              1000, RooFit::ExpectedData(true));
        // Add histos from both particle and anti-particle PDFs to create the
        // final RooHistPdf.
        cr_hist->add(*cr_hist_B_bar);
        RooHistPdf cr_histpdf("cr_histpdf", "cr_histpdf", RooArgSet(*thetat_, *thetab_, *phit_),
                              *cr_hist);

        RooDataHist* scf_hist = scf_angular_pdf_->generateBinned(RooArgSet(*thetat_, *thetab_, *phit_), 1000,
                                                       RooFit::ExpectedData(true));
        RooHistPdf scf_histpdf("scf_histpdf", "scf_histpdf", RooArgSet(*thetat_, *thetab_, *phit_),
                               *scf_hist);

        RooDataHist* bkg_hist = bkg_angular_pdf_->generateBinned(RooArgSet(*thetat_, *thetab_, *phit_), 1000,
                                                       RooFit::ExpectedData(true));
        RooHistPdf bkg_histpdf("bkg_histpdf", "bkg_histpdf", RooArgSet(*thetat_, *thetab_, *phit_),
                               *bkg_hist);

        RooAddPdf all_histpdf("all_histpdf", "all_histpdf", RooArgList(cr_histpdf, scf_histpdf, bkg_histpdf),
                              RooArgList(cr_f, scf_f));

        // Set the current directory back to the one for plots (ugly ROOT stuff)
        if (output_file_) {
            output_file_->cd();
        }

        // const double margin_ap = 0.005;
        // const double margin_apa = 0.25;
        // const double margin_a0 = 0.0015;
        // const double margin_ata = 0.3;
        // SaveLikelihoodScan(sim_pdf, ap_, margin_ap);
        // SaveLikelihoodScan(sim_pdf, apa_, margin_apa);
        // SaveLikelihoodScan(sim_pdf, a0_, margin_a0);
        // SaveLikelihoodScan(sim_pdf, ata_, margin_ata);
        // SaveLikelihoodScan(sim_pdf, ap_, apa_, margin_ap, margin_apa);
        // SaveLikelihoodScan(sim_pdf, ap_, a0_, margin_ap, margin_a0);
        // SaveLikelihoodScan(sim_pdf, ap_, ata_, margin_ap, margin_ata);
        // SaveLikelihoodScan(sim_pdf, apa_, a0_, margin_apa, margin_a0);
        // SaveLikelihoodScan(sim_pdf, apa_, ata_, margin_ap, margin_ata);
        // SaveLikelihoodScan(sim_pdf, a0_, ata_, margin_a0, margin_ata);

        // PlotWithPull(*thetat_, *dataset_, cr_histpdf);
        // PlotWithPull(*thetab_, *dataset_, cr_histpdf);
        // PlotWithPull(*phit_, *dataset_, cr_histpdf);

        // RooHistPdf all_histpdf("all_histpdf", "all_histpdf", RooArgSet(*thetat_, *thetab_,
        // *phit_), *all_hist);
        std::vector<RooAbsPdf*> components;
        components.push_back(&cr_histpdf);
        components.push_back(&scf_histpdf);
        components.push_back(&bkg_histpdf);
        // thetab_->setBins(100);
        PlotWithPull(*thetat_, *dataset_, all_histpdf, components);
        PlotWithPull(*thetab_, *dataset_, all_histpdf, components);
        PlotWithPull(*phit_, *dataset_, all_histpdf, components);

        // We need this to get the exact number of events to be generated for
        // the PDF distribution. This is unnecessary for the above, because we
        // create a RooHistPdf from it and that is normalized to the data when
        // plotted.
        RooDataSet* dataset_B = dynamic_cast<RooDataSet*>(
            dataset_->reduce("decaytype==decaytype::a||decaytype==decaytype::b"));
        RooDataSet* dataset_B_bar = dynamic_cast<RooDataSet*>(
            dataset_->reduce("decaytype==decaytype::ab||decaytype==decaytype::bb"));

        thetat_->setBins(40);
        thetab_->setBins(40);
        phit_->setBins(40);
        RooDataHist hist("hist", "hist", RooArgSet(*thetat_, *thetab_, *phit_), *dataset_);
        tools::PlotVars2D(*thetat_, *thetab_, hist);
        tools::PlotVars2D(*thetat_, *phit_, hist);
        tools::PlotVars2D(*thetab_, *phit_, hist);

        // We change the binning for 2D plots, so we have to generate new binned
        // dataset, to avoid dealing with rebinning.
        delete cr_hist;
        delete cr_hist_B_bar;

        RooDataHist* all_hist = pdf_B.generateBinned(RooArgSet(*thetat_, *thetab_, *phit_),
                                       dataset_B->sumEntries(), RooFit::ExpectedData(true));
        RooDataHist* all_hist_B_bar =
            pdf_B_bar.generateBinned(RooArgSet(*thetat_, *thetab_, *phit_),
                                     dataset_B_bar->sumEntries(), RooFit::ExpectedData(true));
        all_hist->add(*all_hist_B_bar);

        // Set the current directory back to the one for plots (ugly ROOT stuff)
        if (output_file_) {
            output_file_->cd();
        }
        tools::PlotPull2D(*thetat_, *thetab_, hist, *all_hist);
        tools::PlotPull2D(*thetat_, *phit_, hist, *all_hist);
        tools::PlotPull2D(*thetab_, *phit_, hist, *all_hist);

        const double chi2 = Calculate3DChi2(hist, *all_hist);
        std::cout << "Chi2 = " << chi2 << std::endl;

        // SaveChi2Scan(sim_pdf, ap_, margin_apa);
        // SaveChi2Scan(sim_pdf, apa_, margin_apa);
        // SaveChi2Scan(sim_pdf, a0_, margin_a0);
        // SaveChi2Scan(sim_pdf, ata_, margin_ata);

        delete dataset_B;
        delete dataset_B_bar;
    }
}

RooDataSet* FitterCPV::ReduceDataToFitRange(const rapidjson::Document& config) {
    RooDataSet* reduced_dataset = 0;
    std::ostringstream reduce_string;
    bool first = true;

    for (auto var : dataset_vars_) {
        const char* var_name = (**var).GetName();
        if (config["fitRanges"].HasMember(var_name)) {
            if (first == false) {
                reduce_string << " && ";
            } else {
                first = false;
            }

            reduce_string << "(";
            const rapidjson::SizeType num_ranges = config["fitRanges"][var_name].Size();
            for (rapidjson::SizeType i = 0; i < num_ranges; i++) {
                double low, high;
                low = config["fitRanges"][var_name][i][0].GetDouble();
                high = config["fitRanges"][var_name][i][1].GetDouble();

                reduce_string << "(" << var_name << " > " << low << " && " << var_name << " < "
                              << high << ")";
                if (i < num_ranges - 1) {
                    reduce_string << " || ";
                }
            }
            reduce_string << ")";
        }
    }

    reduced_dataset = dynamic_cast<RooDataSet*>(dataset_->reduce(reduce_string.str().c_str()));
    Log::print(Log::info, "Dataset after reduce: %i\n", reduced_dataset->numEntries());
    Log::print(Log::debug, "Reduce string: %s\n", reduce_string.str().c_str());
    return reduced_dataset;
}

/* static */ rapidjson::Document FitterCPV::ReadJSONConfig(const char* filename) {
    std::ifstream filestream(filename);
    std::stringstream buffer;
    buffer << filestream.rdbuf();

    rapidjson::Document config;
    config.Parse(buffer.str().c_str());
    return config;
}

std::string FitterCPV::ApplyJSONConfig(const rapidjson::Document& config) {
    // Store the applied config in the resulting ROOT file for reference
    rapidjson::StringBuffer buffer;
    rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
    config.Accept(writer);
    const char* json_string = buffer.GetString();

    if (config.HasMember("fitRanges")) {
        ChangeFitRanges(config["fitRanges"]);
        // dataset_ = ReduceDataToFitRange(config);
    }

    if (config.HasMember("modelParameters")) {
        ChangeModelParameters(config["modelParameters"]);
    }

    std::string json_config_string = json_string;
    return json_config_string;
}

void FitterCPV::ChangeFitRanges(const rapidjson::GenericValue<rapidjson::UTF8<char>>& config) {
    for (auto var : dataset_vars_) {
        const char* var_name = (**var).GetName();
        if (config.HasMember(var_name)) {
            if (config[var_name].HasMember("min")) {
                (**var).setMin(config[var_name]["min"].GetDouble());
            }
            if (config[var_name].HasMember("max")) {
                (**var).setMax(config[var_name]["max"].GetDouble());
            }
        }
    }
}

void FitterCPV::ChangeModelParameters(const rapidjson::GenericValue<rapidjson::UTF8<char>>& config) {
    for (rapidjson::Value::ConstMemberIterator itr = config.MemberBegin();
         itr != config.MemberEnd(); ++itr) {
        assert(scf_parameters_argset_.find(itr->name.GetString()) ||
               bkg_parameters_argset_.find(itr->name.GetString()));
        scf_parameters_argset_.setRealValue(itr->name.GetString(), itr->value.GetDouble());
        bkg_parameters_argset_.setRealValue(itr->name.GetString(), itr->value.GetDouble());
    }
}

const void FitterCPV::LogCLIArguments(int argc, char* argv[]) {
    // Set the current directory back to the one for plots (ugly ROOT stuff)
    if (output_file_) {
        output_file_->cd();
    }

    std::string str;
    for (int i = 0; i < argc; i++) {
        str += argv[i];
        str += " ";
    }
    // Remove the final space
    str.pop_back();
    TNamed cli_arguments("cli_arguments", str.c_str());
    cli_arguments.Write();
}

const void FitterCPV::LogEnvironmentMetadata() {
    // Set the current directory back to the one for plots (ugly ROOT stuff)
    if (output_file_) {
        output_file_->cd();
    }

    TNamed root_version("root_version", ROOT_RELEASE);
    root_version.Write();

    char buffer[100];
    gethostname(buffer, 100);
    TNamed hostname("hostname", buffer);
    hostname.Write();

    const time_t now = time(0);
    const char* local_time_string = ctime(&now);
    TNamed local_date("local_date", local_time_string);

    tm* gmtm = gmtime(&now);
    const char* utc_time_string = asctime(gmtm);
    TNamed utc_date("utc_date", utc_time_string);

    TNamed git_version("git_version", gitversion);
    git_version.Write();
}

void FitterCPV::GenerateToys(const int num_events, const int num_toys) {
    Log::print(Log::info, "Generating dataset with %i events...\n", num_events);

    // vars with phiweak = phiweak + pi, r = 0.01
    double par_input[] = {
        0.269,  // ap
        0.56,   // apa
        0.941,  // a0
        3.11,   // ata

        0.00816649,   // xp
        0.00532961,   // x0
        0.00829682,   // xt
        -0.00659988,  // yp
        -0.0084614,   // y0
        +0.00373802,  // yt

        -0.0102141,   // xpb
        -0.00845715,  // x0b
        -0.00587413,  // xtb
        +0.00243354,  // ypb
        +0.00533635,  // y0b
        -0.00695015   // ytb
    };

    //  // vars with phiweak = phiweak + pi, r = 0.10
    //  double par_input[] = {
    //          0.269,
    //          0.56,
    //          0.941,
    //          3.11,
    //
    //          0.0816649, // xp
    //          0.0532961, // x0
    //          0.0829682, // xt
    //          -0.0659988, // yp
    //          -0.084614,  // y0
    //          +0.0373802, // yt
    //
    //          -0.102141,  // xpb
    //          -0.0845715, // x0b
    //          -0.0587413, // xtb
    //          +0.0243354, // ypb
    //          +0.0533635, // y0b
    //          -0.0695015  // ytb
    //  };

    RooRealVar ap("ap", "ap", par_input[0], 0, 0.5);
    RooRealVar apa("apa", "apa", par_input[1], 0, 1);
    RooRealVar a0("a0", "a0", par_input[2], 0.8, 1);
    RooRealVar a0a("a0a", "a0a", 0);
    RooFormulaVar at("at", "sqrt(1-ap*ap-a0*a0)", RooArgSet(ap, a0));
    RooRealVar ata("ata", "ata", par_input[3], 2, 4);

    RooRealVar xp("xp", "xp", par_input[4], -0.2, 0.2);
    RooRealVar x0("x0", "x0", par_input[5], -0.2, 0.2);
    RooRealVar xt("xt", "xt", par_input[6], -0.2, 0.2);

    RooRealVar yp("yp", "yp", par_input[7], -0.2, 0.2);
    RooRealVar y0("y0", "y0", par_input[8], -0.2, 0.2);
    RooRealVar yt("yt", "yt", par_input[9], -0.2, 0.2);

    RooRealVar xpb("xpb", "xpb", par_input[10], -0.2, 0.2);
    RooRealVar x0b("x0b", "x0b", par_input[11], -0.2, 0.2);
    RooRealVar xtb("xtb", "xtb", par_input[12], -0.2, 0.2);

    RooRealVar ypb("ypb", "ypb", par_input[13], -0.2, 0.2);
    RooRealVar y0b("y0b", "y0b", par_input[14], -0.2, 0.2);
    RooRealVar ytb("ytb", "ytb", par_input[15], -0.2, 0.2);

    DtCPPDF mixing_pdf_a("mixing_pdf_a", "mixing_pdf_a", false, true, perfect_tagging_,
                         efficiency_model_, efficiency_files_, *thetat_, *thetab_, *phit_, ap, apa, a0, ata, xp, x0,
                         xt, yp, y0, yt, *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_,
                         *benergy_, *mbc_, *vrntrk_, *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_,
                         *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF mixing_pdf_ab("mixing_pdf_ab", "mixing_pdf_ab", true, true, perfect_tagging_,
                          efficiency_model_, efficiency_files_, *thetat_, *thetab_, *phit_, ap, apa, a0, ata, xpb, x0b,
                          xtb, ypb, y0b, ytb, *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_,
                          *shcosthb_, *benergy_, *mbc_, *vrntrk_, *vrzerr_, *vrchi2_, *vrndf_,
                          *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF mixing_pdf_b("mixing_pdf_b", "mixing_pdf_b", false, false, perfect_tagging_,
                         efficiency_model_, efficiency_files_, *thetat_, *thetab_, *phit_, ap, apa, a0, ata, xp, x0,
                         xt, yp, y0, yt, *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_,
                         *benergy_, *mbc_, *vrntrk_, *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_,
                         *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF mixing_pdf_bb("mixing_pdf_bb", "mixing_pdf_bb", true, false, perfect_tagging_,
                          efficiency_model_, efficiency_files_, *thetat_, *thetab_, *phit_, ap, apa, a0, ata, xpb, x0b,
                          xtb, ypb, y0b, ytb, *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_,
                          *shcosthb_, *benergy_, *mbc_, *vrntrk_, *vrzerr_, *vrchi2_, *vrndf_,
                          *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    dt_->setRange(-15, 15);

    RooArgSet varsToGenerate(*dt_, *thetat_, *thetab_, *phit_);
    RooArgSet varsToAdd(*tagwtag_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
                        *vrzerr_, *vrchi2_);
    varsToAdd.add(*vrndf_);
    varsToAdd.add(*vtntrk_);
    varsToAdd.add(*vtzerr_);
    varsToAdd.add(*vtchi2_);
    varsToAdd.add(*vtndf_);
    varsToAdd.add(*vtistagl_);
    varsToAdd.add(*evmcflag_);
    varsToAdd.add(*tagqr_);
    varsToAdd.add(*vrusable_);
    varsToAdd.add(*vtusable_);

    RooRandom::randomGenerator()->SetSeed(0);

    /// For some reason this is the generated fraction of favored decays, have to check why
    const double frac_fav = 0.812;
    const double frac_sup = 1 - frac_fav;

    const int num_fav = num_events * frac_fav / 2.0;
    const int num_sup = num_events * frac_sup / 2.0;

    Log::print(Log::info, "Creating GenSpecs...\n");
    RooAbsPdf::GenSpec* genSpec_a =
        mixing_pdf_a.prepareMultiGen(varsToGenerate, RooFit::NumEvents(num_fav));
    Log::print(Log::info, "1/4 GenSpecs ready\n");
    RooAbsPdf::GenSpec* genSpec_ab =
        mixing_pdf_ab.prepareMultiGen(varsToGenerate, RooFit::NumEvents(num_fav));
    Log::print(Log::info, "2/4 GenSpecs ready\n");
    RooAbsPdf::GenSpec* genSpec_b =
        mixing_pdf_b.prepareMultiGen(varsToGenerate, RooFit::NumEvents(num_sup));
    Log::print(Log::info, "3/4 GenSpecs ready\n");
    RooAbsPdf::GenSpec* genSpec_bb =
        mixing_pdf_bb.prepareMultiGen(varsToGenerate, RooFit::NumEvents(num_sup));
    Log::print(Log::info, "4/4 GenSpecs ready\n");

    RooAbsData::setDefaultStorageType(RooAbsData::Tree);
    RooDataSet* dataset;
    TString filename;

    RooRealVar btagmcli("btagmcli", "btagmcli", -511, 511);
    RooRealVar brecflav("brecflav", "brecflav", -1, 1);

    for (int i = 1; i <= num_toys; i++) {
        //  RooDataSet* dataset_a = mixing_pdf_a.generate(varsToGenerate,
        //  RooFit::NumEvents(num_fav));
        RooDataSet* dataset_a = mixing_pdf_a.generate(*genSpec_a);
        decaytype_->setLabel("a");
        dataset_a->addColumn(*decaytype_);
        brecflav.setVal(1);
        btagmcli.setVal(-511);
        dataset_a->addColumn(brecflav);
        dataset_a->addColumn(btagmcli);
        dataset = (RooDataSet*)dataset_a->Clone();
        delete dataset_a;
        Log::print(Log::info, "1/4 decay type ready.\n");

        //  RooDataSet* dataset_ab = mixing_pdf_ab.generate(varsToGenerate,
        //  RooFit::NumEvents(num_fav));
        RooDataSet* dataset_ab = mixing_pdf_ab.generate(*genSpec_ab);
        decaytype_->setLabel("ab");
        dataset_ab->addColumn(*decaytype_);
        brecflav.setVal(-1);
        btagmcli.setVal(511);
        dataset_ab->addColumn(brecflav);
        dataset_ab->addColumn(btagmcli);
        dataset->append(*dataset_ab);
        delete dataset_ab;
        Log::print(Log::info, "2/4 decay types ready.\n");

        //  RooDataSet* dataset_b = mixing_pdf_b.generate(varsToGenerate,
        //  RooFit::NumEvents(num_sup));
        RooDataSet* dataset_b = mixing_pdf_b.generate(*genSpec_b);
        decaytype_->setLabel("b");
        dataset_b->addColumn(*decaytype_);
        brecflav.setVal(1);
        btagmcli.setVal(511);
        dataset_b->addColumn(brecflav);
        dataset_b->addColumn(btagmcli);
        dataset->append(*dataset_b);
        delete dataset_b;
        Log::print(Log::info, "3/4 decay types ready.\n");

        //  RooDataSet* dataset_bb = mixing_pdf_bb.generate(varsToGenerate,
        //  RooFit::NumEvents(num_sup));
        RooDataSet* dataset_bb = mixing_pdf_bb.generate(*genSpec_bb);
        decaytype_->setLabel("bb");
        dataset_bb->addColumn(*decaytype_);
        brecflav.setVal(-1);
        btagmcli.setVal(-511);
        dataset_bb->addColumn(brecflav);
        dataset_bb->addColumn(btagmcli);
        dataset->append(*dataset_bb);
        delete dataset_bb;
        Log::print(Log::info, "4/4 decay types ready.\n");

        dataset->addColumns(varsToAdd);

        vrzerr_formula_ =
            new RooFormulaVar("vrzerr", "#sigma z_{rec} [cm]", "sqrt(vrerr6)", RooArgSet(*vrerr6_));
        vtzerr_formula_ =
            new RooFormulaVar("vtzerr", "#sigma z_{tag} [cm]", "sqrt(vterr6)", RooArgSet(*vterr6_));

        RooFormulaVar vrvtxz("vrvtxz", "dt*(0.425*0.0299792458/2)", *dt_);
        RooFormulaVar vtvtxz("vtvtxz", "-dt*(0.425*0.0299792458/2)", *dt_);
        RooFormulaVar vrerr6("vrerr6", "vrzerr^2", *vrzerr_);
        RooFormulaVar vterr6("vterr6", "vtzerr^2", *vtzerr_);

        dataset->addColumns(RooArgList(vrvtxz, vtvtxz, vrerr6, vterr6));

        const TTree& tree = ((RooTreeDataStore*)dataset->store())->tree();
        TTree* new_tree = static_cast<TTree*>(tree.Clone("h2000"));

        filename = "test_";
        filename += i;
        filename += ".root";
        TFile f(filename, "RECREATE");
        new_tree->Write();
        f.Close();

        delete dataset;
        // dataset_ = dataset;
    }

    //  return dataset;

    //  result_ = mixing_pdf_a.fitTo(*new_dataset,
    // RooFit::ConditionalObservables(conditional_vars_argset_),
    //          RooFit::Minimizer("Minuit2"), RooFit::Hesse(false), RooFit::Minos(false),
    //          RooFit::Save(true), RooFit::NumCPU(num_CPUs_));
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
void FitterCPV::PlotWithPull(const RooRealVar& var, const RooAbsData& data, const RooAbsPdf& pdf,
                             const std::vector<RooAbsPdf*> components, const char* title) const {
    TString name = pdf.GetName();
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

    // Renormalization required for certain plots, when using multiple CPUs but
    // not with a single CPU; possible RooFit bug
    double norm = 1;
    if (!name.Contains("hist") && num_CPUs_ > 1) {
        norm = 1.0 / data.numEntries();
    }

    // Plot components before the total PDF as the pull plots are made from the
    // last plotted PDF
    const int colors[] = {4, 7, 5};
    int i = 0;
    for (auto component : components) {
        pdf.plotOn(plot, RooFit::ProjWData(conditional_vars_argset_, data, kFALSE),
                   RooFit::NumCPU(num_CPUs_), RooFit::LineColor(colors[i++]),
                   RooFit::Components(*component), RooFit::Normalization(norm));
    }

    pdf.plotOn(plot, RooFit::ProjWData(conditional_vars_argset_, data, kFALSE),
               RooFit::NumCPU(num_CPUs_), RooFit::Normalization(norm));

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
    RooPlot* plot_pull_ = var.frame(RooFit::Title("Pull Distribution"));
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
TPaveText* FitterCPV::CreateStatBox(const double chi2, const int ndof, const bool position_top,
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
void FitterCPV::ReadInFile(std::vector<const char*> file_names, const int& num_events) {
    TChain* input_chain = new TChain("h2000");
    for (auto file_name : file_names) {
        input_chain->Add(file_name);
    }
    
    Log::print(Log::info, "Reading %i input files...\n", file_names.size());
    
    TTree* input_tree;
    if (num_events) {
	if (file_names.size() > 1) {
            Log::print(Log::warning,
                       "You limited the number of events to read, while reading multiple "
		       "files. Since this limiting works sequentially not randomly, you "
		       "will probably not get the result you want!\n");
	}
        input_tree = input_chain->CloneTree(num_events);
    } else {
        input_tree = input_chain->CloneTree();
    }

    TString common_cuts = tools::GetCommonCutsString();

    TString a_cuts;
    TString ab_cuts;
    TString b_cuts;
    TString bb_cuts;
    if (perfect_tagging_) {
        a_cuts = "brecflav==1&&btagmcli<0";
        ab_cuts = "brecflav==-1&&btagmcli>0";
        b_cuts = "brecflav==1&&btagmcli>0";
        bb_cuts = "brecflav==-1&&btagmcli<0";
    } else {
        if (generator_level_) {
            Log::print(Log::warning,
                       "Attempting to use realistic tagging with generator level fit. This will "
                       "probably end badly. Consider using the '--perfect-tag' switch.\n");
        }
        a_cuts = "brecflav==1&&tagqr<0";
        ab_cuts = "brecflav==-1&&tagqr>0";
        b_cuts = "brecflav==1&&tagqr>0";
        bb_cuts = "brecflav==-1&&tagqr<0";
    }

    // A temporary RooDataSet is created from the whole tree and then we apply cuts to get
    // the 4 different B and f flavor datasets, as that is faster then reading the tree 4 times
    RooDataSet* temp_dataset =
        new RooDataSet("dataset", "dataset", input_tree, dataset_vars_argset_, common_cuts);
    Log::print(Log::debug, "Num events passing common cuts: %i\n", temp_dataset->numEntries());

    //  separate conditional vars from stuff like thetat

    // We add an identifying label to each of the 4 categories and then combine it into a single
    // dataset for RooSimultaneous fitting
    RooDataSet* dataset_a = static_cast<RooDataSet*>(temp_dataset->reduce(a_cuts));
    decaytype_->setLabel("a");
    dataset_a->addColumn(*decaytype_);

    RooDataSet* dataset_ab = static_cast<RooDataSet*>(temp_dataset->reduce(ab_cuts));
    decaytype_->setLabel("ab");
    dataset_ab->addColumn(*decaytype_);

    RooDataSet* dataset_b = static_cast<RooDataSet*>(temp_dataset->reduce(b_cuts));
    decaytype_->setLabel("b");
    dataset_b->addColumn(*decaytype_);

    RooDataSet* dataset_bb = static_cast<RooDataSet*>(temp_dataset->reduce(bb_cuts));
    decaytype_->setLabel("bb");
    dataset_bb->addColumn(*decaytype_);

    delete temp_dataset;

    // thetat_->setBins(10);
    // thetab_->setBins(10);
    // phit_->setBins(10);
    // RooDataHist datahist_a("datahist_a", "datahist_a", RooArgSet(*thetat_, *thetab_, *phit_), *dataset_a);
    // RooDataHist datahist_ab("datahist_ab", "datahist_ab", RooArgSet(*thetat_, *thetab_, *phit_), *dataset_ab);
    // RooDataHist datahist_b("datahist_b", "datahist_b", RooArgSet(*thetat_, *thetab_, *phit_), *dataset_b);
    // RooDataHist datahist_bb("datahist_bb", "datahist_bb", RooArgSet(*thetat_, *thetab_, *phit_), *dataset_bb);
    // tools::PlotPull2D(*thetat_, *thetab_, datahist_a, datahist_ab);
    // tools::PlotPull2D(*thetat_, *phit_, datahist_a, datahist_ab);
    // tools::PlotPull2D(*thetab_, *phit_, datahist_a, datahist_ab);
    // tools::PlotPull2D(*thetat_, *thetab_, datahist_b, datahist_bb);
    // tools::PlotPull2D(*thetat_, *phit_, datahist_b, datahist_bb);
    // tools::PlotPull2D(*thetab_, *phit_, datahist_b, datahist_bb);

    // dataset_ = static_cast<RooDataSet*>(dataset_bb->Clone());
    dataset_ = static_cast<RooDataSet*>(dataset_a->Clone());
    delete dataset_a;
    dataset_->append(*dataset_ab);
    delete dataset_ab;
    dataset_->append(*dataset_b);
    delete dataset_b;
    dataset_->append(*dataset_bb);
    delete dataset_bb;

    Log::print(Log::info, "Num events passing all initial cuts: %i\n", dataset_->numEntries());

    delete input_tree;

    vrzerr_ = static_cast<RooRealVar*>(dataset_->addColumn(*vrzerr_formula_));
    vrzerr_->setRange(0, 10000);
    vtzerr_ = static_cast<RooRealVar*>(dataset_->addColumn(*vtzerr_formula_));
    vtzerr_->setRange(0, 10000);

    conditional_vars_argset_.add(*vrzerr_);
    conditional_vars_argset_.add(*vtzerr_);
    conditional_vars_argset_.add(*decaytype_);

    // Bind the variables to the dataset, so that dataset->get(i) changes values of, e.g., expno_
    const RooArgSet* vars = dataset_->get();
    for (RooRealVar** var : dataset_vars_) {
        *var = static_cast<RooRealVar*>(vars->find((*var)->GetName()));
    }

    dataset_->get(1)->Print("v");
}

/**
 * Set the directory to which to ouput plots
 */
void FitterCPV::SetPlotDir(const char* plot_dir) {
    make_plots_ = true;
    tools::SetPlotDir(plot_dir);
}

/**
 * Fix requested parameters.
 *
 * Pars can be composed of a comma separated list of parameter names and/or the
 * following short-hand keywords which are expanded: all, xy, trans, nota0.
 */
bool FitterCPV::FixParameters(const char* pars) {
    std::string input = pars;

    // These are helper short-hand options to save some typing
    input = std::regex_replace(input, std::regex("all"),
                               "ap,apa,a0,ata,xp,x0,xt,yp,y0,yt,xpb,x0b,xtb,ypb,y0b,ytb");
    input = std::regex_replace(input, std::regex("xy"),
                               "xp,x0,xt,yp,y0,yt,xpb,x0b,xtb,ypb,y0b,ytb");
    input = std::regex_replace(input, std::regex("trans"), "ap,apa,a0,ata");
    input = std::regex_replace(input, std::regex("nota0"),
                               "ap,apa,ata,xp,x0,xt,yp,y0,yt,xpb,x0b,xtb,ypb,y0b,ytb");

    std::istringstream ss(input);
    std::string token;
    std::vector<std::string> par_vector;

    while (std::getline(ss, token, ',')) {
        par_vector.push_back(token);
    }

    RooRealVar* rooPar = 0;
    for (auto par : par_vector) {
        bool par_not_found = true;
        TIterator* parIter = parameters_argset_.createIterator();
        while ((rooPar = (RooRealVar*)parIter->Next())) {
            if (rooPar->GetName() == par) {
                rooPar->setConstant();
                par_not_found = false;
            }
        }
        if (par_not_found) {
            Log::print(Log::error, "Parameter %s doesn't exist.\n", par.c_str());
            return true;
        }
    }
    return false;
}

/**
 * Set file in which output will be saved.
 */
void FitterCPV::SetOutputFile(const char* filename) {
    output_file_ = new TFile(filename, "RECREATE");
}

/**
 * Create a string that holds initial values, fit results and errors.
 */
const std::string FitterCPV::CreateResultsString() {
    int numParameters = 54;
    double* parameters = new double[numParameters];

    parameters[0] = par_input_[0];
    parameters[1] = ap_->getVal();
    parameters[2] = ap_->getError();
    parameters[3] = par_input_[1];
    parameters[4] = apa_->getVal();
    parameters[5] = apa_->getError();
    parameters[6] = par_input_[2];
    parameters[7] = a0_->getVal();
    parameters[8] = a0_->getError();
    parameters[9] = 0;
    parameters[10] = a0a_->getVal();
    parameters[11] = a0a_->getError();
    parameters[12] = sqrt(1 - par_input_[0] * par_input_[0] - par_input_[2] * par_input_[2]);
    parameters[13] = at_->getVal();
    parameters[14] = at_->getPropagatedError(*result_);
    parameters[15] = par_input_[3];
    parameters[16] = ata_->getVal();
    parameters[17] = ata_->getError();

    if (!do_time_independent_fit_) {
        parameters[18] = par_input_[4];
        parameters[19] = xp_->getVal();
        parameters[20] = xp_->getError();
        parameters[21] = par_input_[5];
        parameters[22] = x0_->getVal();
        parameters[23] = x0_->getError();
        parameters[24] = par_input_[6];
        parameters[25] = xt_->getVal();
        parameters[26] = xt_->getError();
        parameters[27] = par_input_[7];
        parameters[28] = yp_->getVal();
        parameters[29] = yp_->getError();
        parameters[30] = par_input_[8];
        parameters[31] = y0_->getVal();
        parameters[32] = y0_->getError();
        parameters[33] = par_input_[9];
        parameters[34] = yt_->getVal();
        parameters[35] = yt_->getError();

        parameters[36] = par_input_[10];
        parameters[37] = xpb_->getVal();
        parameters[38] = xpb_->getError();
        parameters[39] = par_input_[11];
        parameters[40] = x0b_->getVal();
        parameters[41] = x0b_->getError();
        parameters[42] = par_input_[12];
        parameters[43] = xtb_->getVal();
        parameters[44] = xtb_->getError();
        parameters[45] = par_input_[13];
        parameters[46] = ypb_->getVal();
        parameters[47] = ypb_->getError();
        parameters[48] = par_input_[14];
        parameters[49] = y0b_->getVal();
        parameters[50] = y0b_->getError();
        parameters[51] = par_input_[15];
        parameters[52] = ytb_->getVal();
        parameters[53] = ytb_->getError();
    } else {
        numParameters = 18;
    }

    Int_t separators[] = {0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2,
                          0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2,
                          0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2};

    char buffer[100];
    std::string results;
    for (Int_t i = 0; i < numParameters; i++) {
        /// These parameters can be both positive and negative therefore the
        /// added + in front of positive numbers keeps the columns aligned
        if (i >= 18 && (i - 20) % 3 != 0) {
            snprintf(buffer, 100, "%+.5f ", parameters[i]);
        } else {
            snprintf(buffer, 100, "%.5f ", parameters[i]);
        }
        results += buffer;

        if (i == numParameters - 1) continue;
        if (separators[i] == 1)
            results += "| ";
        else if (separators[i] == 2)
            results += "|| ";
    }
    delete[] parameters;

    // Replace the final space with a newline
    results.pop_back();
    results.append("\n");
    return results;
}

/**
 * Create a Markdown table with fit results and pulls and return it in a string.
 */
const std::string FitterCPV::CreatePullTableString(const bool asymmetric) {
    std::stringstream ss;
    if (asymmetric) {
        ss << "Var | True    | Fit     | Err Lo  | Err Hi | Pull\n";
        ss << ":---|:-------:|:-------:|:-------:|:------:|-----:\n";
    } else {
        ss << "Var | True    | Fit     | Err    | Pull\n";
        ss << ":---|:-------:|:-------:|:------:|-----:\n";
    }
    bool at_not_processed = true;
    for (unsigned int i = 0; i < parameters_.size(); i++) {
        std::string name = (*parameters_[i])->GetName();

        // Skip 'x' and 'y' pars if we are doing a time-independent fit
        if (do_time_independent_fit_) {
            if (name.front() == 'x' || name.front() == 'y') continue;
        }

        double tru = par_input_[i];
        double fit = (*parameters_[i])->getVal();
        double err = (*parameters_[i])->getError();
        double err_low = (*parameters_[i])->getAsymErrorLo();
        double err_high = (*parameters_[i])->getAsymErrorHi();

        // Parameter 'at' is completely determined by 'ap' and 'a0', so it is
        // treated differently. It's not a member of parameters_ so it must be
        // added by hand.
        if (name == "ata" && at_not_processed) {
            name = "at";
            tru = sqrt(1 - par_input_[0] * par_input_[0] - par_input_[2] * par_input_[2]);
            fit = at_->getVal();
            err = at_->getPropagatedError(*result_);
            err_low = -err;
            err_high = err;
            at_not_processed = false;
            i--;
        }

        ss << std::setw(3) << std::left << name << " | ";
        ss << std::fixed;
        ss << std::setprecision(4);
        ss << std::showpos;
        ss << tru << " | ";
        ss << fit << " | ";
        if (asymmetric) {
            ss << std::noshowpos << err_low << " | ";
            ss << std::noshowpos << err_high << " | ";
        } else {
            ss << std::noshowpos << err << " | ";
        }
        ss << std::setprecision(2) << std::showpos;
        if (asymmetric) {
            if (fit < tru) {
                ss << (fit - tru) / err_high;
            } else {
                ss << (fit - tru) / fabs(err_low);
            }
        } else {
            ss << (fit - tru) / err;
        }
        ss << std::endl;
    }
    return ss.str();
}

/**
 * Create a LaTeX table with fit results and pulls and return it in a string.
 */
const std::string FitterCPV::CreateLatexPullTableString(const bool asymmetric) {
    std::stringstream ss;
    ss << "\\begin{table}\n";
    ss << "\t\\small\n";
    if (asymmetric) {
        ss << "\t\\begin{tabular}[]{lccccr}\n";
    } else {
        ss << "\t\\begin{tabular}[]{lcccr}\n";
    }
    ss << "\t\t\\toprule\n";
    if (asymmetric) {
        ss << "\t\tVar & True & Fit & Err. Low & Err. High & Pull \\\\\n";
    } else {
        ss << "\t\tVar & True & Fit & Err. & Pull \\\\\n";
    }
    ss << "\t\t\\midrule\n";
    bool at_not_processed = true;
    for (unsigned int i = 0; i < parameters_.size(); i++) {
        std::string name = "$";
        name += (*parameters_[i])->GetTitle();
        name += "$";
        std::replace(name.begin(), name.end(), '#', '\\');

        // Skip 'x' and 'y' pars if we are doing a time-independent fit
        if (do_time_independent_fit_) {
            if (name.find('x') != std::string::npos || name.find('y') != std::string::npos)
                continue;
        }

        double tru = par_input_[i];
        double fit = (*parameters_[i])->getVal();
        double err = (*parameters_[i])->getError();
        double err_low = (*parameters_[i])->getAsymErrorLo();
        double err_high = (*parameters_[i])->getAsymErrorHi();

        // Parameter 'at' is completely determined by 'ap' and 'a0', so it is
        // treated differently. It's not a member of parameters_ so it must be
        // added by hand.
        if (name == "$\\arg(a_{\\perp})$" && at_not_processed) {
            name = "$|a_{\\perp}|$";
            tru = sqrt(1 - par_input_[0] * par_input_[0] - par_input_[2] * par_input_[2]);
            fit = at_->getVal();
            err = at_->getPropagatedError(*result_);
            err_low = -err;
            err_high = err;
            at_not_processed = false;
            i--;
        }

        ss << "\t\t";
        ss << std::setw(21) << std::left << name << " & ";
        ss << std::fixed;
        ss << std::setprecision(4);
        ss << std::showpos;
        ss << tru << " & ";
        ss << fit << " & ";
        if (asymmetric) {
            ss << std::noshowpos << err_low << " & ";
            ss << std::noshowpos << err_high << " & ";
        } else {
            ss << std::noshowpos << err << " & ";
        }
        ss << std::setprecision(2) << std::showpos;
        if (asymmetric) {
            if (fit < tru) {
                ss << (fit - tru) / err_high;
            } else {
                ss << (fit - tru) / fabs(err_low);
            }
        } else {
            ss << (fit - tru) / err;
        }
        ss << " \\\\" << std::endl;
    }
    ss << "\t\t\\bottomrule\n";
    ss << "\t\\end{tabular}\n";
    ss << "\\end{table}\n";
    return ss.str();
}

/**
 * Save results, status codes, covariance matrix quality to the ROOT output file
 */
const void FitterCPV::LogResults() {
    // Set the current directory back to the one for plots (ugly ROOT stuff)
    if (output_file_) {
        output_file_->cd();
    }
    result_->Write();

    std::string results_string = CreateResultsString();
    TNamed txt_result("txt_result", results_string);
    txt_result.Write();

    char buffer[100];
    snprintf(buffer, 100, "%i", result_->covQual());
    TNamed cov_qual("cov_qual", buffer);
    cov_qual.Write();

    for (uint i = 0; i < result_->numStatusHistory(); i++) {
        TNamed status(result_->statusLabelHistory(i),
                      std::to_string(result_->statusCodeHistory(i)));
        status.Write();
    }
}

/**
 * Save results and initial values into a plain file.
 */
const void FitterCPV::SaveTXTResults(const char* filename) {
    std::string results_string = CreateResultsString();
    std::ofstream file(filename);
    file << results_string;
    file.close();
}

void FitterCPV::TestEfficiency() {
    Efficiency* eff = new Efficiency();
    TF3* f =
        new TF3("f", eff, &Efficiency::EfficiencyInterface, thetat_->getMin(), thetat_->getMax(),
                thetab_->getMin(), thetab_->getMax(), phit_->getMin(), phit_->getMax(), 0);

    TCanvas* c1 = new TCanvas("c1", "canvas", 500, 500);
    RooAbsReal* rf = RooFit::bindFunction(f, *thetat_, *thetab_, *phit_);
    RooPlot* plot = thetab_->frame();
    rf->plotOn(plot);
    plot->Draw();

    c1->SaveAs("eff_test.pdf");
}

TH3D* FitterCPV::GetBinnedEfficiency(std::vector<const char*> files, const int model) {
    Efficiency eff;
    for (auto file : files) {
        eff.ReadInFile(file);
    }
    TH3D* histo =
        new TH3D("histo", "histo", 100, thetat_->getMin(), thetat_->getMax(), 100,
                 thetab_->getMin(), thetab_->getMax(), 100, phit_->getMin(), phit_->getMax());

    double thetat, thetab, phit;
    for (int x = 1; x <= histo->GetXaxis()->GetNbins(); x++) {
        for (int y = 1; y <= histo->GetYaxis()->GetNbins(); y++) {
            for (int z = 1; z <= histo->GetZaxis()->GetNbins(); z++) {
                thetat = histo->GetXaxis()->GetBinCenter(x);
                thetab = histo->GetYaxis()->GetBinCenter(y);
                phit = histo->GetZaxis()->GetBinCenter(z);
                histo->SetBinContent(histo->GetBin(x, y, z),
                                     eff.GetEfficiency(thetat, thetab, phit, model));
            }
        }
    }

    return histo;
}

void FitterCPV::PlotEfficiency() {
    TH3D* eff_histo = GetBinnedEfficiency(efficiency_files_, efficiency_model_);
    RooArgSet vars(*thetat_, *thetab_, *phit_);
    RooDataHist eff_roohisto("eff", "eff", vars, eff_histo);

    thetat_->setBins(50);
    thetab_->setBins(50);
    phit_->setBins(50);

    tools::PlotVars2D(*thetat_, *thetab_, eff_roohisto);
    tools::PlotVars2D(*thetat_, *phit_, eff_roohisto);
    tools::PlotVars2D(*thetab_, *phit_, eff_roohisto);
}

const void FitterCPV::LogTextFromFile(const char* field_name, const char* filename) {
    // Set the current directory back to the one for plots (ugly ROOT stuff)
    if (output_file_) {
        output_file_->cd();
    }
    std::ifstream file;
    file.open(filename);
    std::stringstream buffer;
    buffer << file.rdbuf();
    TNamed text(field_name, buffer.str());
    text.Write();
}

const void FitterCPV::LogFileCRC(const char* field_name, const char* filename) {
    char buffer[100];
    snprintf(buffer, 100, "%lu", cksum(filename, true));
    TNamed crc(field_name, buffer);
    crc.Write();
}

const void FitterCPV::LogText(const char* field_name, const char* text) {
    TNamed text_field(field_name, text);
    text_field.Write();
}

/**
 * Save -2 log(likelihood) scan of a single variable.
 *
 * @param pdf PDF to be used for likelihood calculation
 * @param var Variable to scan
 * @param margin [optional] Set plot range [var_value - margin, var_value + margin]
 */
const void FitterCPV::SaveLikelihoodScan(RooAbsPdf& pdf, RooRealVar* var, const double margin) {
    TString name;
    name = "nll_";
    name += var->GetName();

    TCanvas canvas(name, name, 500, 500);

    const int steps = 100;
    const double orig_val = var->getVal();
    const double min = margin ? orig_val - margin : var->getMin();
    const double max = margin ? orig_val + margin : var->getMax();
    const double stepsize = (max - min) / steps;

    TH1D h1_nll("h1_" + name, "h1_" + name, steps, min, max);
    RooAbsReal* nll;
    nll = pdf.createNLL(*dataset_, RooFit::NumCPU(num_CPUs_));

    for (Int_t i = 0; i < steps; i++) {
        var->setVal(i * stepsize + min + stepsize / 2);
        Log::print(Log::debug, "Computing %i/%i likelihood function.\n", i + 1, steps);
        const double nll_val = 2 * nll->getVal();
        if (!std::isnan(nll_val)) {
            h1_nll.Fill(var->getVal(), nll_val);
        }
    }

    // Set optimal y-axis range, while ignoring empty bins
    double min_bin_content = h1_nll.GetMaximum();
    for (int i = 1; i < h1_nll.GetSize(); i++) {
        if (min_bin_content > h1_nll.GetBinContent(i) && h1_nll.GetBinContent(i) != 0) {
            min_bin_content = h1_nll.GetBinContent(i);
        }
    }
    const double ymargin = (h1_nll.GetMaximum() - min_bin_content) * 0.05;
    h1_nll.GetYaxis()->SetRangeUser(min_bin_content - ymargin, h1_nll.GetMaximum() + ymargin);

    delete nll;

    // Set the current directory back to the one for plots (ugly ROOT stuff)
    if (output_file_) {
        output_file_->cd();
    }

    h1_nll.GetXaxis()->SetTitle(var->GetTitle());
    h1_nll.GetYaxis()->SetTitle("");
    h1_nll.SetTitle("");
    h1_nll.SetStats(kFALSE);
    h1_nll.Draw("HIST");
    h1_nll.Write();
    canvas.SaveAs(constants::format);
    var->setVal(orig_val);
}

/**
 * Save -2 log(likelihood) scan of a two variables.
 *
 * @param pdf PDF to be used for likelihood calculation
 * @param var1 x-axis variable to scan
 * @param var2 y-axis variable to scan
 * @param margin1 [optional] Set plot range [var1_value - margin1, var1_value + margin1]
 * @param margin2 [optional] Set plot range [var2_value - margin2, var2_value + margin2]
 */
const void FitterCPV::SaveLikelihoodScan(RooAbsPdf& pdf, RooRealVar* var1, RooRealVar* var2,
                                         const double margin1, const double margin2) {
    TString name;
    name = "nll2_";
    name += var1->GetName();
    name += "_";
    name += var2->GetName();

    TCanvas canvas(name, name, 500, 500);

    const int steps1 = 30;
    const double orig_val1 = var1->getVal();
    const double min1 = margin1 ? orig_val1 - margin1 : var1->getMin();
    const double max1 = margin1 ? orig_val1 + margin1 : var1->getMax();
    const double stepsize1 = (max1 - min1) / steps1;

    const int steps2 = 30;
    const double orig_val2 = var2->getVal();
    const double min2 = margin2 ? orig_val2 - margin2 : var2->getMin();
    const double max2 = margin2 ? orig_val2 + margin2 : var2->getMax();
    const double stepsize2 = (max2 - min2) / steps2;

    TH2F h2_nll("h2_" + name, "h2_" + name, steps1, min1, max1, steps2, min2, max2);
    RooAbsReal* nll;
    nll = pdf.createNLL(*dataset_, RooFit::NumCPU(num_CPUs_));

    for (Int_t i = 0; i < steps1; i++) {
        for (Int_t j = 0; j < steps2; j++) {
            var1->setVal(i * stepsize1 + min1 + stepsize1 / 2);
            var2->setVal(j * stepsize2 + min2 + stepsize2 / 2);
            Log::print(Log::debug, "Computing %i/%i likelihood function.\n", i * steps2 + j + 1, steps1 * steps2);
            h2_nll.Fill(var1->getVal(), var2->getVal(), 2 * nll->getVal());
        }
    }

    delete nll;

    // Set the current directory back to the one for plots (ugly ROOT stuff)
    if (output_file_) {
        output_file_->cd();
    }

    canvas.SetRightMargin(0.14);
    h2_nll.GetXaxis()->SetTitle(var1->GetTitle());
    h2_nll.GetYaxis()->SetTitle(var2->GetTitle());
    h2_nll.GetZaxis()->SetTitle("");
    h2_nll.SetTitle("");
    h2_nll.SetStats(kFALSE);
    h2_nll.SetOption("colz");
    h2_nll.Draw();
    h2_nll.Write();
    canvas.SaveAs(constants::format);
    var1->setVal(orig_val1);
    var2->setVal(orig_val2);
}

const double FitterCPV::Calculate3DChi2(const RooDataHist& data, const RooDataHist& pdf) {
    // TODO: Remove the commented histogram lines
    // TH1I h_bin_content("h_bin_content", "Bin Content", 100, 0, 99);
    double chi2 = 0;
    int bins_used = 0;
    for(int i = 0; i < data.numEntries(); i++) {
        data.get(i);
        pdf.get(i);
        double data_bin = data.weight();
        double pdf_bin = pdf.weight();
        // h_bin_content.Fill(data_bin);
        if (data_bin > 5) {
            chi2 += (data_bin - pdf_bin)*(data_bin - pdf_bin)/pdf_bin;
            bins_used++;
        }
    }
    // h_bin_content.GetXaxis()->SetTitle("bin content");
    // h_bin_content.GetYaxis()->SetTitle("");
    // h_bin_content.SetTitle("");
    // h_bin_content.SetStats(kFALSE);
    // h_bin_content.Draw("HIST");
    // h_bin_content.Write();
    // h_bin_content.SaveAs(constants::format);
    Log::print(Log::info, "Bins used for chi2: %i/%i (%.1f%%)\n", bins_used, data.numEntries(), (double) bins_used / data.numEntries() * 100);

    return chi2;
}

/**
 * Save chi2 scan of a single variable.
 */
const void FitterCPV::SaveChi2Scan(RooSimultaneous& pdf, RooRealVar* var, const double margin) {
    TString name;
    name = "chi2_";
    name += var->GetName();

    TCanvas canvas(name, name, 500, 500);

    const int steps = 100;
    const double orig_val = var->getVal();
    const double min = margin ? orig_val - margin : var->getMin();
    const double max = margin ? orig_val + margin : var->getMax();
    const double stepsize = (max - min) / steps;

    RooDataSet* dataset_B = dynamic_cast<RooDataSet*>(
        dataset_->reduce("decaytype==decaytype::a||decaytype==decaytype::b"));
    RooDataSet* dataset_B_bar = dynamic_cast<RooDataSet*>(
        dataset_->reduce("decaytype==decaytype::ab||decaytype==decaytype::bb"));

    RooDataHist hist("hist", "hist", RooArgSet(*thetat_, *thetab_, *phit_), *dataset_);

    RooAbsPdf* pdf_B = pdf.getPdf("a");
    RooAbsPdf* pdf_B_bar = pdf.getPdf("ab");

    TH1D h1_chi2("h1_" + name, "h1_" + name, steps, min, max);

    for (Int_t i = 0; i < steps; i++) {
        var->setVal(i * stepsize + min + stepsize / 2);
        Log::print(Log::debug, "Computing chi2 %i/%i.\n", i + 1, steps);

        RooDataHist* cr_hist = pdf_B->generateBinned(RooArgSet(*thetat_, *thetab_, *phit_),
                                       dataset_B->sumEntries(), RooFit::ExpectedData(true));
        RooDataHist* cr_hist_B_bar =
            pdf_B_bar->generateBinned(RooArgSet(*thetat_, *thetab_, *phit_),
                                     dataset_B_bar->sumEntries(), RooFit::ExpectedData(true));
        cr_hist->add(*cr_hist_B_bar);

        const double chi2 = Calculate3DChi2(hist, *cr_hist);
        if (!std::isnan(chi2)) {
            h1_chi2.Fill(var->getVal(), chi2);
        }
    }

    // Set optimal y-axis range, while ignoring empty bins
    double min_bin_content = h1_chi2.GetMaximum();
    for (int i = 1; i < h1_chi2.GetSize(); i++) {
        if (min_bin_content > h1_chi2.GetBinContent(i) && h1_chi2.GetBinContent(i) != 0) {
            min_bin_content = h1_chi2.GetBinContent(i);
        }
    }
    const double ymargin = (h1_chi2.GetMaximum() - min_bin_content) * 0.05;
    h1_chi2.GetYaxis()->SetRangeUser(min_bin_content - ymargin, h1_chi2.GetMaximum() + ymargin);

    h1_chi2.GetXaxis()->SetTitle(var->GetTitle());
    h1_chi2.GetYaxis()->SetTitle("");
    h1_chi2.SetTitle("");
    h1_chi2.SetStats(kFALSE);
    h1_chi2.Draw("HIST");

    // Set the current directory back to the one for plots (ugly ROOT stuff)
    if (output_file_) {
        output_file_->cd();
    }

    h1_chi2.Write();
    canvas.SaveAs(constants::format);
    var->setVal(orig_val);
}

RooAbsPdf* FitterCPV::CreateAngularBKGPDF() {
    // Background phit model
    RooRealVar* bkg_phit_poly_p2 = new RooRealVar("bkg_phit_poly_p2", "p_(2)", 0.040);
    RooRealVar* bkg_phit_f = new RooRealVar("bkg_phit_f", "f_(poly)", 0.660);
    RooPolynomial* bkg_phit_poly =
        new RooPolynomial("bkg_phit_poly", "bkg_phit_poly", *phit_, *bkg_phit_poly_p2, 2);
    RooRealVar* bkg_phit_offset = new RooRealVar("bkg_phit_offset", "#phi_(t)^(offset)", 0.103);
    RooFormulaVar* bkg_phit_phit =
        new RooFormulaVar("bkg_phit_phit", "bkg_phit_phit", "phit - bkg_phit_offset",
                          RooArgList(*phit_, *bkg_phit_offset));
    RooGenericPdf* bkg_phit_cos = new RooGenericPdf(
        "bkg_phit_cos", "bkg_phit_cos", "cos(bkg_phit_phit)^2", RooArgList(*bkg_phit_phit));
    RooAddPdf* bkg_phit_model =
        new RooAddPdf("bkg_phit_model", "bkg_phit_model", RooArgList(*bkg_phit_poly, *bkg_phit_cos),
                      RooArgList(*bkg_phit_f));

    bkg_parameters_argset_.add(*bkg_phit_poly_p2);
    bkg_parameters_argset_.add(*bkg_phit_f);
    bkg_parameters_argset_.add(*bkg_phit_offset);

    // Background thetat model
    RooRealVar* bkg_thetat_f = new RooRealVar("bkg_thetat_f", "#theta_(t)^(w)", -0.207);
    RooFormulaVar* bkg_thetat_thetat = new RooFormulaVar(
        "bkg_thetat_thetat", "bkg_thetat_thetat", "(thetat - 1.5708)*(1+bkg_thetat_f) + 1.5708",
        RooArgList(*thetat_, *bkg_thetat_f));
    RooGenericPdf* bkg_thetat_model =
        new RooGenericPdf("bkg_thetat_model", "bkg_thetat_model", "sin(bkg_thetat_thetat)^3",
                          RooArgList(*bkg_thetat_thetat));

    bkg_parameters_argset_.add(*bkg_thetat_f);

    // Background thetab model
    RooRealVar* bkg_thetab_gaus_mu = new RooRealVar("bkg_thetab_gaus_mu", "#mu", 2.895);
    RooRealVar* bkg_thetab_gaus_sigma_l =
        new RooRealVar("bkg_thetab_gaus_sigma_l", "#sigma_(L)", 0.902);
    RooRealVar* bkg_thetab_gaus_sigma_r =
        new RooRealVar("bkg_thetab_gaus_sigma_r", "#sigma_(R)", 0.089);
    RooBifurGauss* bkg_thetab_gaus =
        new RooBifurGauss("bkg_thetab_gaus", "bkg_thetab_gaus", *thetab_, *bkg_thetab_gaus_mu,
                          *bkg_thetab_gaus_sigma_l, *bkg_thetab_gaus_sigma_r);
    RooRealVar* bkg_thetab_exp_alpha = new RooRealVar("bkg_thetab_exp_alpha", "#alpha", -2.182);
    RooExponential* bkg_thetab_exp =
        new RooExponential("bkg_thetab_exp", "bkg_thetab_exp", *thetab_, *bkg_thetab_exp_alpha);
    RooRealVar* bkg_thetab_f = new RooRealVar("bkg_thetab_f", "f_(exp)", 0.661);

    RooAddPdf* bkg_thetab_model =
        new RooAddPdf("bkg_thetab_model", "bkg_thetab_model",
                      RooArgList(*bkg_thetab_exp, *bkg_thetab_gaus), RooArgList(*bkg_thetab_f));

    bkg_parameters_argset_.add(*bkg_thetab_gaus_mu);
    bkg_parameters_argset_.add(*bkg_thetab_gaus_sigma_l);
    bkg_parameters_argset_.add(*bkg_thetab_gaus_sigma_r);
    bkg_parameters_argset_.add(*bkg_thetab_exp_alpha);
    bkg_parameters_argset_.add(*bkg_thetab_f);

    RooProdPdf* bkg_pdf = new RooProdPdf(
        "bkg_pdf", "bkg_pdf", RooArgList(*bkg_thetat_model, *bkg_thetab_model, *bkg_phit_model));

    return bkg_pdf;
}

RooAbsPdf* FitterCPV::CreateAngularSCFPDF() {
    // Self-cross-feed phit model
    RooRealVar* scf_phit_poly_p2 = new RooRealVar("scf_phit_poly_p2", "p_(2)", 0.856);
    RooRealVar* scf_phit_f = new RooRealVar("scf_phit_f", "f_(poly)", 0.147);
    RooPolynomial* scf_phit_poly =
        new RooPolynomial("scf_phit_poly", "scf_phit_poly", *phit_, *scf_phit_poly_p2, 2);
    RooRealVar* scf_phit_offset = new RooRealVar("scf_phit_offset", "#phi_(t)^(offset)", 0.056);
    RooFormulaVar* scf_phit_phit =
        new RooFormulaVar("scf_phit_phit", "scf_phit_phit", "phit - scf_phit_offset",
                          RooArgList(*phit_, *scf_phit_offset));
    RooGenericPdf* scf_phit_cos = new RooGenericPdf(
        "scf_phit_cos", "scf_phit_cos", "cos(scf_phit_phit)^2", RooArgList(*scf_phit_phit));
    RooAddPdf* scf_phit_model =
        new RooAddPdf("scf_phit_model", "scf_phit_model", RooArgList(*scf_phit_poly, *scf_phit_cos),
                      RooArgList(*scf_phit_f));

    scf_parameters_argset_.add(*scf_phit_poly_p2);
    scf_parameters_argset_.add(*scf_phit_f);
    scf_parameters_argset_.add(*scf_phit_offset);

    // Self-cross-feed thetat model
    RooRealVar* scf_thetat_f = new RooRealVar("scf_thetat_f", "#theta_(t)^(w)", -0.051);
    RooFormulaVar* scf_thetat_thetat = new RooFormulaVar(
        "scf_thetat_thetat", "scf_thetat_thetat", "(thetat - 1.5708)*(1+scf_thetat_f) + 1.5708",
        RooArgList(*thetat_, *scf_thetat_f));
    RooGenericPdf* scf_thetat_model =
        new RooGenericPdf("scf_thetat_model", "scf_thetat_model", "sin(scf_thetat_thetat)^3",
                          RooArgList(*scf_thetat_thetat));

    scf_parameters_argset_.add(*scf_thetat_f);

    // Self-cross-feed thetab model
    RooRealVar* scf_thetab_gaus_mu = new RooRealVar("scf_thetab_gaus_mu", "#mu", 2.885);
    RooRealVar* scf_thetab_gaus_sigma_l =
        new RooRealVar("scf_thetab_gaus_sigma_l", "#sigma_(L)", 0.411);
    RooRealVar* scf_thetab_gaus_sigma_r =
        new RooRealVar("scf_thetab_gaus_sigma_r", "#sigma_(R)", 0.094);
    RooBifurGauss* scf_thetab_gaus =
        new RooBifurGauss("scf_thetab_gaus", "scf_thetab_gaus", *thetab_, *scf_thetab_gaus_mu,
                          *scf_thetab_gaus_sigma_l, *scf_thetab_gaus_sigma_r);
    RooRealVar* scf_thetab_exp_alpha = new RooRealVar("scf_thetab_exp_alpha", "#alpha", -4.63);
    RooExponential* scf_thetab_exp =
        new RooExponential("scf_thetab_exp", "scf_thetab_exp", *thetab_, *scf_thetab_exp_alpha);
    RooRealVar* scf_thetab_f = new RooRealVar("scf_thetab_f", "f_(exp)", 0.625);

    RooAddPdf* scf_thetab_model =
        new RooAddPdf("scf_thetab_model", "scf_thetab_model",
                      RooArgList(*scf_thetab_exp, *scf_thetab_gaus), RooArgList(*scf_thetab_f));

    scf_parameters_argset_.add(*scf_thetab_gaus_mu);
    scf_parameters_argset_.add(*scf_thetab_gaus_sigma_l);
    scf_parameters_argset_.add(*scf_thetab_gaus_sigma_r);
    scf_parameters_argset_.add(*scf_thetab_exp_alpha);
    scf_parameters_argset_.add(*scf_thetab_f);

    RooProdPdf* scf_pdf = new RooProdPdf(
        "scf_pdf", "scf_pdf", RooArgList(*scf_thetat_model, *scf_thetab_model, *scf_phit_model));

    return scf_pdf;
}

void FitterCPV::SetSCFKDE(const char* file) {
    Log::print(Log::info, "Setting up KDE SCF model from file '%s'\n", file);
    OneDimPhaseSpace phasespace_thetat{"phasespace_thetat", thetat_->getMin(), thetat_->getMax()};
    OneDimPhaseSpace phasespace_thetab{"phasespace_thetab", thetab_->getMin(), thetab_->getMax()};
    OneDimPhaseSpace phasespace_phit{"phasespace_phit", phit_->getMin(), phit_->getMax()};
    CombinedPhaseSpace phasespace{"phasespace", &phasespace_thetat, &phasespace_thetab,
                                  &phasespace_phit};
    BinnedDensity binned_scf_kde("binned_scf_kde", &phasespace, file);

    TH3D* histo =
        new TH3D("histo", "histo", 100, thetat_->getMin(), thetat_->getMax(), 100,
                 thetab_->getMin(), thetab_->getMax(), 100, phit_->getMin(), phit_->getMax());

    double thetat, thetab, phit;
    int num_replaced = 0;
    int num_checked = 0;
    for (int x = 1; x <= histo->GetXaxis()->GetNbins(); x++) {
        for (int y = 1; y <= histo->GetYaxis()->GetNbins(); y++) {
            for (int z = 1; z <= histo->GetZaxis()->GetNbins(); z++) {
                thetat = histo->GetXaxis()->GetBinCenter(x);
                thetab = histo->GetYaxis()->GetBinCenter(y);
                phit = histo->GetZaxis()->GetBinCenter(z);
                std::vector<Double_t> coords(3);
                coords[0] = thetat;
                coords[1] = thetab;
                coords[2] = phit;
                num_checked++;
                if (CloseToEdge(coords, 0.0)) {
                    thetat_->setVal(thetat);
                    thetab_->setVal(thetab);
                    phit_->setVal(phit);
                    RooArgSet set(*thetat_, *thetab_, *phit_);
                    // scf_angular_pdf_->getObservables();
                    // Log::print(Log::debug, "Close to edge, %f replaced by %f at [%f, %f, %f]\n",
                    //            binned_scf_kde.density(coords), scf_angular_pdf_->getVal(),
                    //            thetat_->getVal(), thetab_->getVal(), phit_->getVal());
                    // histo->SetBinContent(histo->GetBin(x, y, z), scf_angular_pdf_->getVal());
                    histo->SetBinContent(histo->GetBin(x, y, z),
                                         binned_scf_kde.density(coords) * 1.5);
                    num_replaced++;
                } else {
                    histo->SetBinContent(histo->GetBin(x, y, z), binned_scf_kde.density(coords));
                }
            }
        }
    }
    Log::print(Log::info, "Replaced %i/%i (%.1f%%) SCF bins\n", num_replaced, num_checked,
               (double)num_replaced / num_checked * 100);

    scf_angular_kde_hist_ = new RooDataHist("scf_angular_kde_hist_", "scf_angular_kde_hist_",
                                            RooArgList(*thetat_, *thetab_, *phit_), histo);
    scf_angular_kde_ =
        new RooHistPdf("scf_angular_kde", "scf_angular_kde", RooArgSet(*thetat_, *thetab_, *phit_),
                       *scf_angular_kde_hist_);

    scf_angular_pdf_ = scf_angular_kde_;

    // OneDimPhaseSpace* phasespace_thetat = new OneDimPhaseSpace{"phasespace_thetat", thetat_->getMin(), thetat_->getMax()};
    // OneDimPhaseSpace* phasespace_thetab = new OneDimPhaseSpace{"phasespace_thetab", thetab_->getMin(), thetab_->getMax()};
    // OneDimPhaseSpace* phasespace_phit = new OneDimPhaseSpace{"phasespace_phit", phit_->getMin(), phit_->getMax()};
    // CombinedPhaseSpace* phasespace = new CombinedPhaseSpace{"phasespace", phasespace_thetat, phasespace_thetab,
    //                               phasespace_phit};
	// BinnedDensity* binned_scf_kde = new BinnedDensity("binned_scf_kde", phasespace, file);
    // RooArgList list(*thetat_, *thetab_, *phit_);
    // RooMeerkatPdf* meerkat_pdf =
    //     new RooMeerkatPdf("meerkat_pdf", "meerkat_pdf", list, binned_scf_kde);
    // scf_angular_pdf_ = meerkat_pdf;
}

void FitterCPV::SetSCFHisto(const char* file) {
    Log::print(Log::info, "Setting up SCF RooHistPdf model from file '%s'\n", file);
    TFile f(file, "READ");
    RooHistPdf* temp_pdf = dynamic_cast<RooHistPdf*>(f.Get("scf_hist_pdf"));

    // We create a RooHistPdfFast (our version of RooHistPdf that implements
    // caching of its "analytical integral") from the original RooHistPdf to
    // avoid the huge performance hit due to a RooFit bug.
    scf_angular_pdf_ =
        new RooHistPdfFast("scf_hist_pdf_fast", "scf_hist_pdf_fast", RooArgSet(*thetat_, *thetab_, *phit_),
                       temp_pdf->dataHist());
}

int FitterCPV::CloseToEdge(const std::vector<Double_t> vals, const double margin) const {
    RooRealVar* vars[3] = {thetat_, thetab_, phit_};
    for (int var_num = 0; var_num < 3; var_num++) {
        const double min = vars[var_num]->getMin();
        const double max = vars[var_num]->getMax();
        const double range = max - min;
        if (vals[var_num] < min + range * margin) {
            return 1;
        } else if (vals[var_num] > max - range * margin) {
            return 2;
        }
    }
    return 0;
}
