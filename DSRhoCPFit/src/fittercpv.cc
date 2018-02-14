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
#include <array>
#include <boost/filesystem.hpp>
#include <sstream>
#include <string>

// Belle includes
#include "tatami/tatami.h"

// ROOT includes
#include "RooArgSet.h"
#include "RooBifurGauss.h"
#include "RooCategory.h"
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
#include "RooTreeDataStore.h"
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
#include "constants.h"
#include "dtcppdf.h"
#include "dtscfpdf.h"

FitterCPV::FitterCPV(std::array<double, 16> par_input) {
    ap_ = new RooRealVar("ap", "ap", par_input[0], 0, 0.5);
    apa_ = new RooRealVar("apa", "apa", par_input[1], 0, 1);
    a0_ = new RooRealVar("a0", "a0", par_input[2], 0.8, 1);
    a0a_ = new RooRealVar("a0a", "a0a", 0);
    at_ = new RooFormulaVar("at", "sqrt(1-ap*ap-a0*a0)", RooArgSet(*ap_, *a0_));
    ata_ = new RooRealVar("ata", "ata", par_input[3], 2, 4);

    xp_ = new RooRealVar("xp", "xp", par_input[4], -0.4, 0.4);
    x0_ = new RooRealVar("x0", "x0", par_input[5], -0.4, 0.4);
    xt_ = new RooRealVar("xt", "xt", par_input[6], -0.4, 0.4);

    yp_ = new RooRealVar("yp", "yp", par_input[7], -0.4, 0.4);
    y0_ = new RooRealVar("y0", "y0", par_input[8], -0.4, 0.4);
    yt_ = new RooRealVar("yt", "yt", par_input[9], -0.4, 0.4);

    xpb_ = new RooRealVar("xpb", "xpb", par_input[10], -0.4, 0.4);
    x0b_ = new RooRealVar("x0b", "x0b", par_input[11], -0.4, 0.4);
    xtb_ = new RooRealVar("xtb", "xtb", par_input[12], -0.4, 0.4);

    ypb_ = new RooRealVar("ypb", "ypb", par_input[13], -0.4, 0.4);
    y0b_ = new RooRealVar("y0b", "y0b", par_input[14], -0.4, 0.4);
    ytb_ = new RooRealVar("ytb", "ytb", par_input[15], -0.4, 0.4);

    thetat_ = new RooRealVar("thetat", "thetat", 0, constants::pi);
    thetab_ = new RooRealVar("thetab", "thetab", 0.5, 2.95);
    phit_ = new RooRealVar("phit", "phit", -constants::pi, constants::pi);

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

    // TODO: beta*gamma should be computed not a constant and c should be taken from constants.cc
    dt_formula_ = new RooFormulaVar("dt", "#Deltat [ps]", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)",
                                    RooArgSet(*vrvtxz_, *vtvtxz_));
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

    // Make a copy of the input parameters for saving results, etc.
    par_input_ = par_input;

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

    num_CPUs_ = 1;
    do_lifetime_fit_ = false;
    do_mixing_fit_ = false;
    make_plots_ = false;
    perfect_tagging_ = false;

    dt_ = NULL;
    vrzerr_ = NULL;
    vtzerr_ = NULL;
}

FitterCPV::~FitterCPV() {
    if (output_file_) {
        if (output_file_->IsOpen()) {
            output_file_->Close();
        }
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

// TODO: Remove/refactor
void FitterCPV::FitSignal() {
    RooDataSet* temp_dataset = static_cast<RooDataSet*>(dataset_->reduce("evmcflag==1"));
    dataset_ = temp_dataset;

    DtCPPDF mixing_pdf_a(
        "mixing_pdf_a", "mixing_pdf_a", false, true, perfect_tagging_, efficiency_model_, *thetat_, *thetab_, *phit_,
        *ap_, *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_,
        *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
        *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF mixing_pdf_ab(
        "mixing_pdf_ab", "mixing_pdf_ab", true, true, perfect_tagging_, efficiency_model_, *thetat_, *thetab_, *phit_,
        *ap_, *apa_, *a0_, *ata_, *xpb_, *x0b_, *xtb_, *ypb_, *y0b_, *ytb_,
        *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
        *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF mixing_pdf_b(
        "mixing_pdf_b", "mixing_pdf_b", false, false, perfect_tagging_, efficiency_model_, *thetat_, *thetab_, *phit_,
        *ap_, *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_,
        *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
        *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF mixing_pdf_bb(
        "mixing_pdf_bb", "mixing_pdf_bb", true, false, perfect_tagging_, efficiency_model_, *thetat_, *thetab_, *phit_,
        *ap_, *apa_, *a0_, *ata_, *xpb_, *x0b_, *xtb_, *ypb_, *y0b_, *ytb_,
        *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
        *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    RooSimultaneous sim_pdf("sim_pdf", "sim_pdf", *decaytype_);
    sim_pdf.addPdf(mixing_pdf_a, "a");
    sim_pdf.addPdf(mixing_pdf_ab, "ab");
    sim_pdf.addPdf(mixing_pdf_b, "b");
    sim_pdf.addPdf(mixing_pdf_bb, "bb");

    dt_->setRange("dtFitRange", -15, 15);

    tau_->setConstant(true);
    dm_->setConstant(true);

    if (do_mixing_fit_) {
        result_ = sim_pdf.fitTo(*dataset_, RooFit::ConditionalObservables(conditional_vars_argset_),
                                RooFit::Minimizer("Minuit2"), RooFit::Range("dtFitRange"),
                                RooFit::Save(true), RooFit::NumCPU(num_CPUs_));
        result_->Print();

        if (make_plots_) {
            RooDataSet* dataset_a =
                static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::a"));
            PlotWithPull(*dt_, *dataset_a, mixing_pdf_a);

            RooDataSet* dataset_b =
                static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::b"));
            PlotWithPull(*dt_, *dataset_b, mixing_pdf_b);

            RooDataSet* dataset_ab =
                static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::ab"));
            PlotWithPull(*dt_, *dataset_ab, mixing_pdf_ab);

            RooDataSet* dataset_bb =
                static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::bb"));
            PlotWithPull(*dt_, *dataset_bb, mixing_pdf_bb);

            // PlotWithPull(*thetat_, *dataset_, mixing_pdf_a);
            // PlotWithPull(*thetab_, *dataset_, mixing_pdf_a);
            // PlotWithPull(*phit_, *dataset_, mixing_pdf_a);
        }
    }
}

// TODO: Remove/refactor
void FitterCPV::FitSCF() {
    RooDataSet* temp_dataset = static_cast<RooDataSet*>(dataset_->reduce("evmcflag!=1"));
    dataset_ = temp_dataset;

    DtSCFPDF mixing_pdf_a(
        "mixing_pdf_a", "mixing_pdf_a", false, true, perfect_tagging_,
        *ap_, *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_,
        *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
        *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtSCFPDF mixing_pdf_ab(
        "mixing_pdf_ab", "mixing_pdf_ab", true, true, perfect_tagging_,
        *ap_, *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_,
        *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
        *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtSCFPDF mixing_pdf_b(
        "mixing_pdf_b", "mixing_pdf_b", false, false, perfect_tagging_,
        *ap_, *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_,
        *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
        *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtSCFPDF mixing_pdf_bb(
        "mixing_pdf_bb", "mixing_pdf_bb", true, false, perfect_tagging_,
        *ap_, *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_,
        *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
        *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    // Self-cross-feed phit model
    RooRealVar scf_phit_poly_p2("scf_phit_poly_p2", "p_(2)", 0.856);
    RooRealVar scf_phit_f("scf_phit_f", "f_(poly)", 0.147);
    RooPolynomial scf_phit_poly("scf_phit_poly", "scf_phit_poly", *phit_, scf_phit_poly_p2, 2);
    RooRealVar scf_phit_offset("scf_phit_offset", "#phi_(t)^(offset)", 0.056);
    RooFormulaVar scf_phit_phit("scf_phit_phit", "scf_phit_phit", "phit - scf_phit_offset",
                                 RooArgList(*phit_, scf_phit_offset));
    RooGenericPdf scf_phit_cos("scf_phit_cos", "scf_phit_cos", "cos(scf_phit_phit)^2",
                                RooArgList(scf_phit_phit));
    RooAddPdf scf_phit_model("scf_phit_model", "scf_phit_model",
                              RooArgList(scf_phit_poly, scf_phit_cos), RooArgList(scf_phit_f));

    // Self-cross-feed thetat model
    RooRealVar scf_thetat_f("scf_thetat_f", "#theta_(t)^(w)", -0.051);
    RooFormulaVar scf_thetat_thetat("scf_thetat_thetat", "scf_thetat_thetat",
                                     "(thetat - 1.5708)*(1+scf_thetat_f) + 1.5708",
                                     RooArgList(*thetat_, scf_thetat_f));
    RooGenericPdf scf_thetat_model("scf_thetat_model", "scf_thetat_model",
                                    "sin(scf_thetat_thetat)^3", RooArgList(scf_thetat_thetat));

    // Self-cross-feed thetab model
    RooRealVar scf_thetab_gaus_mu("scf_thetab_gaus_mu", "#mu", 2.885);
    RooRealVar scf_thetab_gaus_sigma_l("scf_thetab_gaus_sigma_l", "#sigma_(L)", 0.411);
    RooRealVar scf_thetab_gaus_sigma_r("scf_thetab_gaus_sigma_r", "#sigma_(R)", 0.094);
    RooBifurGauss scf_thetab_gaus(
        "scf_thetab_gaus",  "scf_thetab_gaus",       *thetab_,
        scf_thetab_gaus_mu, scf_thetab_gaus_sigma_l, scf_thetab_gaus_sigma_r);
    RooRealVar scf_thetab_exp_alpha("scf_thetab_exp_alpha", "#alpha", -4.63);
    RooExponential scf_thetab_exp("scf_thetab_exp", "scf_thetab_exp", *thetab_,
                                   scf_thetab_exp_alpha);
    RooRealVar scf_thetab_f("scf_thetab_f", "f_(exp)", 0.625);

    RooAddPdf scf_thetab_model("scf_thetab_model", "scf_thetab_model",
                              RooArgList(scf_thetab_exp, scf_thetab_gaus), RooArgList(scf_thetab_f));
    

    RooProdPdf pdf_a("scf_pdf_a", "scf_pdf_a", 
                     RooArgList(mixing_pdf_a, scf_thetat_model, scf_thetab_model, scf_phit_model));
    RooProdPdf pdf_ab("scf_pdf_ab", "scf_pdf_ab", 
                     RooArgList(mixing_pdf_ab, scf_thetat_model, scf_thetab_model, scf_phit_model));
    RooProdPdf pdf_b("scf_pdf_b", "scf_pdf_b", 
                     RooArgList(mixing_pdf_b, scf_thetat_model, scf_thetab_model, scf_phit_model));
    RooProdPdf pdf_bb("scf_pdf_bb", "scf_pdf_bb", 
                     RooArgList(mixing_pdf_bb, scf_thetat_model, scf_thetab_model, scf_phit_model));

    RooSimultaneous sim_pdf("sim_pdf", "sim_pdf", *decaytype_);
    sim_pdf.addPdf(pdf_a, "a");
    sim_pdf.addPdf(pdf_ab, "ab");
    sim_pdf.addPdf(pdf_b, "b");
    sim_pdf.addPdf(pdf_bb, "bb");

    dt_->setRange("dtFitRange", -15, 15);

    tau_->setConstant(true);
    dm_->setConstant(true);

    if (do_mixing_fit_) {
        // result_ = sim_pdf.fitTo(*dataset_, RooFit::ConditionalObservables(conditional_vars_argset_),
        //                         RooFit::Minimizer("Minuit2"), RooFit::Range("dtFitRange"),
        //                         RooFit::Save(true), RooFit::NumCPU(num_CPUs_));
        // result_->Print();

        if (make_plots_) {
            RooDataSet* dataset_a =
                static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::a"));
            PlotWithPull(*dt_, *dataset_a, pdf_a);

            RooDataSet* dataset_b =
                static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::b"));
            PlotWithPull(*dt_, *dataset_b, pdf_b);

            RooDataSet* dataset_ab =
                static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::ab"));
            PlotWithPull(*dt_, *dataset_ab, pdf_ab);

            RooDataSet* dataset_bb =
                static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::bb"));
            PlotWithPull(*dt_, *dataset_bb, pdf_bb);

            PlotWithPull(*thetat_, *dataset_, pdf_a);
            PlotWithPull(*thetab_, *dataset_, pdf_a);
            PlotWithPull(*phit_, *dataset_, pdf_a);
        }
    }
}

void FitterCPV::FitAll() {
    DtCPPDF cr_pdf_a(
        "cr_pdf_a", "cr_pdf_a", false, true, perfect_tagging_, efficiency_model_, *thetat_, *thetab_, *phit_,
        *ap_, *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_,
        *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
        *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF cr_pdf_ab(
        "cr_pdf_ab", "cr_pdf_ab", true, true, perfect_tagging_, efficiency_model_, *thetat_, *thetab_, *phit_,
        *ap_, *apa_, *a0_, *ata_, *xpb_, *x0b_, *xtb_, *ypb_, *y0b_, *ytb_,
        *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
        *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF cr_pdf_b(
        "cr_pdf_b", "cr_pdf_b", false, false, perfect_tagging_, efficiency_model_, *thetat_, *thetab_, *phit_,
        *ap_, *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_,
        *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
        *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF cr_pdf_bb(
        "cr_pdf_bb", "cr_pdf_bb", true, false, perfect_tagging_, efficiency_model_, *thetat_, *thetab_, *phit_,
        *ap_, *apa_, *a0_, *ata_, *xpb_, *x0b_, *xtb_, *ypb_, *y0b_, *ytb_,
        *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
        *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);


    DtSCFPDF scf_dt_pdf_a(
        "scf_dt_pdf_a", "scf_dt_pdf_a", false, true, perfect_tagging_,
        *ap_, *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_,
        *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
        *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtSCFPDF scf_dt_pdf_ab(
        "scf_dt_pdf_ab", "scf_dt_pdf_ab", true, true, perfect_tagging_,
        *ap_, *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_,
        *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
        *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtSCFPDF scf_dt_pdf_b(
        "scf_dt_pdf_b", "scf_dt_pdf_b", false, false, perfect_tagging_,
        *ap_, *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_,
        *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
        *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtSCFPDF scf_dt_pdf_bb(
        "scf_dt_pdf_bb", "scf_dt_pdf_bb", true, false, perfect_tagging_,
        *ap_, *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_,
        *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
        *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    // Self-cross-feed phit model
    RooRealVar scf_phit_poly_p2("scf_phit_poly_p2", "p_(2)", 0.856);
    RooRealVar scf_phit_f("scf_phit_f", "f_(poly)", 0.147);
    RooPolynomial scf_phit_poly("scf_phit_poly", "scf_phit_poly", *phit_, scf_phit_poly_p2, 2);
    RooRealVar scf_phit_offset("scf_phit_offset", "#phi_(t)^(offset)", 0.056);
    RooFormulaVar scf_phit_phit("scf_phit_phit", "scf_phit_phit", "phit - scf_phit_offset",
                                 RooArgList(*phit_, scf_phit_offset));
    RooGenericPdf scf_phit_cos("scf_phit_cos", "scf_phit_cos", "cos(scf_phit_phit)^2",
                                RooArgList(scf_phit_phit));
    RooAddPdf scf_phit_model("scf_phit_model", "scf_phit_model",
                              RooArgList(scf_phit_poly, scf_phit_cos), RooArgList(scf_phit_f));

    // Self-cross-feed thetat model
    RooRealVar scf_thetat_f("scf_thetat_f", "#theta_(t)^(w)", -0.051);
    RooFormulaVar scf_thetat_thetat("scf_thetat_thetat", "scf_thetat_thetat",
                                     "(thetat - 1.5708)*(1+scf_thetat_f) + 1.5708",
                                     RooArgList(*thetat_, scf_thetat_f));
    RooGenericPdf scf_thetat_model("scf_thetat_model", "scf_thetat_model",
                                    "sin(scf_thetat_thetat)^3", RooArgList(scf_thetat_thetat));

    // Self-cross-feed thetab model
    RooRealVar scf_thetab_gaus_mu("scf_thetab_gaus_mu", "#mu", 2.885);
    RooRealVar scf_thetab_gaus_sigma_l("scf_thetab_gaus_sigma_l", "#sigma_(L)", 0.411);
    RooRealVar scf_thetab_gaus_sigma_r("scf_thetab_gaus_sigma_r", "#sigma_(R)", 0.094);
    RooBifurGauss scf_thetab_gaus(
        "scf_thetab_gaus",  "scf_thetab_gaus",       *thetab_,
        scf_thetab_gaus_mu, scf_thetab_gaus_sigma_l, scf_thetab_gaus_sigma_r);
    RooRealVar scf_thetab_exp_alpha("scf_thetab_exp_alpha", "#alpha", -4.63);
    RooExponential scf_thetab_exp("scf_thetab_exp", "scf_thetab_exp", *thetab_,
                                   scf_thetab_exp_alpha);
    RooRealVar scf_thetab_f("scf_thetab_f", "f_(exp)", 0.625);

    RooAddPdf scf_thetab_model("scf_thetab_model", "scf_thetab_model",
                              RooArgList(scf_thetab_exp, scf_thetab_gaus), RooArgList(scf_thetab_f));
    

    RooProdPdf scf_pdf_a("scf_pdf_a", "scf_pdf_a", 
                         RooArgList(scf_dt_pdf_a, scf_thetat_model, scf_thetab_model, scf_phit_model));
    RooProdPdf scf_pdf_ab("scf_pdf_ab", "scf_pdf_ab", 
                         RooArgList(scf_dt_pdf_ab, scf_thetat_model, scf_thetab_model, scf_phit_model));
    RooProdPdf scf_pdf_b("scf_pdf_b", "scf_pdf_b", 
                         RooArgList(scf_dt_pdf_b, scf_thetat_model, scf_thetab_model, scf_phit_model));
    RooProdPdf scf_pdf_bb("scf_pdf_bb", "scf_pdf_bb", 
                     RooArgList(scf_dt_pdf_bb, scf_thetat_model, scf_thetab_model, scf_phit_model));
    

    RooRealVar cr_scf_f("cr_scf_f", "f_{cr}", 0.8);

    RooAddPdf pdf_a("pdf_a", "pdf_a", RooArgList(cr_pdf_a, scf_pdf_a), RooArgList(cr_scf_f));
    RooAddPdf pdf_ab("pdf_ab", "pdf_ab", RooArgList(cr_pdf_ab, scf_pdf_ab), RooArgList(cr_scf_f));
    RooAddPdf pdf_b("pdf_b", "pdf_b", RooArgList(cr_pdf_b, scf_pdf_b), RooArgList(cr_scf_f));
    RooAddPdf pdf_bb("pdf_bb", "pdf_bb", RooArgList(cr_pdf_bb, scf_pdf_bb), RooArgList(cr_scf_f));

    RooSimultaneous sim_pdf("sim_pdf", "sim_pdf", *decaytype_);
    sim_pdf.addPdf(pdf_a, "a");
    sim_pdf.addPdf(pdf_ab, "ab");
    sim_pdf.addPdf(pdf_b, "b");
    sim_pdf.addPdf(pdf_bb, "bb");

    dt_->setRange("dtFitRange", -15, 15);

    tau_->setConstant(true);
    dm_->setConstant(true);

    if (do_mixing_fit_) {
        result_ = sim_pdf.fitTo(*dataset_, RooFit::ConditionalObservables(conditional_vars_argset_),
                                RooFit::Minimizer("Minuit2"), RooFit::Save(true), RooFit::NumCPU(num_CPUs_));
        result_->Print();

        if (make_plots_) {
            RooDataSet* cr_dataset_ =
                static_cast<RooDataSet*>(dataset_->reduce("evmcflag==1"));
            RooDataSet* scf_dataset_ =
                static_cast<RooDataSet*>(dataset_->reduce("evmcflag!=1"));


            RooDataSet* cr_dataset_a =
                static_cast<RooDataSet*>(cr_dataset_->reduce("decaytype==decaytype::a"));
            PlotWithPull(*dt_, *cr_dataset_a, cr_pdf_a);

            RooDataSet* cr_dataset_b =
                static_cast<RooDataSet*>(cr_dataset_->reduce("decaytype==decaytype::b"));
            PlotWithPull(*dt_, *cr_dataset_b, cr_pdf_b);

            RooDataSet* cr_dataset_ab =
                static_cast<RooDataSet*>(cr_dataset_->reduce("decaytype==decaytype::ab"));
            PlotWithPull(*dt_, *cr_dataset_ab, cr_pdf_ab);

            RooDataSet* cr_dataset_bb =
                static_cast<RooDataSet*>(cr_dataset_->reduce("decaytype==decaytype::bb"));
            PlotWithPull(*dt_, *cr_dataset_bb, cr_pdf_bb);


            RooDataSet* scf_dataset_a =
                static_cast<RooDataSet*>(scf_dataset_->reduce("decaytype==decaytype::a"));
            PlotWithPull(*dt_, *scf_dataset_a, scf_pdf_a);

            RooDataSet* scf_dataset_b =
                static_cast<RooDataSet*>(scf_dataset_->reduce("decaytype==decaytype::b"));
            PlotWithPull(*dt_, *scf_dataset_b, scf_pdf_b);

            RooDataSet* scf_dataset_ab =
                static_cast<RooDataSet*>(scf_dataset_->reduce("decaytype==decaytype::ab"));
            PlotWithPull(*dt_, *scf_dataset_ab, scf_pdf_ab);

            RooDataSet* scf_dataset_bb =
                static_cast<RooDataSet*>(scf_dataset_->reduce("decaytype==decaytype::bb"));
            PlotWithPull(*dt_, *scf_dataset_bb, scf_pdf_bb);


            // PlotWithPull(*thetat_, *cr_dataset_, cr_pdf_a);
            // PlotWithPull(*thetab_, *cr_dataset_, cr_pdf_a);
            // PlotWithPull(*phit_, *cr_dataset_, cr_pdf_a);
            PlotWithPull(*thetat_, *scf_dataset_, scf_pdf_a);
            PlotWithPull(*thetab_, *scf_dataset_, scf_pdf_a);
            PlotWithPull(*phit_, *scf_dataset_, scf_pdf_a);
            // thetab_->setBins(300);
            RooDataHist* cr_hist = cr_pdf_a.generateBinned(RooArgSet(*thetat_, *thetab_, *phit_), 1000*0.8, RooFit::ExpectedData(true));
            RooDataHist* scf_hist = scf_pdf_a.generateBinned(RooArgSet(*thetat_, *thetab_, *phit_), 1000*0.2, RooFit::ExpectedData(true));

            // RooHistPdf all_histpdf("all_histpdf", "all_histpdf", RooArgSet(*thetat_, *thetab_, *phit_), *all_hist);
            RooHistPdf cr_histpdf("cr_histpdf", "cr_histpdf", RooArgSet(*thetat_, *thetab_, *phit_), *cr_hist);
            RooHistPdf scf_histpdf("scf_histpdf", "scf_histpdf", RooArgSet(*thetat_, *thetab_, *phit_), *scf_hist);
            RooAddPdf all_histpdf("all_histpdf", "all_histpdf", RooArgList(cr_histpdf, scf_histpdf), cr_scf_f);
            std::vector<RooAbsPdf*> components;
            components.push_back(&cr_histpdf);
            components.push_back(&scf_histpdf);
            // thetab_->setBins(100);
            PlotWithPull(*thetat_, *dataset_, all_histpdf, components);
            PlotWithPull(*thetab_, *dataset_, all_histpdf, components);
            PlotWithPull(*phit_, *dataset_, all_histpdf, components);
        }
    }
}

void FitterCPV::GenerateToys(const int num_events, const int num_toys) {
    printf("INFO: Generating dataset with %i events...\n", num_events);

    // vars with phiweak = phiweak + pi, r = 0.01
    double par_input[] = {
        0.269,
        0.56,
        0.941,
        3.11,

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

    DtCPPDF mixing_pdf_a(
        "mixing_pdf_a", "mixing_pdf_a", false, true, perfect_tagging_, efficiency_model_, *thetat_, *thetab_, *phit_,
        ap, apa, a0, ata, xp, x0, xt, yp, y0, yt,
        *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
        *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF mixing_pdf_ab(
        "mixing_pdf_ab", "mixing_pdf_ab", true, true, perfect_tagging_,  efficiency_model_,*thetat_, *thetab_, *phit_,
        ap, apa, a0, ata, xpb, x0b, xtb, ypb, y0b, ytb,
        *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
        *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF mixing_pdf_b(
        "mixing_pdf_b", "mixing_pdf_b", false, false, perfect_tagging_, efficiency_model_, *thetat_, *thetab_, *phit_,
        ap, apa, a0, ata, xp, x0, xt, yp, y0, yt,
        *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
        *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF mixing_pdf_bb(
        "mixing_pdf_bb", "mixing_pdf_bb", true, false, perfect_tagging_, efficiency_model_, *thetat_, *thetab_, *phit_,
        ap, apa, a0, ata, xpb, x0b, xtb, ypb, y0b, ytb,
        *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
        *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

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

    printf("INFO: Creating GenSpecs...\n");
    RooAbsPdf::GenSpec* genSpec_a =
        mixing_pdf_a.prepareMultiGen(varsToGenerate, RooFit::NumEvents(num_fav));
    printf("INFO: 1/4 GenSpecs ready\n");
    RooAbsPdf::GenSpec* genSpec_ab =
        mixing_pdf_ab.prepareMultiGen(varsToGenerate, RooFit::NumEvents(num_fav));
    printf("INFO: 2/4 GenSpecs ready\n");
    RooAbsPdf::GenSpec* genSpec_b =
        mixing_pdf_b.prepareMultiGen(varsToGenerate, RooFit::NumEvents(num_sup));
    printf("INFO: 3/4 GenSpecs ready\n");
    RooAbsPdf::GenSpec* genSpec_bb =
        mixing_pdf_bb.prepareMultiGen(varsToGenerate, RooFit::NumEvents(num_sup));
    printf("INFO: 4/4 GenSpecs ready\n");

    RooAbsData::setDefaultStorageType(RooAbsData::Tree);
    RooDataSet* dataset;
    TString filename;

    RooRealVar btagmcli("btagmcli", "btagmcli", -511, 511);
    RooRealVar brecflav("brecflav", "brecflav", -1, 1);

    for (int i = 1; i <= num_toys; i++) {
        //  RooDataSet* dataset_a = mixing_pdf_a.generate(varsToGenerate, RooFit::NumEvents(num_fav));
        RooDataSet* dataset_a = mixing_pdf_a.generate(*genSpec_a);
        decaytype_->setLabel("a");
        dataset_a->addColumn(*decaytype_);
        brecflav.setVal(1);
        btagmcli.setVal(-511);
        dataset_a->addColumn(brecflav);
        dataset_a->addColumn(btagmcli);
        dataset = (RooDataSet*)dataset_a->Clone();
        delete dataset_a;
        printf("INFO: 1/4 decay type ready.\n");

        //  RooDataSet* dataset_ab = mixing_pdf_ab.generate(varsToGenerate, RooFit::NumEvents(num_fav));
        RooDataSet* dataset_ab = mixing_pdf_ab.generate(*genSpec_ab);
        decaytype_->setLabel("ab");
        dataset_ab->addColumn(*decaytype_);
        brecflav.setVal(-1);
        btagmcli.setVal(511);
        dataset_ab->addColumn(brecflav);
        dataset_ab->addColumn(btagmcli);
        dataset->append(*dataset_ab);
        delete dataset_ab;
        printf("INFO: 2/4 decay types ready.\n");

        //  RooDataSet* dataset_b = mixing_pdf_b.generate(varsToGenerate, RooFit::NumEvents(num_sup));
        RooDataSet* dataset_b = mixing_pdf_b.generate(*genSpec_b);
        decaytype_->setLabel("b");
        dataset_b->addColumn(*decaytype_);
        brecflav.setVal(1);
        btagmcli.setVal(511);
        dataset_b->addColumn(brecflav);
        dataset_b->addColumn(btagmcli);
        dataset->append(*dataset_b);
        delete dataset_b;
        printf("INFO: 3/4 decay types ready.\n");

        //  RooDataSet* dataset_bb = mixing_pdf_bb.generate(varsToGenerate, RooFit::NumEvents(num_sup));
        RooDataSet* dataset_bb = mixing_pdf_bb.generate(*genSpec_bb);
        decaytype_->setLabel("bb");
        dataset_bb->addColumn(*decaytype_);
        brecflav.setVal(-1);
        btagmcli.setVal(-511);
        dataset_bb->addColumn(brecflav);
        dataset_bb->addColumn(btagmcli);
        dataset->append(*dataset_bb);
        delete dataset_bb;
        printf("INFO: 4/4 decay types ready.\n");

        dataset->addColumns(varsToAdd);

        dt_formula_ =
            new RooFormulaVar("dt", "#Deltat [ps]", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)",
                              RooArgSet(*vrvtxz_, *vtvtxz_));
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
    //RooFit::ConditionalObservables(conditional_vars_argset_),
    //          RooFit::Minimizer("Minuit2"), RooFit::Hesse(false), RooFit::Minos(false),
    //          RooFit::Range("dtFitRange"), RooFit::Save(true), RooFit::NumCPU(num_CPUs_));
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

    RooPlot* plot = var.frame(RooFit::Range("dtFitRange"));
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
    if (!name.Contains("hist") && num_CPUs_ > 1){
        norm = 1.0 / data.numEntries();
    }

    // Plot components before the total PDF as the pull plots are made from the
    // last plotted PDF
    const int colors[] = {4, 7, 5};
    int i = 0;
    for (auto component : components) {
        pdf.plotOn(plot, RooFit::ProjWData(conditional_vars_argset_, data, kFALSE),
                   RooFit::NumCPU(num_CPUs_), RooFit::LineColor(colors[i++]),
                   RooFit::NormRange("dtFitRange"), RooFit::Components(*component),
                   RooFit::Normalization(norm));
    }

    pdf.plotOn(plot, RooFit::ProjWData(conditional_vars_argset_, data, kFALSE),
               RooFit::NumCPU(num_CPUs_), RooFit::NormRange("dtFitRange"),
               RooFit::Normalization(norm));
    
    plot->GetXaxis()->SetTitle("");
    plot->GetXaxis()->SetLabelSize(0);

    // This line makes sure the 0 is not drawn as it would overlap with the lower pad
    plot->GetYaxis()->SetRangeUser(0.001, plot->GetMaximum());
    plot->SetTitle("");
    plot->GetYaxis()->SetTitle(title);
    plot->GetYaxis()->SetTitleOffset(1.60);
    plot->Draw();

    const double chi2 = plot->chiSquare();
    TPaveText* stat_box = CreateStatBox(chi2, true, true);
    if (stat_box) {
        stat_box->Draw();
    }

    pad_pull->cd();
    pad_pull->SetTopMargin(0.0);
    pad_pull->SetBottomMargin(0.35);
    pad_pull->SetLeftMargin(0.12);

    // Create a new frame to draw the pull distribution and add the distribution to the frame
    RooPlot* plot_pull_ =
        var.frame(RooFit::Title("Pull Distribution"), RooFit::Range("dtFitRange"));
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
    plot_pull_->addPlotable( hpull_clone, "P");  

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
 * Create and return statbox with fit results to overlay on plots. Returns a NULL pointer
 * if no fit result exist.
 *
 * @param chi2 Reduced chi2 of the projection plot
 * @param position_top Should the box be displayed at the top or bottom of the plot
 * @param position_left Should the box be displayed at the left or right of the plot
 */
TPaveText* FitterCPV::CreateStatBox(const double chi2, const bool position_top,
                                    const bool position_left) const {
    // If no fit result exists return a null pointer
    if (!result_) {
        printf("WARNING: No result exists, can't create stat box!\n");
        return NULL;
    }

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

/**
 * Construct and return a cut string common for all for categories
 */
TString FitterCPV::GetCommonCutsString() const {
    TString common_cuts("vrusable==1&&vtusable==1&&");
    common_cuts += "((vrchi2/vrndf)<";
    common_cuts += constants::cuts::sig_vtx_h;
    common_cuts += "||vrntrk==1)&&";
    common_cuts += "((vtchi2/vtndf)<";
    common_cuts += constants::cuts::tag_vtx_h;
    common_cuts += "||vtntrk==1)&&";
    common_cuts += "((sqrt(vrerr6)<";
    common_cuts += constants::cuts::sig_vtx_multitrack_sigma_z;
    common_cuts += "&&vrntrk>1)||(sqrt(vrerr6)<";
    common_cuts += constants::cuts::sig_vtx_singletrack_sigma_z;
    common_cuts += "&&vrntrk==1))";
    common_cuts += "&&";
    common_cuts += "((sqrt(vterr6)<";
    common_cuts += constants::cuts::tag_vtx_multitrack_sigma_z;
    common_cuts += "&&vtntrk>1)||(sqrt(vterr6)<";
    common_cuts += constants::cuts::tag_vtx_singletrack_sigma_z;
    common_cuts += "&&vtntrk==1))";
    return common_cuts;
}

/**
 * Reads in data from a ROOT file. Constructs separate datasets for the 4 categories.
 * Binds the variables to the dataset, so that dataset->get(i) changes values of, e.g., expno_
 *
 * @param file_path Path to the ROOT file
 * @param num_events [optional] Maximum number of events to use (0 to read all)
 */
void FitterCPV::ReadInFile(const char* file_path, const int& num_events) {
    TFile* input_file = new TFile(file_path);
    TTree* input_tree = dynamic_cast<TTree*>(input_file->Get("h2000"));
    //  TTree* input_tree = dynamic_cast<TTree*>(input_file->Get("mixing_pdf_aData"));

    if (num_events) {
        TTree* temp_tree = input_tree;
        input_tree = temp_tree->CloneTree(num_events);
        delete temp_tree;
    }

    TString common_cuts = GetCommonCutsString();

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
        a_cuts = "brecflav==1&&tagqr<0";
        ab_cuts = "brecflav==-1&&tagqr>0";
        b_cuts = "brecflav==1&&tagqr>0";
        bb_cuts = "brecflav==-1&&tagqr<0";
    }

    // A temporary RooDataSet is created from the whole tree and then we apply cuts to get
    // the 4 different B and f flavor datasets, as that is faster then reading the tree 4 times
    RooDataSet* temp_dataset =
        new RooDataSet("dataset", "dataset", input_tree, dataset_vars_argset_, common_cuts);
    
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

    dataset_ = static_cast<RooDataSet*>(dataset_a->Clone());
    delete dataset_a;
    dataset_->append(*dataset_ab);
    delete dataset_ab;
    dataset_->append(*dataset_b);
    delete dataset_b;
    dataset_->append(*dataset_bb);
    delete dataset_bb;

    delete input_tree;
    input_file->Close();

    vrzerr_ = static_cast<RooRealVar*>(dataset_->addColumn(*vrzerr_formula_));
    vrzerr_->setRange(0, 10000);
    vtzerr_ = static_cast<RooRealVar*>(dataset_->addColumn(*vtzerr_formula_));
    vtzerr_->setRange(0, 10000);
    dt_ = static_cast<RooRealVar*>(dataset_->addColumn(*dt_formula_));
    dt_->setRange(-15, 15);

    conditional_vars_argset_.add(*vrzerr_);
    conditional_vars_argset_.add(*vtzerr_);
    conditional_vars_argset_.add(*decaytype_);

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
}

/**
 * Set the directory to which to ouput plots
 */
void FitterCPV::SetPlotDir(const char* plot_dir) {
    make_plots_ = true;
    if (!boost::filesystem::is_directory(plot_dir)) {
        boost::filesystem::create_directories(plot_dir);
    }
    gEnv->SetValue("Canvas.PrintDirectory", plot_dir);
    printf("print dir: %s\n", gEnv->GetValue("Canvas.PrintDirectory", "not found"));
    output_file_ = new TFile(TString(plot_dir) + "/plots.root", "RECREATE");
}

/**
 * Fix requested parameters.
 */
bool FitterCPV::FixParameters(const char* pars) {
    std::string input = pars;
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
            printf("ERROR: Parameter %s doesn't exist.\n", par.c_str());
            return true;
        }
    }
    return false;
}


/**
 * Save results and initial values into a file.
 */
 bool FitterCPV::SaveResults(const char* file) {
     // const RooArgSet* args = dataSet->get();
     // RooPlot* frame = 0;

     const int numParameters = 54;
     double* parameters = new double[numParameters];

     FILE* pFile;
     pFile = fopen(file, "w");
     if (pFile == NULL) {
         printf("ERROR: couldn't open file %s for writing!\n", file);
         delete[] parameters;
         return 1;
     }

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

     Int_t separators[] = {0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2,
                           0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2,
                           0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2};

     for (Int_t i = 0; i < numParameters; i++) {
         /// These parameters can be both positive and negative therefore the
         /// added + in front of positive numbers keeps the columns aligned
         if (i >= 18 && (i - 20) % 3 != 0) {
             fprintf(pFile, "%+.5f ", parameters[i]);
         } else {
             fprintf(pFile, "%.5f ", parameters[i]);
         }
         if (i == numParameters - 1) continue;
         if (separators[i] == 1)
             fprintf(pFile, "| ");
         else if (separators[i] == 2)
             fprintf(pFile, "|| ");
     }
     fprintf(pFile, "\n");
     fclose(pFile);
     delete[] parameters;

     return 0;
 }
