/**
 *  @file    fittercpv.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2016-11-03
 *
 *  Class that performs the CP fit itself as well as plotting
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
#include "RooHistPdfFast.h"
#include "angularpdf.h"
#include "config.h"
#include "constants.h"
#include "dtcppdf.h"
#include "dtscfpdf.h"
#include "log.h"
#include "tools.h"

FitterCPV::FitterCPV(nlohmann::json config) {
    do_lifetime_fit_ = false;
    make_plots_ = false;
    perfect_tagging_ = false;
    generator_level_ = false;
    efficiency_model_ = -1;

    vrzerr_ = nullptr;
    vtzerr_ = nullptr;

    const int var_bins = config.contains("plotBins1D") ? config["plotBins1D"].get<int>() : 100;
    if (config.contains("initialPars")) {
        InitVars(ToParInputArray(config["initialPars"]), var_bins);
    } else {
        InitVars(constants::par_input, var_bins);
    }

    if (config.contains("fitRanges")) {
        ChangeFitRanges(config["fitRanges"]);
    }
    if (config.contains("fixedParameters")) {
        FixParameters(config["fixedParameters"].get<std::string>().c_str());
    }

    data_ = GetData(config);
    pdf_ = CreatePDF(config);
}

FitterCPV::~FitterCPV() {
    if (output_file_) {
        if (output_file_->IsOpen()) {
            output_file_->Close();
        }
    }
}

/**
 * Initialize all the member variables to reasonable values.
 */
void FitterCPV::InitVars(std::array<double, 16> par_input, const int var_bins) {
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
        thetab_ = new RooRealVar("thetabg", "#theta_{b}^{g} [rad]", constants::cuts::thetab_low,
                                 constants::cuts::thetab_high);
        phit_ = new RooRealVar("phitg", "#phi_{t}^{g} [rad]", -TMath::Pi(), TMath::Pi());
    } else {
        thetat_ = new RooRealVar("thetat", "#theta_{t} [rad]", 0, TMath::Pi());
        thetab_ = new RooRealVar("thetab", "#theta_{b} [rad]", constants::cuts::thetab_low,
                                 constants::cuts::thetab_high);
        phit_ = new RooRealVar("phit", "#phi_{t} [rad]", -TMath::Pi(), TMath::Pi());
    }

    thetat_->setBins(var_bins);
    thetab_->setBins(var_bins);
    phit_->setBins(var_bins);
    dt_->setBins(var_bins);

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

    vrzerr_formula_ =
        new RooFormulaVar("vrzerr", "#sigma z_{rec} [cm]", "sqrt(vrerr6)", RooArgSet(*vrerr6_));
    vtzerr_formula_ =
        new RooFormulaVar("vtzerr", "#sigma z_{tag} [cm]", "sqrt(vterr6)", RooArgSet(*vterr6_));

    tau_ = new RooRealVar("tau", "#tau", constants::tau - 1, constants::tau + 1);
    dm_ = new RooRealVar("dm", "#Deltam", constants::dm - 1, constants::dm + 1);

    decaytype_ = new RooCategory("decaytype", "decaytype");
    decaytype_->defineType("FB", 1);
    decaytype_->defineType("FA", 2);
    decaytype_->defineType("SB", 3);
    decaytype_->defineType("SA", 4);

    channel_cat_ = new RooCategory("channel_cat", "channel_cat");

    // Make a copy of the input parameters for saving results, etc.
    par_input_ = par_input;

    PrepareVarArgsets();
}

/**
 * Populate various vectors and argsets with all parameters.
 *
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

void FitterCPV::CreateDtCPPDFs(DtCPPDF*& cr_pdf_FB, DtCPPDF*& cr_pdf_FA, DtCPPDF*& cr_pdf_SB,
                               DtCPPDF*& cr_pdf_SA, const std::string channel_name,
                               const nlohmann::json common_config,
                               const nlohmann::json channel_config) const {
    const bool perfect_tagging = common_config.contains("perfectTagging");
    const int efficiency_model = channel_config["efficiencyModel"].get<int>();

    // The possibility to supply more than one efficiency file is a remnant of
    // hybrid histo-KDE models. I'm reluctant to remove it as it might be useful
    // in the systematics estimation.
    std::vector<std::string> efficiency_files;
    efficiency_files.push_back(channel_config["efficiencyFile"]);

    cr_pdf_FB = new DtCPPDF(
        "cr_pdf_FB", "cr_pdf_FB", false, true, perfect_tagging, efficiency_model, efficiency_files,
        *thetat_, *thetab_, *phit_, *ap_, *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_,
        *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
        *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    cr_pdf_FA = new DtCPPDF("cr_pdf_FA", "cr_pdf_FA", true, true, perfect_tagging, efficiency_model,
                            efficiency_files, *thetat_, *thetab_, *phit_, *ap_, *apa_, *a0_, *ata_,
                            *xpb_, *x0b_, *xtb_, *ypb_, *y0b_, *ytb_, *tagwtag_, *dt_, *tau_, *dm_,
                            *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_, *vrzerr_,
                            *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    cr_pdf_SB = new DtCPPDF(
        "cr_pdf_SB", "cr_pdf_SB", false, false, perfect_tagging, efficiency_model, efficiency_files,
        *thetat_, *thetab_, *phit_, *ap_, *apa_, *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_,
        *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
        *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    cr_pdf_SA = new DtCPPDF(
        "cr_pdf_SA", "cr_pdf_SA", true, false, perfect_tagging, efficiency_model, efficiency_files,
        *thetat_, *thetab_, *phit_, *ap_, *apa_, *a0_, *ata_, *xpb_, *x0b_, *xtb_, *ypb_, *y0b_,
        *ytb_, *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_,
        *vrntrk_, *vrzerr_, *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);
}

void FitterCPV::CreateDtSCFPDFs(DtSCFPDF*& scf_dt_pdf_FB, DtSCFPDF*& scf_dt_pdf_FA,
                                DtSCFPDF*& scf_dt_pdf_SB, DtSCFPDF*& scf_dt_pdf_SA,
                                const std::string channel_name) const {
    scf_dt_pdf_FB =
        new DtSCFPDF("scf_dt_pdf_FB", "scf_dt_pdf_FB", false, true, perfect_tagging_, *ap_, *apa_,
                     *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_, *tagwtag_, *dt_, *tau_, *dm_,
                     *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_, *vrzerr_, *vrchi2_,
                     *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    scf_dt_pdf_FA =
        new DtSCFPDF("scf_dt_pdf_FA", "scf_dt_pdf_FA", true, true, perfect_tagging_, *ap_, *apa_,
                     *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_, *tagwtag_, *dt_, *tau_, *dm_,
                     *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_, *vrzerr_, *vrchi2_,
                     *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    scf_dt_pdf_SB =
        new DtSCFPDF("scf_dt_pdf_SB", "scf_dt_pdf_SB", false, false, perfect_tagging_, *ap_, *apa_,
                     *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_, *tagwtag_, *dt_, *tau_, *dm_,
                     *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_, *vrzerr_, *vrchi2_,
                     *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    scf_dt_pdf_SA =
        new DtSCFPDF("scf_dt_pdf_SA", "scf_dt_pdf_SA", true, false, perfect_tagging_, *ap_, *apa_,
                     *a0_, *ata_, *xp_, *x0_, *xt_, *yp_, *y0_, *yt_, *tagwtag_, *dt_, *tau_, *dm_,
                     *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_, *vrzerr_, *vrchi2_,
                     *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);
}

/**
 * Create a Voigtian + Gaussian PDF.
 *
 * A common prefix is prepended to all object names.
 *
 * @param prefix Text to be prepended to the ROOT name
 *
 */
RooAddPdf* FitterCPV::CreateVoigtGaussDtPdf(const std::string prefix) {
    TString pre(prefix);

    RooRealVar* voigt_mu = new RooRealVar(pre + "_voigt_mu", "v_{#mu}", -0.249);
    RooRealVar* voigt_sigma = new RooRealVar(pre + "_voigt_sigma_", "v_{#sigma}", 1.910);
    RooRealVar* voigt_width = new RooRealVar(pre + "_voigt_width_", "v_{w}", 0.684);
    RooVoigtian* voigt = new RooVoigtian(pre + "_voigt", pre + "_voigt", *dt_, *voigt_mu,
                                         *voigt_width, *voigt_sigma);

    RooRealVar* gaus_mu = new RooRealVar(pre + "_gaus_mu", "g_{#mu}", -0.105);
    RooRealVar* gaus_sigma = new RooRealVar(pre + "_gaus_sigma_", "g_{#sigma}", 0.881);
    RooGaussian* gaus = new RooGaussian(pre + "_gaus", pre + "_gaus", *dt_, *gaus_mu, *gaus_sigma);

    RooRealVar* f = new RooRealVar(pre + "_f", "f_{v/g}", 0.720);
    RooAddPdf* model =
        new RooAddPdf(pre + "_model", pre + "_model", RooArgList(*voigt, *gaus), RooArgList(*f));

    return model;
}

/**
 * Create SCF or BKG PDFs with a function-based dt distribution model.
 *
 * @param pdf_FB PDF for B0 -> D*- rho+
 * @param pdf_FA PDF for anti-B0 -> D*+ rho-
 * @param pdf_SB PDF for B0 -> D*+ rho-
 * @param pdf_SA PDF for anti-B0 -> D*- rho+
 * @param channel_name Name of the channel the PDF is meant for
 * @param channel_config JSON object with configuration details (useful for SCF PDF)
 * @param scf Whether the PDF to be created is for SCF (true) or BKG (false)
 *
 */
void FitterCPV::CreateFunctionalDtSCFBKGPDFs(RooProdPdf*& pdf_FB, RooProdPdf*& pdf_FA,
                                             RooProdPdf*& pdf_SB, RooProdPdf*& pdf_SA,
                                             const std::string channel_name,
                                             const nlohmann::json channel_config, const bool scf) {
    Log::print(Log::debug, "Creating SCF dt PDF\n");

    RooAbsPdf* angular_pdf;
    if (scf) {
        angular_pdf = CreateSCFPDF(channel_name, channel_config);
    } else {
        angular_pdf = CreateAngularSCFBKGPDF(channel_name + "_bkg_");
    }

    const std::string type = scf ? "scf_" : "bkg_";
    RooAbsPdf* dt_cf_pdf = CreateVoigtGaussDtPdf(channel_name + "_" + type + "dt_cf");
    RooAbsPdf* dt_dcs_pdf = CreateVoigtGaussDtPdf(channel_name + "_" + type + "dt_dcs");

    TString pfx = channel_name;
    pfx += "_";
    pfx += type;

    pdf_FB = new RooProdPdf(pfx + "pdf_FB", pfx + "pdf_FB", RooArgList(*dt_cf_pdf, *angular_pdf));
    pdf_FA = new RooProdPdf(pfx + "pdf_FA", pfx + "pdf_FA", RooArgList(*dt_cf_pdf, *angular_pdf));
    pdf_SB = new RooProdPdf(pfx + "pdf_SB", pfx + "pdf_SB", RooArgList(*dt_dcs_pdf, *angular_pdf));
    pdf_SA = new RooProdPdf(pfx + "pdf_SA", pfx + "pdf_SA", RooArgList(*dt_dcs_pdf, *angular_pdf));
}

/**
 * Create a time-independent PDF with the requested components
 *
 * @param scf Whether to add self-crossfeed component
 * @param bkg Whether to add background component
 * @return RooSimultaneous* The complete PDF
 */
RooSimultaneous* FitterCPV::CreateAngularPDF(const std::string name_prefix, const bool scf,
                                             const bool bkg, const nlohmann::json channel_config) {
    RooArgList B_pdfs;
    RooArgList B_bar_pdfs;

    TString prefix(name_prefix);
    prefix += "_";

    auto efficiency_model = channel_config["efficiencyModel"].get<int>();

    // The possibility to supply more than one efficiency file is a remnant of
    // hybrid histo-KDE models. I'm reluctant to remove it as it might be useful
    // in the systematics estimation.
    std::vector<std::string> efficiency_files;
    if (efficiency_model > 4) {
        efficiency_files.push_back(channel_config["efficiencyFile"]);
    }

    AngularPDF* cr_pdf_B = new AngularPDF(
        (prefix + "cr_angular_pdf_B").Data(), (prefix + "cr_angular_pdf_B").Data(), false,
        efficiency_model, efficiency_files, *thetat_, *thetab_, *phit_, *ap_, *apa_, *a0_, *ata_);
    AngularPDF* cr_pdf_B_bar = new AngularPDF(
        (prefix + "cr_angular_pdf_B_bar").Data(), (prefix + "cr_angular_pdf_B_bar").Data(), true,
        efficiency_model, efficiency_files, *thetat_, *thetab_, *phit_, *ap_, *apa_, *a0_, *ata_);

    B_pdfs.add(*cr_pdf_B);
    B_bar_pdfs.add(*cr_pdf_B_bar);

    if (scf) {
        RooAbsPdf* scf_angular_pdf = CreateSCFPDF(name_prefix, channel_config);
        B_pdfs.add(*scf_angular_pdf);
        B_bar_pdfs.add(*scf_angular_pdf);
    }

    if (bkg) {
        RooAbsPdf* bkg_angular_pdf = CreateAngularSCFBKGPDF(name_prefix + "_bkg_");
        B_pdfs.add(*bkg_angular_pdf);
        B_bar_pdfs.add(*bkg_angular_pdf);
    }

    RooArgList fractions;
    if (scf && !bkg) {
        RooRealVar* cr_scf_f_ = new RooRealVar(prefix + "cr_scf_f", "f_{cr}", constants::fraction_cr_of_crscf, 0.80, 0.99);
        cr_scf_f_->setConstant();
        fractions.add(*cr_scf_f_);
    } else if (scf && bkg) {
        RooRealVar* cr_f_ = new RooRealVar(prefix + "cr_f", "f_{cr}", constants::fraction_cr_of_crscfbkg, 0.10, 0.99);
        RooRealVar* scf_f_ = new RooRealVar(prefix + "scf_f", "f_{scf}", constants::fraction_scf_of_crscfbkg, 0.10, 0.99);
        cr_f_->setConstant();
        scf_f_->setConstant();
        fractions.add(*cr_f_);
        fractions.add(*scf_f_);
    }

    RooAddPdf* pdf_B =
        new RooAddPdf(prefix + "angular_pdf_B", prefix + "angular_pdf_B", B_pdfs, fractions);
    RooAddPdf* pdf_B_bar = new RooAddPdf(prefix + "angular_pdf_B_bar", prefix + "angular_pdf_B_bar",
                                         B_bar_pdfs, fractions);

    RooSimultaneous* sim_pdf =
        new RooSimultaneous(prefix + "angular_pdf", prefix + "angular_pdf", *decaytype_);
    sim_pdf->addPdf(*pdf_B, "FB");
    sim_pdf->addPdf(*pdf_B_bar, "FA");
    sim_pdf->addPdf(*pdf_B, "SB");
    sim_pdf->addPdf(*pdf_B_bar, "SA");

    return sim_pdf;
}

/**
 * Create a time-dependent PDF with the requested components
 *
 * @param scf Whether to add self-crossfeed component
 * @param bkg Whether to add background component
 * @return RooSimultaneous* The complete PDF
 */
RooSimultaneous* FitterCPV::CreateTimeDependentPDF(const std::string channel_name,
                                                   const nlohmann::json common_config,
                                                   const nlohmann::json channel_config) {
    RooArgList FB_pdfs;
    RooArgList FA_pdfs;
    RooArgList SB_pdfs;
    RooArgList SA_pdfs;

    bool scf = false;
    bool bkg = false;
    if (common_config["components"] == "CRSCF") {
        scf = true;
    } else if (common_config["components"] == "all") {
        scf = true;
        bkg = true;
    }

    DtCPPDF* cr_pdf_FB = 0;
    DtCPPDF* cr_pdf_FA = 0;
    DtCPPDF* cr_pdf_SB = 0;
    DtCPPDF* cr_pdf_SA = 0;
    CreateDtCPPDFs(cr_pdf_FB, cr_pdf_FA, cr_pdf_SB, cr_pdf_SA, channel_name, common_config,
                   channel_config);
    FB_pdfs.add(*cr_pdf_FB);
    FA_pdfs.add(*cr_pdf_FA);
    SB_pdfs.add(*cr_pdf_SB);
    SA_pdfs.add(*cr_pdf_SA);

    RooProdPdf* scf_pdf_FB = 0;
    RooProdPdf* scf_pdf_FA = 0;
    RooProdPdf* scf_pdf_SB = 0;
    RooProdPdf* scf_pdf_SA = 0;
    if (scf) {
        CreateFunctionalDtSCFBKGPDFs(scf_pdf_FB, scf_pdf_FA, scf_pdf_SB, scf_pdf_SA, channel_name,
                                     channel_config, true);
        FB_pdfs.add(*scf_pdf_FB);
        FA_pdfs.add(*scf_pdf_FA);
        SB_pdfs.add(*scf_pdf_SB);
        SA_pdfs.add(*scf_pdf_SA);
    }

    RooProdPdf* bkg_pdf_FB = 0;
    RooProdPdf* bkg_pdf_FA = 0;
    RooProdPdf* bkg_pdf_SB = 0;
    RooProdPdf* bkg_pdf_SA = 0;
    if (bkg) {
        CreateFunctionalDtSCFBKGPDFs(bkg_pdf_FB, bkg_pdf_FA, bkg_pdf_SB, bkg_pdf_SA, channel_name,
                                     channel_config, false);
        FB_pdfs.add(*bkg_pdf_FB);
        FA_pdfs.add(*bkg_pdf_FA);
        SB_pdfs.add(*bkg_pdf_SB);
        SA_pdfs.add(*bkg_pdf_SA);
    }

    TString prefix = channel_name.c_str();
    prefix += "_";

    RooArgList fractions;
    if (scf && !bkg) {
        RooRealVar* cr_scf_f_ = new RooRealVar(prefix + "cr_scf_f", "f_{cr}", constants::fraction_cr_of_crscf, 0.80, 0.99);
        cr_scf_f_->setConstant();
        fractions.add(*cr_scf_f_);
    } else if (scf && bkg) {
        RooRealVar* cr_f_ = new RooRealVar(prefix + "cr_f", "f_{cr}", constants::fraction_cr_of_crscfbkg, 0.10, 0.99);
        RooRealVar* scf_f_ = new RooRealVar(prefix + "scf_f", "f_{scf}", constants::fraction_scf_of_crscfbkg, 0.10, 0.99);
        cr_f_->setConstant();
        scf_f_->setConstant();
        fractions.add(*cr_f_);
        fractions.add(*scf_f_);
    }

    RooAddPdf* pdf_FB = new RooAddPdf(prefix + "pdf_FB", prefix + "pdf_FB", FB_pdfs, fractions);
    RooAddPdf* pdf_FA = new RooAddPdf(prefix + "pdf_FA", prefix + "pdf_FA", FA_pdfs, fractions);
    RooAddPdf* pdf_SB = new RooAddPdf(prefix + "pdf_SB", prefix + "pdf_SB", SB_pdfs, fractions);
    RooAddPdf* pdf_SA = new RooAddPdf(prefix + "pdf_SA", prefix + "pdf_SA", SA_pdfs, fractions);

    RooSimultaneous* sim_pdf = new RooSimultaneous(prefix + "sim_pdf", prefix + "sim_pdf", *decaytype_);
    sim_pdf->addPdf(*pdf_FB, "FB");
    sim_pdf->addPdf(*pdf_FA, "FA");
    sim_pdf->addPdf(*pdf_SB, "SB");
    sim_pdf->addPdf(*pdf_SA, "SA");

    return sim_pdf;
}

/**
 * Fit specified PDF, show the results, and create plots if config demands it
 *
 * @param timedep Whether a time-dependent fit should be carried out
 * @param scf Add self-crossfeed component to the fit
 * @param bkg Add background component to the fit
 */
void FitterCPV::Fit(const nlohmann::json config) {
    Log::print(Log::info, "Fitting %i events.\n", data_->numEntries());

    Log::LogLine(Log::info) << "Running a fit with the following components: "
                            << config["components"];

    tau_->setConstant(true);
    dm_->setConstant(true);

    result_ = pdf_->fitTo(*data_, RooFit::ConditionalObservables(conditional_vars_argset_),
                          RooFit::Hesse(false), RooFit::Minos(false), RooFit::Minimizer("Minuit2"),
                          RooFit::Save(true), RooFit::NumCPU(config["numCPUs"]));

    if (result_) {
        result_->Print();
    }
}

/**
 * Create all relevant plots
 */
void FitterCPV::CreatePlots(const nlohmann::json config) const {
    for (auto& chan : config["channels"].items()) {
        const std::string channel_name = chan.key().c_str();
        auto channel_config = chan.value();

        auto common_config = config;
        common_config.erase("channels");

        PlotChannel(common_config, channel_config, channel_name);
    }
}

/**
 * Plot the PDF and data including all components plotted separately
 *
 * This function is rather convoluted to get around a RooFit bug which causes
 * problems when trying to plot complicated PDFs with conditional observables;
 * see
 * https://root-forum.cern.ch/t/rooaddpdf-and-integration-over-conditional-observables/27891
 * Instead of normally plotting the PDFs, we generate Asimov datasets and plot those.
 *
 * @param common_config Part of the config common to all channels
 * @param channel_config Part of the config specific for this channel
 * @param channel_name Name of this channel
 */
void FitterCPV::PlotChannel(const nlohmann::json common_config, const nlohmann::json channel_config,
                            const std::string channel_name) const {
    Log::LogLine(Log::info) << "Creating plots for channel " << channel_name;

    std::vector<RooAbsPdf*> components;
    RooAbsPdf* all_histpdf = CreateHistPdf(common_config, channel_name, components);

    std::string channel_cut = "channel_cat==channel_cat::" + channel_name;
    RooDataSet* channel_data = dynamic_cast<RooDataSet*>(data_->reduce(channel_cut.c_str()));

    std::vector<RooCmdArg> plot_options = {RooFit::Slice(*channel_cat_, channel_name.c_str())};
    RooArgSet* observables = pdf_->getObservables(data_);
    for (auto observable : tools::ToVector<RooRealVar*>(*observables)) {
        tools::PlotWithPull(*observable, conditional_vars_argset_, *channel_data, *all_histpdf, result_, components,
                     common_config["numCPUs"], channel_name, "", plot_options);
    }

    thetat_->setBins(40);
    thetab_->setBins(40);
    phit_->setBins(40);

    RooDataHist hist("hist", "hist", *observables, *data_);
    tools::PlotVars2D(*thetat_, *thetab_, hist, channel_name);
    tools::PlotVars2D(*thetat_, *phit_, hist, channel_name);
    tools::PlotVars2D(*thetab_, *phit_, hist, channel_name);

    // We change the binning for 2D plots, so we have to generate new binned
    // dataset, to avoid dealing with rebinning.

    // We need this to get the exact number of events to be generated for
    // the PDF distribution. This is unnecessary for the above, because we
    // create a RooHistPdf from it and that is normalized to the data when
    // plotted.
    RooDataSet* dataset_B = dynamic_cast<RooDataSet*>(
        channel_data->reduce("decaytype==decaytype::FB||decaytype==decaytype::SB"));
    RooDataSet* dataset_A = dynamic_cast<RooDataSet*>(
        channel_data->reduce("decaytype==decaytype::FA||decaytype==decaytype::SA"));

    TString chan_name = channel_name;
    RooAddPdf* pdf_FB = dynamic_cast<RooAddPdf*>(pdf_->getPdf("{" + chan_name + ";FB}"));
    RooAddPdf* pdf_FA = dynamic_cast<RooAddPdf*>(pdf_->getPdf("{" + chan_name + ";FA}"));
    RooDataHist* all_hist =
        pdf_FB->generateBinned(*observables, dataset_B->sumEntries(), RooFit::ExpectedData(true));
    RooDataHist* all_hist_FA =
        pdf_FA->generateBinned(*observables, dataset_A->sumEntries(), RooFit::ExpectedData(true));
    all_hist->add(*all_hist_FA);

    // // Set the current directory back to the one for plots (ugly ROOT stuff)
    // if (output_file_) {
    //     output_file_->cd();
    // }
    tools::PlotPull2D(*thetat_, *thetab_, hist, *all_hist, channel_name);
    tools::PlotPull2D(*thetat_, *phit_, hist, *all_hist, channel_name);
    tools::PlotPull2D(*thetab_, *phit_, hist, *all_hist, channel_name);
    // const double chi2 = Calculate3DChi2(hist, *all_hist);
    // std::cout << "Chi2 = " << chi2 << std::endl;

    // // SaveChi2Scan(sim_pdf, ap_, margin_apa);
    // // SaveChi2Scan(sim_pdf, apa_, margin_apa);
    // // SaveChi2Scan(sim_pdf, a0_, margin_a0);
    // // SaveChi2Scan(sim_pdf, ata_, margin_ata);

    delete dataset_B;
    delete dataset_A;
}

RooAbsPdf* FitterCPV::CreateHistPdf(const nlohmann::json common_config,
                                    const std::string channel_name,
                                    std::vector<RooAbsPdf*>& components) const {
    // PDF for B_bar differs, so we have to generate them separately. (One
    // could also use *extended* RooSimultaneous or 'ProtoData' for
    // RooSimultaneous.)
    RooAbsPdf* cr_pdf_B;
    RooAbsPdf* cr_pdf_B_bar;
    RooAbsPdf* scf_pdf;
    RooAbsPdf* bkg_pdf;

    TString chan_name = channel_name;
    RooAddPdf* pdf_FB = dynamic_cast<RooAddPdf*>(pdf_->getPdf("{" + chan_name + ";FB}"));
    RooAddPdf* pdf_FA = dynamic_cast<RooAddPdf*>(pdf_->getPdf("{" + chan_name + ";FA}"));
    RooAddPdf* pdf_SB = dynamic_cast<RooAddPdf*>(pdf_->getPdf("{" + chan_name + ";SB}"));
    RooAddPdf* pdf_SA = dynamic_cast<RooAddPdf*>(pdf_->getPdf("{" + chan_name + ";SA}"));
    // pdf_FB->Print();
    // pdf_FA->Print();
    // pdf_SB->Print();
    // pdf_SA->Print();

    // TODO: Make this work for TD PDFs as well
    cr_pdf_B = dynamic_cast<RooAbsPdf*>(pdf_FB->pdfList().find(chan_name + "_cr_angular_pdf_B"));
    assert(cr_pdf_B != nullptr);
    // cr_pdf_B->Print();

    cr_pdf_B_bar = dynamic_cast<RooAbsPdf*>(pdf_FA->pdfList().find(chan_name + "_cr_angular_pdf_B_bar"));
    assert(cr_pdf_B_bar != nullptr);
    // cr_pdf_B_bar->Print();

    bool scf = false;
    bool bkg = false;
    if (common_config["components"] == "CRSCF") {
        scf = true;
    } else if (common_config["components"] == "all") {
        scf = true;
        bkg = true;
    }

    RooArgSet* observables = pdf_->getObservables(data_);

    RooDataHist* cr_hist = cr_pdf_B->generateBinned(*observables, 1000, RooFit::ExpectedData(true));
    RooDataHist* cr_hist_B_bar =
        cr_pdf_B_bar->generateBinned(*observables, 1000, RooFit::ExpectedData(true));

    // Add histos from both particle and anti-particle PDFs to create the
    // final RooHistPdf.
    cr_hist->add(*cr_hist_B_bar);
    RooHistPdf* cr_histpdf = new RooHistPdf("cr_histpdf", "cr_histpdf", *observables, *cr_hist);

    RooHistPdf* scf_histpdf;
    if (scf) {
        const int scf_pdf_FB_index = pdf_FB->pdfList().index(chan_name + "_scf_pdf");
        scf_pdf = dynamic_cast<RooAbsPdf*>(pdf_FB->pdfList().at(scf_pdf_FB_index));
        RooDataHist* scf_hist =
            scf_pdf->generateBinned(*observables, 1000, RooFit::ExpectedData(true));
        scf_histpdf = new RooHistPdf("scf_histpdf", "scf_histpdf", *observables, *scf_hist);
    }

    RooHistPdf* bkg_histpdf;
    if (bkg) {
        const int bkg_pdf_FB_index = pdf_FB->pdfList().index(chan_name + "_bkg_pdf");
        bkg_pdf = dynamic_cast<RooAbsPdf*>(pdf_FB->pdfList().at(bkg_pdf_FB_index));
        RooDataHist* bkg_hist =
            bkg_pdf->generateBinned(*observables, 1000, RooFit::ExpectedData(true));
        bkg_histpdf = new RooHistPdf("bkg_histpdf", "bkg_histpdf", *observables, *bkg_hist);
    }

    RooAbsArg* cr_scf_f = pdf_->getVariables()->find(chan_name + "_cr_scf_f");
    RooAbsArg* cr_f = pdf_->getVariables()->find(chan_name + "_cr_f");
    RooAbsArg* scf_f = pdf_->getVariables()->find(chan_name + "_scf_f");
    RooAbsPdf* all_histpdf;
    if ((!scf) && (!bkg)) {
        all_histpdf = cr_histpdf;
    } else if (scf && (!bkg)) {
        all_histpdf = new RooAddPdf("crscf_histpdf", "crscf_histpdf",
                                    RooArgList(*cr_histpdf, *scf_histpdf), *cr_scf_f);
    } else if (scf && bkg) {
        all_histpdf = new RooAddPdf("all_histpdf", "all_histpdf",
                                    RooArgList(*cr_histpdf, *scf_histpdf, *bkg_histpdf),
                                    RooArgList(*cr_f, *scf_f));
    }

    // // Set the current directory back to the one for plots (ugly ROOT stuff)
    // if (output_file_) {
    //     output_file_->cd();
    // }

    // // const double margin_ap = 0.005;
    // // const double margin_apa = 0.25;
    // // const double margin_a0 = 0.0015;
    // // const double margin_ata = 0.3;
    // // SaveLikelihoodScan(sim_pdf, ap_, margin_ap);
    // // SaveLikelihoodScan(sim_pdf, apa_, margin_apa);
    // // SaveLikelihoodScan(sim_pdf, a0_, margin_a0);
    // // SaveLikelihoodScan(sim_pdf, ata_, margin_ata);
    // // SaveLikelihoodScan(sim_pdf, ap_, apa_, margin_ap, margin_apa);
    // // SaveLikelihoodScan(sim_pdf, ap_, a0_, margin_ap, margin_a0);
    // // SaveLikelihoodScan(sim_pdf, ap_, ata_, margin_ap, margin_ata);
    // // SaveLikelihoodScan(sim_pdf, apa_, a0_, margin_apa, margin_a0);
    // // SaveLikelihoodScan(sim_pdf, apa_, ata_, margin_ap, margin_ata);
    // // SaveLikelihoodScan(sim_pdf, a0_, ata_, margin_a0, margin_ata);

    components.push_back(cr_histpdf);
    if (scf) {
        components.push_back(scf_histpdf);
    }
    if (bkg) {
        components.push_back(bkg_histpdf);
    }

    return all_histpdf;
}

/**
 * Read in a JSON file, parse it, and return a json object.
 *
 * @param filename Path to the file to be read
 */
/* static */ nlohmann::json FitterCPV::ReadJSONConfig(const char* filename) {
    std::ifstream filestream(filename);
    if (!filestream.good()) {
        Log::print(Log::error, "Specified config file '%s' doesn't exist!\n", filename);
        exit(5);
    }
    nlohmann::json config;
    filestream >> config;
    return config;
}

/**
 * Take a JSON config with fit ranges and set the ranges of all variables of
 * the fitter corresponding to dataset variables according to the config
 *
 * @param config JSON config with the fit ranges
 */
void FitterCPV::ChangeFitRanges(const nlohmann::json& config) {
    // The following line is C++17 or newer and we work around it to be able to
    // compile on Travis which uses ancient gcc.
    // for (auto& [key, value] : config.items()) {
    for (auto& element : config.items()) {
        const char* name = element.key().c_str();
        nlohmann::json value = element.value();

        if (dataset_vars_argset_.find(name)) {
            if (value.contains("min")) {
                Log::print(Log::debug, "Changing var %s min to %f\n", name,
                           value["min"].get<float>());
                dynamic_cast<RooRealVar&>(dataset_vars_argset_[name]).setMin(value["min"]);
            }

            if (value.contains("max")) {
                Log::print(Log::debug, "Changing var %s max to %f\n", name,
                           value["max"].get<float>());
                dynamic_cast<RooRealVar&>(dataset_vars_argset_[name]).setMax(value["max"]);
            }
        } else {
            Log::print(Log::warning, "Variable '%s' not found\n", name);
            exit(9);
        }
    }
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

    DtCPPDF mixing_pdf_FB("mixing_pdf_FB", "mixing_pdf_FB", false, true, perfect_tagging_,
                          efficiency_model_, efficiency_files_, *thetat_, *thetab_, *phit_, ap, apa,
                          a0, ata, xp, x0, xt, yp, y0, yt, *tagwtag_, *dt_, *tau_, *dm_, *expmc_,
                          *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_, *vrzerr_, *vrchi2_,
                          *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF mixing_pdf_FA("mixing_pdf_FA", "mixing_pdf_FA", true, true, perfect_tagging_,
                          efficiency_model_, efficiency_files_, *thetat_, *thetab_, *phit_, ap, apa,
                          a0, ata, xpb, x0b, xtb, ypb, y0b, ytb, *tagwtag_, *dt_, *tau_, *dm_,
                          *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_, *vrzerr_,
                          *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF mixing_pdf_SB("mixing_pdf_SB", "mixing_pdf_SB", false, false, perfect_tagging_,
                          efficiency_model_, efficiency_files_, *thetat_, *thetab_, *phit_, ap, apa,
                          a0, ata, xp, x0, xt, yp, y0, yt, *tagwtag_, *dt_, *tau_, *dm_, *expmc_,
                          *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_, *vrzerr_, *vrchi2_,
                          *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

    DtCPPDF mixing_pdf_SA("mixing_pdf_SA", "mixing_pdf_SA", true, false, perfect_tagging_,
                          efficiency_model_, efficiency_files_, *thetat_, *thetab_, *phit_, ap, apa,
                          a0, ata, xpb, x0b, xtb, ypb, y0b, ytb, *tagwtag_, *dt_, *tau_, *dm_,
                          *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_, *vrzerr_,
                          *vrchi2_, *vrndf_, *vtntrk_, *vtzerr_, *vtchi2_, *vtndf_, *vtistagl_);

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
    RooAbsPdf::GenSpec* genSpec_FB =
        mixing_pdf_FB.prepareMultiGen(varsToGenerate, RooFit::NumEvents(num_fav));
    Log::print(Log::info, "1/4 GenSpecs ready\n");
    RooAbsPdf::GenSpec* genSpec_FA =
        mixing_pdf_FA.prepareMultiGen(varsToGenerate, RooFit::NumEvents(num_fav));
    Log::print(Log::info, "2/4 GenSpecs ready\n");
    RooAbsPdf::GenSpec* genSpec_SB =
        mixing_pdf_SB.prepareMultiGen(varsToGenerate, RooFit::NumEvents(num_sup));
    Log::print(Log::info, "3/4 GenSpecs ready\n");
    RooAbsPdf::GenSpec* genSpec_SA =
        mixing_pdf_SA.prepareMultiGen(varsToGenerate, RooFit::NumEvents(num_sup));
    Log::print(Log::info, "4/4 GenSpecs ready\n");

    RooAbsData::setDefaultStorageType(RooAbsData::Tree);
    RooDataSet* dataset;
    TString filename;

    RooRealVar btagmcli("btagmcli", "btagmcli", -511, 511);
    RooRealVar brecflav("brecflav", "brecflav", -1, 1);

    for (int i = 1; i <= num_toys; i++) {
        //  RooDataSet* dataset_FB = mixing_pdf_FB.generate(varsToGenerate,
        //  RooFit::NumEvents(num_fav));
        RooDataSet* dataset_FB = mixing_pdf_FB.generate(*genSpec_FB);
        decaytype_->setLabel("FB");
        dataset_FB->addColumn(*decaytype_);
        brecflav.setVal(1);
        btagmcli.setVal(-511);
        dataset_FB->addColumn(brecflav);
        dataset_FB->addColumn(btagmcli);
        dataset = (RooDataSet*)dataset_FB->Clone();
        delete dataset_FB;
        Log::print(Log::info, "1/4 decay type ready.\n");

        //  RooDataSet* dataset_FA = mixing_pdf_FA.generate(varsToGenerate,
        //  RooFit::NumEvents(num_fav));
        RooDataSet* dataset_FA = mixing_pdf_FA.generate(*genSpec_FA);
        decaytype_->setLabel("FA");
        dataset_FA->addColumn(*decaytype_);
        brecflav.setVal(-1);
        btagmcli.setVal(511);
        dataset_FA->addColumn(brecflav);
        dataset_FA->addColumn(btagmcli);
        dataset->append(*dataset_FA);
        delete dataset_FA;
        Log::print(Log::info, "2/4 decay types ready.\n");

        //  RooDataSet* dataset_SB = mixing_pdf_SB.generate(varsToGenerate,
        //  RooFit::NumEvents(num_sup));
        RooDataSet* dataset_SB = mixing_pdf_SB.generate(*genSpec_SB);
        decaytype_->setLabel("SB");
        dataset_SB->addColumn(*decaytype_);
        brecflav.setVal(1);
        btagmcli.setVal(511);
        dataset_SB->addColumn(brecflav);
        dataset_SB->addColumn(btagmcli);
        dataset->append(*dataset_SB);
        delete dataset_SB;
        Log::print(Log::info, "3/4 decay types ready.\n");

        //  RooDataSet* dataset_SA = mixing_pdf_SA.generate(varsToGenerate,
        //  RooFit::NumEvents(num_sup));
        RooDataSet* dataset_SA = mixing_pdf_SA.generate(*genSpec_SA);
        decaytype_->setLabel("SA");
        dataset_SA->addColumn(*decaytype_);
        brecflav.setVal(-1);
        btagmcli.setVal(-511);
        dataset_SA->addColumn(brecflav);
        dataset_SA->addColumn(btagmcli);
        dataset->append(*dataset_SA);
        delete dataset_SA;
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

    //  result_ = mixing_pdf_FB.fitTo(*new_dataset,
    // RooFit::ConditionalObservables(conditional_vars_argset_),
    //          RooFit::Minimizer("Minuit2"), RooFit::Hesse(false), RooFit::Minos(false),
    //          RooFit::Save(true), RooFit::NumCPU(num_CPUs_));
}


/**
 * Set the directory to which to ouput plots.
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
    input =
        std::regex_replace(input, std::regex("xy"), "xp,x0,xt,yp,y0,yt,xpb,x0b,xtb,ypb,y0b,ytb");
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
void FitterCPV::SetOutputFile(TFile* file) { output_file_ = file; }

/**
 * Create a string that holds initial values, fit results and errors.
 */
const std::string FitterCPV::CreateResultsString() const {
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

    if (IsTimeDependent()) {
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
        if (!IsTimeDependent()) {
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
        if (!IsTimeDependent()) {
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
 * Save results, status codes, covariance matrix quality to the ROOT output file.
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

/**
 * Read efficiency files, create a model and return a histogram with binned
 * values sampled from the model. In principle, no files need to be supplied
 * for functional efficiency models.
 *
 * @param files Files from which a model is to be built
 * @param model Specification of which model is to be built; see the Efficiency class
 *
 * @return TH3D* Binned efficiency model
 */
TH3D* FitterCPV::GetBinnedEfficiency(std::vector<std::string> files, const int model) {
    Efficiency eff;
    for (auto file : files) {
        eff.ReadInFile(file.c_str());
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

/**
 * Create 2D plots of the current efficiency model.
 */
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

/**
 * Save -2 log(likelihood) scan of a single variable.
 *
 * @param pdf PDF to be used for likelihood calculation
 * @param var Variable to scan
 * @param margin [optional] Set plot range [var_value - margin, var_value + margin]
 */
// const void FitterCPV::SaveLikelihoodScan(RooAbsPdf& pdf, RooRealVar* var, const double margin) {
//     TString name;
//     name = "nll_";
//     name += var->GetName();

//     TCanvas canvas(name, name, 500, 500);

//     const int steps = 100;
//     const double orig_val = var->getVal();
//     const double min = margin ? orig_val - margin : var->getMin();
//     const double max = margin ? orig_val + margin : var->getMax();
//     const double stepsize = (max - min) / steps;

//     TH1D h1_nll("h1_" + name, "h1_" + name, steps, min, max);
//     RooAbsReal* nll;
//     nll = pdf.createNLL(*dataset_, RooFit::NumCPU(num_CPUs_));

//     for (Int_t i = 0; i < steps; i++) {
//         var->setVal(i * stepsize + min + stepsize / 2);
//         Log::print(Log::debug, "Computing %i/%i likelihood function.\n", i + 1, steps);
//         const double nll_val = 2 * nll->getVal();
//         if (!std::isnan(nll_val)) {
//             h1_nll.Fill(var->getVal(), nll_val);
//         }
//     }

//     // Set optimal y-axis range, while ignoring empty bins
//     double min_bin_content = h1_nll.GetMaximum();
//     for (int i = 1; i < h1_nll.GetSize(); i++) {
//         if (min_bin_content > h1_nll.GetBinContent(i) && h1_nll.GetBinContent(i) != 0) {
//             min_bin_content = h1_nll.GetBinContent(i);
//         }
//     }
//     const double ymargin = (h1_nll.GetMaximum() - min_bin_content) * 0.05;
//     h1_nll.GetYaxis()->SetRangeUser(min_bin_content - ymargin, h1_nll.GetMaximum() + ymargin);

//     delete nll;

//     // Set the current directory back to the one for plots (ugly ROOT stuff)
//     if (output_file_) {
//         output_file_->cd();
//     }

//     h1_nll.GetXaxis()->SetTitle(var->GetTitle());
//     h1_nll.GetYaxis()->SetTitle("");
//     h1_nll.SetTitle("");
//     h1_nll.SetStats(kFALSE);
//     h1_nll.Draw("HIST");
//     h1_nll.Write();
//     canvas.SaveAs(constants::format);
//     var->setVal(orig_val);
// }

/**
 * Save -2 log(likelihood) scan of a two variables.
 *
 * @param pdf PDF to be used for likelihood calculation
 * @param var1 x-axis variable to scan
 * @param var2 y-axis variable to scan
 * @param margin1 [optional] Set plot range [var1_value - margin1, var1_value + margin1]
 * @param margin2 [optional] Set plot range [var2_value - margin2, var2_value + margin2]
 */
// const void FitterCPV::SaveLikelihoodScan(RooAbsPdf& pdf, RooRealVar* var1, RooRealVar* var2,
//                                          const double margin1, const double margin2) {
//     TString name;
//     name = "nll2_";
//     name += var1->GetName();
//     name += "_";
//     name += var2->GetName();

//     TCanvas canvas(name, name, 500, 500);

//     const int steps1 = 30;
//     const double orig_val1 = var1->getVal();
//     const double min1 = margin1 ? orig_val1 - margin1 : var1->getMin();
//     const double max1 = margin1 ? orig_val1 + margin1 : var1->getMax();
//     const double stepsize1 = (max1 - min1) / steps1;

//     const int steps2 = 30;
//     const double orig_val2 = var2->getVal();
//     const double min2 = margin2 ? orig_val2 - margin2 : var2->getMin();
//     const double max2 = margin2 ? orig_val2 + margin2 : var2->getMax();
//     const double stepsize2 = (max2 - min2) / steps2;

//     TH2F h2_nll("h2_" + name, "h2_" + name, steps1, min1, max1, steps2, min2, max2);
//     RooAbsReal* nll;
//     nll = pdf.createNLL(*dataset_, RooFit::NumCPU(num_CPUs_));

//     for (Int_t i = 0; i < steps1; i++) {
//         for (Int_t j = 0; j < steps2; j++) {
//             var1->setVal(i * stepsize1 + min1 + stepsize1 / 2);
//             var2->setVal(j * stepsize2 + min2 + stepsize2 / 2);
//             Log::print(Log::debug, "Computing %i/%i likelihood function.\n", i * steps2 + j + 1,
//                        steps1 * steps2);
//             h2_nll.Fill(var1->getVal(), var2->getVal(), 2 * nll->getVal());
//         }
//     }

//     delete nll;

//     // Set the current directory back to the one for plots (ugly ROOT stuff)
//     if (output_file_) {
//         output_file_->cd();
//     }

//     canvas.SetRightMargin(0.14);
//     h2_nll.GetXaxis()->SetTitle(var1->GetTitle());
//     h2_nll.GetYaxis()->SetTitle(var2->GetTitle());
//     h2_nll.GetZaxis()->SetTitle("");
//     h2_nll.SetTitle("");
//     h2_nll.SetStats(kFALSE);
//     h2_nll.SetOption("colz");
//     h2_nll.Draw();
//     h2_nll.Write();
//     canvas.SaveAs(constants::format);
//     var1->setVal(orig_val1);
//     var2->setVal(orig_val2);
// }

/**
 * Calculate a chi2 value from a histogram and a PDF, assuming Poisson
 * distribution of data in the bins.
 *
 * @param data
 * @param pdf
 * @return const double
 */
const double FitterCPV::Calculate3DChi2(const RooDataHist& data, const RooDataHist& pdf) {
    // TODO: Remove the commented histogram lines
    // TH1I h_bin_content("h_bin_content", "Bin Content", 100, 0, 99);
    double chi2 = 0;
    int bins_used = 0;
    for (int i = 0; i < data.numEntries(); i++) {
        data.get(i);
        pdf.get(i);
        double data_bin = data.weight();
        double pdf_bin = pdf.weight();
        // h_bin_content.Fill(data_bin);
        if (data_bin > 5) {
            chi2 += (data_bin - pdf_bin) * (data_bin - pdf_bin) / pdf_bin;
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
    Log::print(Log::info, "Bins used for chi2: %i/%i (%.1f%%)\n", bins_used, data.numEntries(),
               (double)bins_used / data.numEntries() * 100);

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
        data_->reduce("decaytype==decaytype::FB||decaytype==decaytype::SB"));
    RooDataSet* dataset_B_bar = dynamic_cast<RooDataSet*>(
        data_->reduce("decaytype==decaytype::FA||decaytype==decaytype::SA"));

    RooDataHist hist("hist", "hist", RooArgSet(*thetat_, *thetab_, *phit_), *data_);

    RooAbsPdf* pdf_B = pdf.getPdf("FB");
    RooAbsPdf* pdf_B_bar = pdf.getPdf("FA");

    TH1D h1_chi2("h1_" + name, "h1_" + name, steps, min, max);

    for (Int_t i = 0; i < steps; i++) {
        var->setVal(i * stepsize + min + stepsize / 2);
        Log::print(Log::debug, "Computing chi2 %i/%i.\n", i + 1, steps);

        RooDataHist* cr_hist =
            pdf_B->generateBinned(RooArgSet(*thetat_, *thetab_, *phit_), dataset_B->sumEntries(),
                                  RooFit::ExpectedData(true));
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

/**
 * Create a functional model of the background angular distribution.
 */
RooAbsPdf* FitterCPV::CreateAngularSCFBKGPDF(const std::string prefix) const {
    TString pfx = prefix;

    // Background phit model
    RooRealVar* phit_poly_p2 = new RooRealVar(pfx + "phit_poly_p2", "p_(2)", 0.040);
    RooRealVar* phit_f = new RooRealVar(pfx + "phit_f", "f_(poly)", 0.660);
    RooPolynomial* phit_poly =
        new RooPolynomial(pfx + "phit_poly", "phit_poly", *phit_, *phit_poly_p2, 2);
    RooRealVar* phit_offset = new RooRealVar(pfx + "phit_offset", "#phi_(t)^(offset)", 0.103);
    RooFormulaVar* phit_phit =
        new RooFormulaVar(pfx + "phit_phit", "phit_phit", "phit - " + pfx + "phit_offset",
                          RooArgList(*phit_, *phit_offset));
    RooGenericPdf* phit_cos = new RooGenericPdf(
        pfx + "phit_cos", "phit_cos", "cos(" + pfx + "phit_phit)^2", RooArgList(*phit_phit));
    RooAddPdf* phit_model = new RooAddPdf(pfx + "phit_model", "phit_model",
                                          RooArgList(*phit_poly, *phit_cos), RooArgList(*phit_f));

    // Background thetat model
    RooRealVar* thetat_p1 = new RooRealVar(pfx + "thetat_p1", "p_{1}", 0);
    RooRealVar* thetat_p2 = new RooRealVar(pfx + "thetat_p2", "p_{2}", -1.147);
    RooRealVar* thetat_p3 = new RooRealVar(pfx + "thetat_p3", "p_{3}", 0);
    RooRealVar* thetat_p4 = new RooRealVar(pfx + "thetat_p4", "p_{4}", 0.174);
    RooRealVar* thetat_p5 = new RooRealVar(pfx + "thetat_p5", "p_{5}", 0);
    RooRealVar* thetat_p6 = new RooRealVar(pfx + "thetat_p6", "p_{6}", -0.029);
    RooChebychev* thetat_model = new RooChebychev(
        pfx + "thetat_model", "thetat_model", *thetat_,
        RooArgList(*thetat_p1, *thetat_p2, *thetat_p3, *thetat_p4, *thetat_p5, *thetat_p6));

    // Background thetab model
    RooRealVar* thetab_gaus_mu = new RooRealVar(pfx + "thetab_gaus_mu", "#mu", 2.895);
    RooRealVar* thetab_gaus_sigma_l =
        new RooRealVar(pfx + "thetab_gaus_sigma_l", "#sigma_(L)", 0.902);
    RooRealVar* thetab_gaus_sigma_r =
        new RooRealVar(pfx + "thetab_gaus_sigma_r", "#sigma_(R)", 0.089);
    RooBifurGauss* thetab_gaus =
        new RooBifurGauss(pfx + "thetab_gaus", "thetab_gaus", *thetab_, *thetab_gaus_mu,
                          *thetab_gaus_sigma_l, *thetab_gaus_sigma_r);
    RooRealVar* thetab_exp_alpha = new RooRealVar(pfx + "thetab_exp_alpha", "#alpha", -2.182);
    RooExponential* thetab_exp =
        new RooExponential(pfx + "thetab_exp", "thetab_exp", *thetab_, *thetab_exp_alpha);
    RooRealVar* thetab_f = new RooRealVar(pfx + "thetab_f", "f_(exp)", 0.661);

    RooAddPdf* thetab_model =
        new RooAddPdf(pfx + "thetab_model", "thetab_model", RooArgList(*thetab_exp, *thetab_gaus),
                      RooArgList(*thetab_f));

    RooProdPdf* pdf =
        new RooProdPdf(pfx + "pdf", pfx + "pdf", RooArgList(*thetat_model, *thetab_model, *phit_model));

    return pdf;
}

/**
 * WARNING: Experimental approach, should not be used without tweaking the code. Set the
 * self-crossfeed model to histogram constructed from a KDE model from a supplied file.
 *
 * @param file File holding a Meerkat KDE model
 */
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

    // scf_angular_pdf_ = scf_angular_kde_;

    // OneDimPhaseSpace* phasespace_thetat = new OneDimPhaseSpace{"phasespace_thetat",
    // thetat_->getMin(), thetat_->getMax()}; OneDimPhaseSpace* phasespace_thetab = new
    // OneDimPhaseSpace{"phasespace_thetab", thetab_->getMin(), thetab_->getMax()};
    // OneDimPhaseSpace* phasespace_phit = new OneDimPhaseSpace{"phasespace_phit", phit_->getMin(),
    // phit_->getMax()}; CombinedPhaseSpace* phasespace = new CombinedPhaseSpace{"phasespace",
    // phasespace_thetat, phasespace_thetab,
    //                               phasespace_phit};
    // BinnedDensity* binned_scf_kde = new BinnedDensity("binned_scf_kde", phasespace, file);
    // RooArgList list(*thetat_, *thetab_, *phit_);
    // RooMeerkatPdf* meerkat_pdf =
    //     new RooMeerkatPdf("meerkat_pdf", "meerkat_pdf", list, binned_scf_kde);
    // scf_angular_pdf_ = meerkat_pdf;
}

/**
 * Get histogram efficiency from supplied file. The name of the histogram must be 'scf_hist_pdf'.
 *
 * @param file File with the efficiency histogram
 */
RooAbsPdf* FitterCPV::GetHistoSCF(const std::string filename) const {
    Log::LogLine(Log::info) << "Setting up SCF RooHistPdf model from file '" << filename << "'";
    TFile f(filename.c_str(), "READ");
    RooHistPdf* temp_pdf = dynamic_cast<RooHistPdf*>(f.Get("scf_hist_pdf"));

    // We create a RooHistPdfFast (our version of RooHistPdf that implements
    // caching of its "analytical integral") from the original RooHistPdf to
    // avoid the huge performance hit due to a RooFit bug.
    return new RooHistPdfFast("scf_pdf", "scf_pdf", RooArgSet(*thetat_, *thetab_, *phit_),
                              temp_pdf->dataHist());
}
/**
 * Determine whether angular coordinates are close to (any) edge of the phasespace.
 *
 * @param vals Coordinates to be checked
 * @param margin Distance (fraction of phasespace range in each axis) to be taken as 'close'
 *
 * @return int 0 for not close, 1 for close to lower edge, 2 for close to upper edge
 */
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

/**
 * Create the complete PDF based on a supplied config.
 */
RooSimultaneous* FitterCPV::CreatePDF(const nlohmann::json config) {
    std::map<std::string, RooAbsPdf*> pdf_map;
    for (auto& chan : config["channels"].items()) {
        std::string channel_name = chan.key();
        auto channel_config = chan.value();
        Log::LogLine(Log::debug) << "Creating PDF for channel " << channel_name;

        auto common_config = config;
        common_config.erase("channels");

        RooSimultaneous* channel_pdf =
            CreateChannelPDF(channel_name, channel_config, common_config);
        pdf_map[channel_name] = channel_pdf;
    }

    RooSimultaneous* pdf = new RooSimultaneous("pdf", "pdf", pdf_map, *channel_cat_);
    return pdf;
}

/**
 * Create the complete PDF for a single channel based on a supplied config.
 */
RooSimultaneous* FitterCPV::CreateChannelPDF(const std::string channel_name,
                                             const nlohmann::json channel_config,
                                             const nlohmann::json common_config) {
    bool scf = false;
    bool bkg = false;
    if (common_config["components"] == "CRSCF") {
        scf = true;
    } else if (common_config["components"] == "all") {
        scf = true;
        bkg = true;
    }
    RooSimultaneous* channel_pdf =
        common_config.contains("timeIndependent")
            ? CreateAngularPDF(channel_name, scf, bkg, channel_config)
            : CreateTimeDependentPDF(channel_name, common_config, channel_config);

    if (common_config.contains("modelParameters")) {
        tools::ChangeModelParameters(channel_pdf, common_config["modelParameters"], channel_name + "_");
    }
    if (channel_config.contains("modelParameters")) {
        tools::ChangeModelParameters(channel_pdf, channel_config["modelParameters"], channel_name + "_");
    }

    return channel_pdf;
}

/**
 * Create a dataset with all the categories needed for the complete fit.
 */
RooDataSet* FitterCPV::GetData(const nlohmann::json config) {
    std::vector<RooDataSet*> datasets;

    for (auto& chan : config["channels"].items()) {
        channel_cat_->defineType(chan.key().c_str());
    }

    for (auto& chan : config["channels"].items()) {
        const std::string channel_name = chan.key().c_str();
        auto channel_config = chan.value();
        Log::LogLine(Log::debug) << "Constructing dataset for channel " << channel_name;

        auto common_config = config;
        common_config.erase("channels");

        RooDataSet* channel_dataset = GetChannelData(channel_name, channel_config, common_config);
        channel_cat_->setLabel(channel_name.c_str());
        channel_dataset->addColumn(*channel_cat_);
        datasets.push_back(channel_dataset);
    }

    RooDataSet* dataset = nullptr;
    for (auto channel_dataset : datasets) {
        if (dataset == nullptr) {
            dataset = static_cast<RooDataSet*>(channel_dataset->Clone());
            dataset->SetName("combined_dataset");
            dataset->SetTitle("Combined Dataset");
            continue;
        }
        dataset->append(*channel_dataset);
    }

    vrzerr_ = static_cast<RooRealVar*>(dataset->addColumn(*vrzerr_formula_));
    vrzerr_->setRange(0, 10000);
    vtzerr_ = static_cast<RooRealVar*>(dataset->addColumn(*vtzerr_formula_));
    vtzerr_->setRange(0, 10000);

    conditional_vars_argset_.add(*vrzerr_);
    conditional_vars_argset_.add(*vtzerr_);
    conditional_vars_argset_.add(*decaytype_);

    // Bind the variables to the dataset, so that dataset->get(i) changes values of, e.g., expno_
    const RooArgSet* vars = dataset->get();
    for (RooRealVar** var : dataset_vars_) {
        *var = static_cast<RooRealVar*>(vars->find((*var)->GetName()));
    }

    Log::LogLine(Log::info) << "Total dataset events: " << dataset->numEntries();

    return dataset;
}

RooDataSet* FitterCPV::GetChannelData(const std::string channel_name,
                                      const nlohmann::json channel_config,
                                      const nlohmann::json common_config) {
    TChain* chain = new TChain("h2000");
    int num_files = 0;
    if (common_config["MC"].get<bool>() == false) {
        for (auto& file_name : channel_config["inputFiles"]["data"].items()) {
            if (!chain->Add(file_name.value().get<std::string>().c_str(), 0)) {
                Log::LogLine(Log::error) << "File " << file_name.value() << " not found!";
                exit(9);
            }
            num_files++;
        }
    } else {
        for (auto& file_name : channel_config["inputFiles"]["signalMC"].items()) {
            if (!chain->Add(file_name.value().get<std::string>().c_str(), 0)) {
                Log::LogLine(Log::error) << "File " << file_name.value() << " not found!";
                exit(9);
            }
            num_files++;
        }

        if (common_config["components"] == "all") {
            for (auto& file_name : channel_config["inputFiles"]["genericMC"].items()) {
                if (!chain->Add(file_name.value().get<std::string>().c_str(), 0)) {
                    Log::LogLine(Log::error) << "File " << file_name.value() << " not found!";
                    exit(9);
                }
                num_files++;
            }
        }
    }

    Log::print(Log::info, "Reading %i input files for channel %s\n", num_files, channel_name.c_str());

    TTree* tree = nullptr;
    if (common_config.contains("events")) {
        double fraction_events = common_config["events"].get<double>() / chain->GetEntries();
        std::string fraction_string = "LocalEntry$<LocalEntries$*";
        fraction_string += std::to_string(fraction_events);
        Log::print(
            Log::info,
            "Taking a total of %i events (%.1f%%) properly distributed across all input files\n",
            common_config["events"].get<int>(), fraction_events * 100);
        tree = chain->CopyTree(fraction_string.c_str());
    }

    TString common_cuts = tools::GetCommonCutsString();

    if (common_config["components"] == "CR") {
        common_cuts += "&&evmcflag==1";
    }

    TString FB_cuts;
    TString FA_cuts;
    TString SB_cuts;
    TString SA_cuts;
    // TODO: Treat case that key exists but is false
    if (common_config.contains("perfectTagging")) {
        FB_cuts = "brecflav==1&&btagmcli<0";
        FA_cuts = "brecflav==-1&&btagmcli>=0";
        SB_cuts = "brecflav==1&&btagmcli>=0";
        SA_cuts = "brecflav==-1&&btagmcli<0";
    } else {
        // TODO: Treat case that key exists but is false
        if (common_config.contains("generatorLevel")) {
            Log::print(Log::warning,
                       "Attempting to use realistic tagging with generator level fit. This will "
                       "probably end badly. Consider using the '--perfect-tag' switch.\n");
        }
        FB_cuts = "brecflav==1&&tagqr<0";
        FA_cuts = "brecflav==-1&&tagqr>=0";
        SB_cuts = "brecflav==1&&tagqr>=0";
        SA_cuts = "brecflav==-1&&tagqr<0";
    }

    // A temporary RooDataSet is created from the whole tree and then we apply cuts to get
    // the 4 different B and f flavor datasets, as that is faster then reading the tree 4 times
    std::string dataset_name = channel_name;
    dataset_name += "_dataset";
    RooDataSet* temp_dataset =
        new RooDataSet(dataset_name.c_str(), dataset_name.c_str(), tree == nullptr ? chain : tree,
                       dataset_vars_argset_, common_cuts);
    Log::print(Log::debug, "Num events passing common cuts: %i\n", temp_dataset->numEntries());

    // We add an identifying label to each of the 4 categories and then combine it into a single
    // dataset for RooSimultaneous fitting
    RooDataSet* dataset_FB = static_cast<RooDataSet*>(temp_dataset->reduce(FB_cuts));
    decaytype_->setLabel("FB");
    dataset_FB->addColumn(*decaytype_);

    RooDataSet* dataset_FA = static_cast<RooDataSet*>(temp_dataset->reduce(FA_cuts));
    decaytype_->setLabel("FA");
    dataset_FA->addColumn(*decaytype_);

    RooDataSet* dataset_SB = static_cast<RooDataSet*>(temp_dataset->reduce(SB_cuts));
    decaytype_->setLabel("SB");
    dataset_SB->addColumn(*decaytype_);

    RooDataSet* dataset_SA = static_cast<RooDataSet*>(temp_dataset->reduce(SA_cuts));
    decaytype_->setLabel("SA");
    dataset_SA->addColumn(*decaytype_);

    delete temp_dataset;

    RooDataSet* dataset = static_cast<RooDataSet*>(dataset_FB->Clone());
    delete dataset_FB;
    dataset->append(*dataset_FA);
    delete dataset_FA;
    dataset->append(*dataset_SB);
    delete dataset_SB;
    dataset->append(*dataset_SA);
    delete dataset_SA;

    Log::print(Log::info, "Num events passing all initial cuts: %i\n", dataset->numEntries());

    return dataset;
}

RooAbsPdf* FitterCPV::CreateSCFPDF(const std::string channel_name,
                                   const nlohmann::json channel_config) const {
    if (channel_config.contains("scfHisto")) {
        RooAbsPdf* histo_pdf = GetHistoSCF(channel_config["scfHisto"]);
        std::string name = channel_name;
        name += "_";
        name += histo_pdf->GetName();
        histo_pdf->SetName(name.c_str());
        return histo_pdf;
    } else {
        std::string name = channel_name;
        name += "_scf_";
        return CreateAngularSCFBKGPDF(name);
    }
}

/**
 * Check whether the current fitter is time-dependent or not
 *
 * @return true Time-dependent
 * @return false Not time-dependent
 */
bool FitterCPV::IsTimeDependent() const {
    return pdf_->getObservables(data_)->contains(*dt_);
}

/**
 * Take JSON config with initial parameter values and output ordered std::array
 *
 * @param initial_pars JSON config with the initial parameter values
 * @return std::array<double, 16> Ordered array readable by InitVars()
 */
std::array<double, 16> FitterCPV::ToParInputArray(nlohmann::json initial_pars) {
    const char* names[16] = {"ap", "apa", "a0",  "ata", "xp",  "x0",  "xt",  "yp",
                           "y0", "yt",  "xpb", "x0b", "xtb", "ypb", "y0b", "ytb"};
    std::array<double, 16> pars;
    for (int i = 0; i < 16; i++) {
        if (!initial_pars.contains(names[i])) {
            Log::LogLine(Log::error) << "'" << names[i] << "' not present in initialPars!";
        }
        pars[i] = initial_pars[names[i]].get<double>();
    }
}
