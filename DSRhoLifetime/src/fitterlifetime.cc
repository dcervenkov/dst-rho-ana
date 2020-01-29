/**
 *  @file    fitterlifetime.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2016-04-05
 *
 *  @brief This class performs the lifetime fitting itself as well as plotting
 *
 */

#include "fitterlifetime.h"

// Standard includes
#include <fstream>

// BASF includes
#include "tatami/tatami.h"

// ROOT includes
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooHist.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooVoigtian.h"
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
#include "dtbkg.h"
#include "dtpdf.h"
#include "log.h"
#include "nlohmann/json.hpp"
#include "tools.h"

FitterLifetime::FitterLifetime(const nlohmann::json config) {
    config_ = config;

    dt_ = new RooRealVar("dt", "#Deltat [ps]", constants::cuts::dt_low, constants::cuts::dt_high);
    thetat_ = new RooRealVar("thetat", "#theta_{t} [rad]", constants::cuts::thetat_low, constants::cuts::thetat_high);
    thetab_ = new RooRealVar("thetab", "#theta_{b} [rad]", constants::cuts::thetab_low,
                             constants::cuts::thetab_high);
    phit_ = new RooRealVar("phit", "#phi_{t} [rad]", constants::cuts::phit_low, constants::cuts::phit_high);

    vrusable_ = new RooRealVar("vrusable", "vrusable", 0, 1);
    vrvtxz_ = new RooRealVar("vrvtxz", "vrvtxz", -10, 10);
    vrerr6_ = new RooRealVar("vrerr6", "vrerr6", -1, 1);
    vrchi2_ = new RooRealVar("vrchi2", "vrchi2", 0, 10000000);
    vreffxi_ = new RooRealVar("vreffxi", "vreffxi", 0, 10000000);
    vrndf_ = new RooRealVar("vrndf", "vrndf", 0, 100);
    vreffndf_ = new RooRealVar("vreffndf", "vreffndf", 0, 100);
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

    evmcflag_ = new RooRealVar("evmcflag", "evmcflag", -1, 10);
    brecflav_ = new RooRealVar("brecflav", "brecflav", -1, 1);
    btagmcli_ = new RooRealVar("btagmcli", "btagmcli", -1000, 1000);
    tagqr_ = new RooRealVar("tagqr", "tagqr", -1, 1);
    tagwtag_ = new RooRealVar("tagwtag", "tagwtag", 0, 0.5);

    benergy_ = new RooRealVar("benergy", "benergy", 0, 100);
    mbc_ = new RooRealVar("mbc", "mbc", 5, 6);
    de_ = new RooRealVar("de", "de", -1, 1);
    csbdtg_ = new RooRealVar("csbdtg", "csbdtg", -1, 1);

    shcosthb_ = new RooRealVar("shcosthb", "shcosthb", -1, 1);

    tau_ = new RooRealVar("tau", "#tau", constants::tau - 1, constants::tau + 1);
    dm_ = new RooRealVar("dm", "#Deltam", constants::dm - 1, constants::dm + 1);

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
    conditional_vars_.push_back(&vreffxi_);
    conditional_vars_.push_back(&vrndf_);
    conditional_vars_.push_back(&vreffndf_);
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
    dataset_vars_.push_back(&phit_);
    dataset_vars_.push_back(&thetab_);
    dataset_vars_.push_back(&thetat_);

    // The variables present in the ntuple are added to an RooArgSet that will be needed
    // when we create a RooDataSet from the input_tree
    for (auto var : conditional_vars_) {
        conditional_argset_.add(**var);
    }

    for (auto var : dataset_vars_) {
        dataset_argset_.add(**var);
    }

    num_CPUs_ = 1;
    do_lifetime_fit_ = false;
    do_mixing_fit_ = false;
    make_plots_ = false;
    perfect_tagging_ = false;

    vars = nullptr;
}

FitterLifetime::~FitterLifetime() {
    if (output_file_) {
        if (output_file_->IsOpen()) {
            output_file_->Close();
        }
    }
}

void FitterLifetime::Test() {
    // tau_->setConstant(true);
    //	dm_->setConstant(true);

    if (do_lifetime_fit_) {
        Log::print(Log::info, "Fitting %i events\n", dataset_->numEntries());
        std::vector<RooAbsPdf*> components;
        RooAbsPdf* lifetime_pdf = CreateLifetimePDF(components, scf_, bkg_, use_physical_pdf_);
        result_ = lifetime_pdf->fitTo(*dataset_, RooFit::ConditionalObservables(conditional_argset_),
                                     RooFit::Minimizer("Minuit2"),
                                     RooFit::Save(true), RooFit::NumCPU(num_CPUs_));
        if (make_plots_) {
            tools::PlotWithPull(*dt_, conditional_argset_, *dataset_, *lifetime_pdf, result_, components);
        }
    }

    if (do_mixing_fit_) {
        // Small pars (r = 0.01)
        //	RooRealVar S("S", "S", 0.0144908);
        //	RooRealVar A("A", "A", -0.999801);

        // Large pars (r = 0.1)
        RooRealVar S("S", "S", 0.154719, -1, +1);
        RooRealVar A("A", "A", -0.980200, -1, +1);
        RooRealVar bS("bS", "bS", -0.0912793, -1, +1);
        RooRealVar bA("bA", "bA", -0.980200, -1, +1);

        S.setConstant(true);
        A.setConstant(true);
        bS.setConstant(true);
        bA.setConstant(true);

        DtPDF* cr_mixing_pdf_FB = new DtPDF("cr_mixing_pdf_FB", "cr_mixing_pdf_FB", true, perfect_tagging_, S, A, *tagwtag_,
                        *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
                        *vrerr6_, *vrchi2_, *vrndf_, *vtntrk_, *vterr6_, *vtchi2_, *vtndf_,
                        *vtistagl_);

        DtPDF* cr_mixing_pdf_FA = new DtPDF("cr_mixing_pdf_FA", "cr_mixing_pdf_FA", true, perfect_tagging_, bS, bA, *tagwtag_,
                            *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
                            *vrerr6_, *vrchi2_, *vrndf_, *vtntrk_, *vterr6_, *vtchi2_, *vtndf_,
                            *vtistagl_);

        DtPDF* cr_mixing_pdf_SB = new DtPDF("cr_mixing_pdf_SB", "cr_mixing_pdf_SB", false, perfect_tagging_, S, A, *tagwtag_,
                        *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_, *mbc_, *vrntrk_,
                        *vrerr6_, *vrchi2_, *vrndf_, *vtntrk_, *vterr6_, *vtchi2_, *vtndf_,
                        *vtistagl_);

        DtPDF* cr_mixing_pdf_SA = new DtPDF("cr_mixing_pdf_SA", "cr_mixing_pdf_SA", false, perfect_tagging_, bS, bA,
                            *tagwtag_, *dt_, *tau_, *dm_, *expmc_, *expno_, *shcosthb_, *benergy_,
                            *mbc_, *vrntrk_, *vrerr6_, *vrchi2_, *vrndf_, *vtntrk_, *vterr6_, *vtchi2_,
                            *vtndf_, *vtistagl_);

        RooArgList mixing_pdfs_FB;
        RooArgList mixing_pdfs_FA;
        RooArgList mixing_pdfs_SB;
        RooArgList mixing_pdfs_SA;
        mixing_pdfs_FB.add(*cr_mixing_pdf_FB);
        mixing_pdfs_FA.add(*cr_mixing_pdf_FA);
        mixing_pdfs_SB.add(*cr_mixing_pdf_SB);
        mixing_pdfs_SA.add(*cr_mixing_pdf_SA);

        if (scf_) {
            RooAddPdf* scf_dt_pdf_F = CreateVoigtGaussDtPdf("scf_dt_cf");
            RooAddPdf* scf_dt_pdf_S = CreateVoigtGaussDtPdf("scf_dt_dcs");
            mixing_pdfs_FB.add(*scf_dt_pdf_F);
            mixing_pdfs_FA.add(*scf_dt_pdf_F);
            mixing_pdfs_SB.add(*scf_dt_pdf_S);
            mixing_pdfs_SA.add(*scf_dt_pdf_S);
        }

        if (bkg_) {
            RooAddPdf* bkg_dt_pdf_F = CreateVoigtGaussDtPdf("bkg_dt_cf");
            RooAddPdf* bkg_dt_pdf_S = CreateVoigtGaussDtPdf("bkg_dt_dcs");
            mixing_pdfs_FB.add(*bkg_dt_pdf_F);
            mixing_pdfs_FA.add(*bkg_dt_pdf_F);
            mixing_pdfs_SB.add(*bkg_dt_pdf_S);
            mixing_pdfs_SA.add(*bkg_dt_pdf_S);
        }

        TString prefix = "";
        RooArgList fractions;
        if (scf_ && !bkg_) {
            RooRealVar* cr_scf_f_ = new RooRealVar(prefix + "cr_scf_f", "f_{cr}", constants::fraction_cr_of_crscf, 0.80, 0.99);
            cr_scf_f_->setConstant();
            fractions.add(*cr_scf_f_);
        } else if (scf_ && bkg_) {
            RooRealVar* cr_f_ = new RooRealVar(prefix + "cr_f", "f_{cr}", constants::fraction_cr_of_crscfbkg, 0.10, 0.99);
            RooRealVar* scf_f_ = new RooRealVar(prefix + "scf_f", "f_{scf}", constants::fraction_scf_of_crscfbkg, 0.10, 0.99);
            cr_f_->setConstant();
            scf_f_->setConstant();
            fractions.add(*cr_f_);
            fractions.add(*scf_f_);
        }

        RooAddPdf* mixing_pdf_FB = new RooAddPdf(prefix + "mixing_pdf_FB", prefix + "mixing_pdf_FB", mixing_pdfs_FB, fractions);
        RooAddPdf* mixing_pdf_FA = new RooAddPdf(prefix + "mixing_pdf_FA", prefix + "mixing_pdf_FA", mixing_pdfs_FA, fractions);
        RooAddPdf* mixing_pdf_SB = new RooAddPdf(prefix + "mixing_pdf_SB", prefix + "mixing_pdf_SB", mixing_pdfs_SB, fractions);
        RooAddPdf* mixing_pdf_SA = new RooAddPdf(prefix + "mixing_pdf_SA", prefix + "mixing_pdf_SA", mixing_pdfs_SA, fractions);

        RooSimultaneous sim_pdf("sim_pdf", "sim_pdf", *decaytype_);
        sim_pdf.addPdf(*mixing_pdf_FB, "FB");
        sim_pdf.addPdf(*mixing_pdf_FA, "FA");
        sim_pdf.addPdf(*mixing_pdf_SB, "SB");
        sim_pdf.addPdf(*mixing_pdf_SA, "SA");

        if (config_.contains("modelParameters")) {
            Log::LogLine(Log::debug) << "Global parameters:";
            tools::ChangeModelParameters(&sim_pdf, config_["modelParameters"]);
        }
        if (config_.contains("channels")) {
            Log::LogLine(Log::debug) << "Channel parameters:";
            tools::ChangeModelParameters(&sim_pdf, config_["channels"][channel_]["modelParameters"]);
        }

        Log::print(Log::info, "Fitting %i events\n", dataset_->numEntries());
        result_ = sim_pdf.fitTo(*dataset_, RooFit::ConditionalObservables(conditional_argset_),
                                RooFit::Minimizer("Minuit2"),
                                RooFit::Save(true), RooFit::NumCPU(num_CPUs_));

        if (make_plots_) {
            RooDataSet* dataset_FB =
                static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::FB"));
            tools::PlotWithPull(*dt_, conditional_argset_, *dataset_FB, *mixing_pdf_FB, result_,
                                tools::ToVector<RooAbsPdf*>(mixing_pdfs_FB));

            RooDataSet* dataset_FA =
                static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::FA"));
            tools::PlotWithPull(*dt_, conditional_argset_, *dataset_FA, *mixing_pdf_FA, result_,
                                tools::ToVector<RooAbsPdf*>(mixing_pdfs_FA));

            RooDataSet* dataset_SB =
                static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::SB"));
            tools::PlotWithPull(*dt_, conditional_argset_, *dataset_SB, *mixing_pdf_SB, result_,
                                tools::ToVector<RooAbsPdf*>(mixing_pdfs_SB));

            RooDataSet* dataset_SA =
                static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::SA"));
            tools::PlotWithPull(*dt_, conditional_argset_, *dataset_SA, *mixing_pdf_SA, result_,
                                tools::ToVector<RooAbsPdf*>(mixing_pdfs_SA));
        }
    }
}

void FitterLifetime::ReadInFile(const std::vector<const char*> file_names, const int& num_events) {
    TChain* chain = new TChain("h2000");
    int num_files = 0;
	for (auto& file_name : file_names) {
		if (!chain->Add(file_name, 0)) {
			Log::LogLine(Log::error) << "File " << file_name << " not found!";
			exit(9);
		}
		num_files++;
	}

    Log::print(Log::info, "Reading %i input files\n", num_files);

    TTree* tree = nullptr;
    if (num_events) {
        double fraction_events = (double) num_events / chain->GetEntries();
        std::string fraction_string = "LocalEntry$<LocalEntries$*";
        fraction_string += std::to_string(fraction_events);
        Log::print(
            Log::info,
            "Taking a total of %i events (%.1f%%) properly distributed across all input files\n",
            num_events, fraction_events * 100);
        tree = chain->CopyTree(fraction_string.c_str());
    }

    TString common_cuts = tools::GetCommonCutsString();
    if (scf_ == false && bkg_ == false) {
        Log::LogLine(Log::debug) << "Using evmcflag==1";
        common_cuts += "&&evmcflag==1";
    }

    TString FB_cuts;
    TString FA_cuts;
    TString SB_cuts;
    TString SA_cuts;
    if (perfect_tagging_) {
        FB_cuts = "brecflav==1&&btagmcli<0";
        FA_cuts = "brecflav==-1&&btagmcli>=0";
        SB_cuts = "brecflav==1&&btagmcli>=0";
        SA_cuts = "brecflav==-1&&btagmcli<0";
    } else {
        FB_cuts = "brecflav==1&&tagqr<0";
        FA_cuts = "brecflav==-1&&tagqr>=0";
        SB_cuts = "brecflav==1&&tagqr>=0";
        SA_cuts = "brecflav==-1&&tagqr<0";
    }

    // A temporary RooDataSet is created from the whole tree and then we apply cuts to get
    // the 4 different B and f flavor datasets, as that is faster then reading the tree 4 times
    RooDataSet* temp_dataset =
        new RooDataSet("dataset", "dataset", tree == nullptr ? chain : tree, dataset_argset_, common_cuts);

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

    dataset_ = static_cast<RooDataSet*>(dataset_FB->Clone());
    delete dataset_FB;
    dataset_->append(*dataset_FA);
    delete dataset_FA;
    dataset_->append(*dataset_SB);
    delete dataset_SB;
    dataset_->append(*dataset_SA);
    delete dataset_SA;

    dataset_argset_.add(*decaytype_);
    conditional_argset_.add(*decaytype_);

    // Bind the variables to the dataset, so that dataset->get(i) changes values of, e.g., expno_
    vars = dataset_->get();
    for (RooRealVar** var : conditional_vars_) {
        *var = static_cast<RooRealVar*>(vars->find((*var)->GetName()));
    }
}

/**
 * Create a Voigtian + Gaussian PDF.
 *
 * A common prefix is prepended to all object names.
 *
 * @param prefix Text to be prepended to the ROOT name
 *
 */
RooAddPdf* FitterLifetime::CreateVoigtGaussDtPdf(const std::string prefix) const {
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
 * Create a Physics conv. Resolution PDF.
 *
 * A common prefix is prepended to all object names.
 *
 * @param prefix Text to be prepended to the ROOT name
 *
 */
RooAbsPdf* FitterLifetime::CreatePhysicsBkgDtPdf(const std::string prefix) const {
    TString pre(prefix);

    RooRealVar* tau = new RooRealVar(pre + "_tau", "#tau_{bkg}", 1.5);
    RooRealVar* f_delta = new RooRealVar(pre + "_f_delta", "f_{d}", 0);
    RooRealVar* mu_delta = new RooRealVar(pre + "_mu_delta", "#mu_{d}", 0);
    RooRealVar* mu_lifetime = new RooRealVar(pre + "_mu_lifetime", "#mu_{l}", 0);
    RooRealVar* f_tail = new RooRealVar(pre + "_f_tail", "f_{t}", 0);
    RooRealVar* S_main = new RooRealVar(pre + "_S_main", "S_{m}", 1);
    RooRealVar* S_tail = new RooRealVar(pre + "_S_tail", "S_{t}", 1);

    DtBKG* model = new DtBKG(pre + "model", pre + "model", *dt_, *vrerr6_, *vterr6_, *tau,
                              *f_delta, *mu_delta, *mu_lifetime, *f_tail, *S_main, *S_tail);

    return model;
}

/**
 * Read a JSON config from a specified filename.
 */
nlohmann::json FitterLifetime::ReadInJSONFile(const char* filename) const {
    std::ifstream filestream(filename);
    if (!filestream.good()) {
        Log::print(Log::error, "Specified config file '%s' doesn't exist!\n", filename);
        exit(5);
    }
    Log::print(Log::debug, "Reading JSON config file '%s'\n", filename);

    nlohmann::json json;
    filestream >> json;
    return json;
}

RooAbsPdf* FitterLifetime::CreateLifetimePDF(std::vector<RooAbsPdf*>& components, const bool scf,
                                             const bool bkg, const bool physical_pdf) const {
    DtPDF* lifetime_cp_pdf = new DtPDF("lifetime_cp_pdf", "lifetime_cp_pdf", *dt_, *tau_, *expmc_, *expno_, *shcosthb_,
                       *benergy_, *mbc_, *vrntrk_, *vrerr6_, *vrchi2_, *vrndf_, *vtntrk_, *vterr6_,
                       *vtchi2_, *vtndf_, *vtistagl_);

    RooArgList lifetime_pdfs;
    lifetime_pdfs.add(*lifetime_cp_pdf);

    if (scf) {
        RooAbsPdf* lifetime_scf_pdf =
            (physical_pdf ? CreatePhysicsBkgDtPdf("scf_dt") : CreateVoigtGaussDtPdf("scf_dt"));
        lifetime_pdfs.add(*lifetime_scf_pdf);
    }

    if (bkg) {
        RooAbsPdf* lifetime_bkg_pdf =
            (physical_pdf ? CreatePhysicsBkgDtPdf("bkg_dt") : CreateVoigtGaussDtPdf("bkg_dt"));
        lifetime_pdfs.add(*lifetime_bkg_pdf);
    }

    TString prefix = "";
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

    RooAddPdf* lifetime_pdf = new RooAddPdf(prefix + "lifetime_pdf", prefix + "lifetime_pdf", lifetime_pdfs, fractions);

    if (config_.contains("modelParameters")) {
        Log::LogLine(Log::debug) << "Global parameters:";
        tools::ChangeModelParameters(lifetime_pdf, config_["modelParameters"]);
    }
    if (config_.contains("channels")) {
        Log::LogLine(Log::debug) << "Channel parameters:";
        tools::ChangeModelParameters(lifetime_pdf, config_["channels"][channel_]["modelParameters"]);
    }

    components = tools::ToVector<RooAbsPdf*>(lifetime_pdfs);
    return lifetime_pdf;
}

void FitterLifetime::SaveTXTResults(const char* results_file) const {
    std::stringstream buffer;
    if (do_lifetime_fit_ || do_mixing_fit_) {
        buffer << tau_->getVal() << ", " << tau_->getError(); 
    }
    if (do_mixing_fit_) {
        buffer << ", " << dm_->getVal() << ", " << dm_->getError(); 
    }
    buffer << std::endl;

    tools::SaveTextToFile(std::string(results_file), buffer.str());
}

void FitterLifetime::SetComponents(const std::string components) {
    if (components == "CRSCF") {
        scf_ = true;
    } else if (components == "all") {
        scf_ = true;
        bkg_ = true;
    }
}
