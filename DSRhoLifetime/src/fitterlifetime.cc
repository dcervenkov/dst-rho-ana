/**
 *  @file    fitterlifetime.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2016-04-05
 *
 *  @brief This class performs the lifetime fitting itself as well as plotting
 *
 */

#include "fitterlifetime.h"

// BASF includes
#include "tatami/tatami.h"

// ROOT includes
#include "TEnv.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooExtendPdf.h"
#include "RooHistPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"

// Local includes
#include "constants.h"
#include "dtpdf.h"
#include "tools.h"

FitterLifetime::FitterLifetime() {
	vrusable_ = new RooRealVar( "vrusable", "vrusable", 0, 1 );
	vrvtxz_ = new RooRealVar( "vrvtxz", "vrvtxz", -10, 10 );
	vrerr6_ = new RooRealVar( "vrerr6", "vrerr6", -1, 1 );
	vrchi2_ = new RooRealVar( "vrchi2", "vrchi2", 0, 10000000 );
	vreffxi_ = new RooRealVar( "vreffxi", "vreffxi", 0, 10000000 );
	vrndf_ = new RooRealVar( "vrndf", "vrndf", 0, 100 );
	vreffndf_ = new RooRealVar( "vreffndf", "vreffndf", 0, 100 );
	vrntrk_ = new RooRealVar( "vrntrk", "vrntrk", 0, 100 );

	vtusable_ = new RooRealVar( "vtusable", "vtusable", 1 );
	vtvtxz_ = new RooRealVar( "vtvtxz", "vtvtxz", -10, 10 );
	vtchi2_ = new RooRealVar( "vtchi2", "vtchi2", 1.6 );
	vtndf_ = new RooRealVar( "vtndf", "vtndf", 4 );
	vterr6_ = new RooRealVar( "vterr6", "vterr6", -1, 1 );
	vtchi2_ = new RooRealVar( "vtchi2", "vtchi2", 0, 10000000 );
	vtndf_ = new RooRealVar( "vtndf", "vtndf", 0, 100 );
	vtntrk_ = new RooRealVar( "vtntrk", "vtntrk", 0, 100 );
	vtistagl_ = new RooRealVar( "vtistagl", "vtistagl", 0, 100 );

	expno_ = new RooRealVar("expno", "expno", 0, 100);
	expmc_ = new RooRealVar( "expmc", "expmc", 0, 10 );

	evmcflag_ = new RooRealVar( "evmcflag", "evmcflag", 0, 10 );
	brecflav_ = new RooRealVar( "brecflav", "brecflav", -1, 1 );
	btagmcli_ = new RooRealVar( "btagmcli", "btagmcli", -1000, 1000 );
	tagqr_ = new RooRealVar( "tagqr", "tagqr", -1, 1 );
	tagwtag_ = new RooRealVar( "tagwtag", "tagwtag", 0, 0.5 );

	benergy_ = new RooRealVar( "benergy", "benergy", 0, 100 );
	mbc_ = new RooRealVar( "mbc", "mbc", 5, 6 );

	shcosthb_ = new RooRealVar( "shcosthb", "shcosthb", -1, 1 );

	// TODO: beta*gamma should be computed not a constant and c should be taken from constants.cc
	dt_formula_ = new RooFormulaVar( "dt", "#Deltat [ps]", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)", RooArgSet(*vrvtxz_, *vtvtxz_));
	vrzerr_formula_ = new RooFormulaVar( "vrzerr", "#sigma z_{rec} [cm]", "sqrt(vrerr6)", RooArgSet(*vrerr6_));
	vtzerr_formula_ = new RooFormulaVar( "vtzerr", "#sigma z_{tag} [cm]", "sqrt(vterr6)", RooArgSet(*vterr6_));

	tau_ = new RooRealVar ("tau", "#tau", constants::tau - 1, constants::tau + 1);
	dm_ = new RooRealVar ("dm", "#Deltam", constants::dm - 1, constants::dm + 1);

	decaytype_ = new RooCategory("decaytype", "decaytype");
	decaytype_->defineType("a",1);
	decaytype_->defineType("ab",2);
	decaytype_->defineType("b",3);
	decaytype_->defineType("bb",4);

	// All the variables present in the ntuple are added to a vector so we can later
	// iterate through them for convenience
	vars_.push_back(&expno_);
	vars_.push_back(&expmc_);

	vars_.push_back(&evmcflag_);
	vars_.push_back(&brecflav_);
	vars_.push_back(&btagmcli_);
	vars_.push_back(&tagqr_);
	vars_.push_back(&tagwtag_);

	vars_.push_back(&benergy_);
	vars_.push_back(&mbc_);

	vars_.push_back(&shcosthb_);

	vars_.push_back(&vrusable_);
	vars_.push_back(&vrvtxz_);
	vars_.push_back(&vrerr6_);
	vars_.push_back(&vrchi2_);
	vars_.push_back(&vreffxi_);
	vars_.push_back(&vrndf_);
	vars_.push_back(&vreffndf_);
	vars_.push_back(&vrntrk_);

	vars_.push_back(&vtusable_);
	vars_.push_back(&vtvtxz_);
	vars_.push_back(&vterr6_);
	vars_.push_back(&vtchi2_);
	vars_.push_back(&vtndf_);
	vars_.push_back(&vtntrk_);

	vars_.push_back(&vtistagl_);

	// The variables present in the ntuple are added to an RooArgSet that will be needed
	// when we create a RooDataSet from the input_tree
	for(auto var : vars_) {
		argset_.add(**var);
	}

	num_CPUs_ = 1;
	do_lifetime_fit_ = false;
	do_mixing_fit_ = false;
	make_plots_ = false;
	perfect_tagging_ = false;

	dt_ = NULL;
	vrzerr_ = NULL;
	vtzerr_ = NULL;
	vars = NULL;

}

FitterLifetime::~FitterLifetime() {
	if(output_file_) {
		if(output_file_->IsOpen()) {
			output_file_->Close();
		}
	}
}

void FitterLifetime::PlotVar(RooRealVar& var) const {
	TCanvas* canvas = new TCanvas(var.GetName(), var.GetTitle(), 500, 500);

	double range_min;
	double range_max;
	dataset_->getRange(var, range_min, range_max, 0.05);

	RooPlot* plot = var.frame(range_min, range_max);

	dataset_->plotOn(plot);
	plot->Draw();

	plot->SetTitle("");
	plot->GetXaxis()->SetTitle(TString(var.GetTitle()));
	plot->GetYaxis()->SetTitle("");

	canvas->Write();
	canvas->SaveAs(constants::format);
}

void FitterLifetime::Test() {
	DtPDF lifetime_pdf("lifetime_pdf", "lifetime_pdf",
			*dt_,
			*tau_,
			*expmc_,
			*expno_,
			*shcosthb_,
			*benergy_,
			*mbc_,
			*vrntrk_,
			*vrzerr_,
			*vrchi2_,
			*vrndf_,
			*vtntrk_,
			*vtzerr_,
			*vtchi2_,
			*vtndf_,
			*vtistagl_);

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

	DtPDF mixing_pdf_a("mixing_pdf_a", "mixing_pdf_a", true, perfect_tagging_,
			S,
			A,
			*tagwtag_,
			*dt_,
			*tau_,
			*dm_,
			*expmc_,
			*expno_,
			*shcosthb_,
			*benergy_,
			*mbc_,
			*vrntrk_,
			*vrzerr_,
			*vrchi2_,
			*vrndf_,
			*vtntrk_,
			*vtzerr_,
			*vtchi2_,
			*vtndf_,
			*vtistagl_);

	DtPDF mixing_pdf_ab("mixing_pdf_ab", "mixing_pdf_ab", true, perfect_tagging_,
			bS,
			bA,
			*tagwtag_,
			*dt_,
			*tau_,
			*dm_,
			*expmc_,
			*expno_,
			*shcosthb_,
			*benergy_,
			*mbc_,
			*vrntrk_,
			*vrzerr_,
			*vrchi2_,
			*vrndf_,
			*vtntrk_,
			*vtzerr_,
			*vtchi2_,
			*vtndf_,
			*vtistagl_);

	DtPDF mixing_pdf_b("mixing_pdf_b", "mixing_pdf_b", false, perfect_tagging_,
			S,
			A,
			*tagwtag_,
			*dt_,
			*tau_,
			*dm_,
			*expmc_,
			*expno_,
			*shcosthb_,
			*benergy_,
			*mbc_,
			*vrntrk_,
			*vrzerr_,
			*vrchi2_,
			*vrndf_,
			*vtntrk_,
			*vtzerr_,
			*vtchi2_,
			*vtndf_,
			*vtistagl_);

	DtPDF mixing_pdf_bb("mixing_pdf_bb", "mixing_pdf_bb", false, perfect_tagging_,
			bS,
			bA,
			*tagwtag_,
			*dt_,
			*tau_,
			*dm_,
			*expmc_,
			*expno_,
			*shcosthb_,
			*benergy_,
			*mbc_,
			*vrntrk_,
			*vrzerr_,
			*vrchi2_,
			*vrndf_,
			*vtntrk_,
			*vtzerr_,
			*vtchi2_,
			*vtndf_,
			*vtistagl_);

	RooSimultaneous sim_pdf("sim_pdf", "sim_pdf", *decaytype_);
	sim_pdf.addPdf(mixing_pdf_a, "a");
	sim_pdf.addPdf(mixing_pdf_ab, "ab");
	sim_pdf.addPdf(mixing_pdf_b, "b");
	sim_pdf.addPdf(mixing_pdf_bb, "bb");

	dt_->setRange("dtFitRange", constants::cuts::dt_low, constants::cuts::dt_high);

//	tau_->setConstant(true);
//	dm_->setConstant(true);

	if (do_lifetime_fit_) {
		result_ = lifetime_pdf.fitTo(*dataset_, RooFit::ConditionalObservables(argset_), RooFit::Minimizer("Minuit2"),
				RooFit::Range("dtFitRange"), RooFit::Save(true), RooFit::NumCPU(num_CPUs_));

		if (make_plots_) {
			PlotWithPull(*dt_, *dataset_, lifetime_pdf);
		}
	}

	if (do_mixing_fit_) {
	result_ = sim_pdf.fitTo(*dataset_, RooFit::ConditionalObservables(argset_), RooFit::Minimizer("Minuit2"),
			RooFit::Range("dtFitRange"), RooFit::Save(true), RooFit::NumCPU(num_CPUs_));

		if (make_plots_) {
			RooDataSet* dataset_a = static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::a"));
			PlotWithPull(*dt_, *dataset_a, mixing_pdf_a);

			RooDataSet* dataset_b = static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::b"));
			PlotWithPull(*dt_, *dataset_b, mixing_pdf_b);

			RooDataSet* dataset_ab = static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::ab"));
			PlotWithPull(*dt_, *dataset_ab, mixing_pdf_ab);

			RooDataSet* dataset_bb = static_cast<RooDataSet*>(dataset_->reduce("decaytype==decaytype::bb"));
			PlotWithPull(*dt_, *dataset_bb, mixing_pdf_bb);
		}
	}

}

void FitterLifetime::PlotWithPull(const RooRealVar& var, const RooAbsData& data, const RooAbsPdf& pdf, const char* title) const {
	TCanvas canvas(pdf.GetName(), var.GetTitle(), 500, 500);
	RooPlot* plot = var.frame(RooFit::Range("dtFitRange"));
	TPad *pad_var;
	TPad *pad_pull;
	pad_var = new TPad("pad_var", "pad_var", 0, 0.25, 1, 1);
	pad_pull = new TPad("pad_pull", "pad_pull", 0, 0, 1, 0.25);
	pad_var->Draw();
	pad_pull->Draw();

	pad_var->cd();
	pad_var->SetBottomMargin(0.0);
	pad_var->SetLeftMargin(0.12);

	data.plotOn(plot);
	// TODO: Arguments should not be hardcoded
	// Can't use more CPUs because a bug (?) in ROOT causes wrong normalization
	pdf.plotOn(plot,RooFit::ProjWData(argset_, data, kFALSE), RooFit::NumCPU(num_CPUs_),
			RooFit::NormRange("dtFitRange"), RooFit::Normalization(1.0/data.numEntries()));
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
	RooPlot* plot_pull_ = var.frame(RooFit::Title("Pull Distribution"), RooFit::Range("dtFitRange"));
	plot_pull_->SetTitle("");
	RooHist* hpull = plot->pullHist();
	hpull->SetFillColor(kGray);
	// The only working way to get rid of error bars; HIST draw option doesn't work with RooPlot
	for (int i = 0; i < hpull->GetN(); i++) {
		hpull->SetPointError(i, 0,0,0,0);
	}
	plot_pull_->addPlotable(hpull, "B");
	// We plot again without bars, so the points are not half covered by bars as in case of "BP" draw option.
	// We need to create and plot a clone, because the ROOT object ownership is transfered to the RooPlot by addPlotable().
	// If we just added the hpull twice, we would get a segfault.
	RooHist* hpull_clone = dynamic_cast<RooHist*>(hpull->Clone());
	plot_pull_->addPlotable(hpull_clone, "P"); // We plot again without bars, so the points are not half covered by bars

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

TPaveText* FitterLifetime::CreateStatBox(const double chi2, const bool position_top, const bool position_left) const {
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

	TPaveText *stat_box = new TPaveText(x_left, y_bottom, x_right, y_top, "NDC");
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

void FitterLifetime::ReadInFile(const char* file_path, const int& num_events) {
	TFile* input_file = new TFile(file_path);
	TTree* input_tree = dynamic_cast<TTree*>(input_file->Get("h2000"));

	if (num_events) {
		TTree* temp_tree = input_tree;
		input_tree = temp_tree->CloneTree(num_events);
		delete temp_tree;
	}

	TString common_cuts = tools::GetCommonCutsString();
	common_cuts += "&&evmcflag==1";

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
	RooDataSet* temp_dataset = new RooDataSet("dataset", "dataset", input_tree, argset_, common_cuts);

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

	// TODO: Shouldn't dt_ be added to argset_?
	vrzerr_ = static_cast<RooRealVar*>(dataset_->addColumn(*vrzerr_formula_));
	vtzerr_ = static_cast<RooRealVar*>(dataset_->addColumn(*vtzerr_formula_));
	dt_ = static_cast<RooRealVar*>(dataset_->addColumn(*dt_formula_));

	argset_.add(*vrzerr_);
	argset_.add(*vtzerr_);
	argset_.add(*decaytype_);

	// Bind the variables to the dataset, so that dataset->get(i) changes values of, e.g., expno_
	vars = dataset_->get();
	for (RooRealVar** var : vars_) {
		*var = static_cast<RooRealVar*>(vars->find((*var)->GetName()));
	}
}

void FitterLifetime::SetOutputDir(const char* output_dir) {
	gEnv->SetValue("Canvas.PrintDirectory", output_dir);
	printf("print dir: %s\n",gEnv->GetValue("Canvas.PrintDirectory","not found"));
	if (make_plots_) {
		output_file_ = new TFile(TString(output_dir) + "/plots.root", "RECREATE");
	}
}

