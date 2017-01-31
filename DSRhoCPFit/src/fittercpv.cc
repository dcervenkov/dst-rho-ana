/**
 *  @file    fittercpv.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2016-11-03
 *
 *  @brief Class that performs the CP fit itself as well as plotting
 *
 */

#include "fittercpv.h"

// Belle includes
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
#include "RooRandom.h"
#include "RooTreeDataStore.h"

// Local includes
#include "constants.h"
#include "dtcppdf.h"

FitterCPV::FitterCPV() {
	thetat_ = new RooRealVar( "thetat", "thetat", 0, constants::pi );
	thetab_ = new RooRealVar( "thetab", "thetab", 0, constants::pi );
	phit_ = new RooRealVar( "phit", "phit", -constants::pi, constants::pi );

	vrusable_ = new RooRealVar( "vrusable", "vrusable", 0, 1 );
	vrvtxz_ = new RooRealVar( "vrvtxz", "vrvtxz", -10, 10 );
	vrerr6_ = new RooRealVar( "vrerr6", "vrerr6", -1, 1 );
	vrchi2_ = new RooRealVar( "vrchi2", "vrchi2", 0, 10000000 );
//	vreffxi_ = new RooRealVar( "vreffxi", "vreffxi", 0, 10000000 );
	vrndf_ = new RooRealVar( "vrndf", "vrndf", 0, 100 );
//	vreffndf_ = new RooRealVar( "vreffndf", "vreffndf", 0, 100 );
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
//	conditional_vars_.push_back(&vreffxi_);
	conditional_vars_.push_back(&vrndf_);
//	conditional_vars_.push_back(&vreffndf_);
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


	// The variables present in the ntuple are added to an RooArgSet that will be needed
	// when we create a RooDataSet from the input_tree
	for(auto var : conditional_vars_) {
		conditional_vars_argset_.add(**var);
	}

	// TODO: Comments
	for(auto var : dataset_vars_) {
		dataset_vars_argset_.add(**var);
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
	if(output_file_) {
		if(output_file_->IsOpen()) {
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
void FitterCPV::Test() {
	// vars with phiweak = phiweak
//	double par_input[] = {
//			0.269,
//			0.56,
//			0.941,
//			3.11,
//
//			0.0816649,
//			0.0532961,
//			0.0829682,
//			0.0659988,
//			0.084614,
//			-0.0373802,
//
//			-0.102141,
//			-0.0845715,
//			-0.0587413,
//			-0.0243354,
//			-0.0533635,
//			0.0695015
//
//	};


	//// vars with phiweak = phiweak + pi, r = 0.10
	//double par_input[] = {
	//		0.269,
	//		0.56,
	//		0.941,
	//		3.11,

	//		0.0816649, // xp
	//		0.0532961, // x0
	//		0.0829682, // xt
	//		-0.0659988, // yp
	//		-0.084614,  // y0
	//		+0.0373802, // yt

	//		-0.102141,  // xpb
	//		-0.0845715, // x0b
	//		-0.0587413, // xtb
	//		+0.0243354, // ypb
	//		+0.0533635, // y0b
	//		-0.0695015  // ytb
	//};

	// vars with phiweak = phiweak + pi, r = 0.01
	double par_input[] = {
			0.269,
			0.56,
			0.941,
			3.11,

			0.00816649, // xp
			0.00532961, // x0
			0.00829682, // xt
			-0.00659988, // yp
			-0.0084614,  // y0
			+0.00373802, // yt

			-0.0102141,  // xpb
			-0.00845715, // x0b
			-0.00587413, // xtb
			+0.00243354, // ypb
			+0.00533635, // y0b
			-0.00695015  // ytb
	};

//	// vars with phiweak = phiweak + pi, x,y=0
//	double par_input[] = {
//			0.269,
//			0.56,
//			0.941,
//			3.11,
//
//			0, // xp
//			0, // x0
//			0, // xt
//			0, // yp
//			0, // y0
//			0, // yt
//
//			0, // xpb
//			0, // x0b
//			0, // xtb
//			0, // ypb
//			0, // y0b
//			0  // ytb
//	};

	RooRealVar ap ("ap","ap", par_input[0],0,0.5);
	RooRealVar apa("apa","apa", par_input[1],0,1);
	RooRealVar a0 ("a0","a0", par_input[2],0.8,1);
	RooRealVar a0a("a0a","a0a",0);
	RooFormulaVar at("at","sqrt(1-ap*ap-a0*a0)",RooArgSet(ap,a0));
	RooRealVar ata("ata","ata", par_input[3],2,4);

	RooRealVar xp ("xp","xp",par_input[4],-0.2,0.2);
	RooRealVar x0 ("x0","x0",par_input[5],-0.2,0.2);
	RooRealVar xt ("xt","xt",par_input[6],-0.2,0.2);

	RooRealVar yp ("yp","yp",par_input[7],-0.2,0.2);
	RooRealVar y0 ("y0","y0",par_input[8],-0.2,0.2);
	RooRealVar yt ("yt","yt",par_input[9],-0.2,0.2);

	RooRealVar xpb("xpb","xpb",par_input[10],-0.2,0.2);
	RooRealVar x0b("x0b","x0b",par_input[11],-0.2,0.2);
	RooRealVar xtb("xtb","xtb",par_input[12],-0.2,0.2);

	RooRealVar ypb("ypb","ypb",par_input[13],-0.2,0.2);
	RooRealVar y0b("y0b","y0b",par_input[14],-0.2,0.2);
	RooRealVar ytb("ytb","ytb",par_input[15],-0.2,0.2);

//	printf("DBG: Test\n");
//	printf("ap  = %f\n"
//		"apa = %f\n"
//		"a0  = %f\n"
//		"a0a = %f\n"
//		"at  = %f\n"
//		"ata = %f\n"
//		"xp  = %f\n"
//		"x0  = %f\n"
//		"xt  = %f\n"
//		"yp  = %f\n"
//		"y0  = %f\n"
//		"yt  = %f\n",
//		(double) ap.getVal(),
//		(double) apa.getVal(),
//		(double) a0.getVal(),
//		(double) a0a.getVal(),
//		(double) at.getVal(),
//		(double) ata.getVal(),
//		(double) xp.getVal(),
//		(double) x0.getVal(),
//		(double) xt.getVal(),
//		(double) yp.getVal(),
//		(double) y0.getVal(),
//		(double) yt.getVal());

	ap.setConstant();
	apa.setConstant();
	a0.setConstant();
	a0a.setConstant();
	ata.setConstant();

	//xp.setConstant();
	//x0.setConstant();
	//xt.setConstant();

	//yp.setConstant();
	//y0.setConstant();
	//yt.setConstant();

	//xpb.setConstant();
	//x0b.setConstant();
	//xtb.setConstant();

	//ypb.setConstant();
	//y0b.setConstant();
	//ytb.setConstant();


	DtCPPDF mixing_pdf_a("mixing_pdf_a", "mixing_pdf_a", true, perfect_tagging_,
			*thetat_,
			*thetab_,
			*phit_,

			ap,
			apa,
			a0,
			ata,
			xp,
			x0,
			xt,
			yp,
			y0,
			yt,

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

	DtCPPDF mixing_pdf_ab("mixing_pdf_ab", "mixing_pdf_ab", true, perfect_tagging_,
			*thetat_,
			*thetab_,
			*phit_,

			ap,
			apa,
			a0,
			ata,
			xpb,
			x0b,
			xtb,
			ypb,
			y0b,
			ytb,

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

	DtCPPDF mixing_pdf_b("mixing_pdf_b", "mixing_pdf_b", false, perfect_tagging_,
			*thetat_,
			*thetab_,
			*phit_,

			ap,
			apa,
			a0,
			ata,
			xp,
			x0,
			xt,
			yp,
			y0,
			yt,

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

	DtCPPDF mixing_pdf_bb("mixing_pdf_bb", "mixing_pdf_bb", false, perfect_tagging_,
			*thetat_,
			*thetab_,
			*phit_,

			ap,
			apa,
			a0,
			ata,
			xpb,
			x0b,
			xtb,
			ypb,
			y0b,
			ytb,

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

	dt_->setRange("dtFitRange", -15, 15);

	tau_->setConstant(true);
	dm_->setConstant(true);

	if (do_mixing_fit_) {
		result_ = sim_pdf.fitTo(*dataset_, RooFit::ConditionalObservables(conditional_vars_argset_), RooFit::Minimizer("Minuit2"),
			RooFit::Range("dtFitRange"), RooFit::Save(true), RooFit::NumCPU(num_CPUs_));
		result_->Print();

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

void FitterCPV::GenerateToys(const int num_events, const int num_toys) {
	printf("INFO: Generating dataset with %i events...\n", num_events);

	// vars with phiweak = phiweak + pi, r = 0.01
	double par_input[] = {
			0.269,
			0.56,
			0.941,
			3.11,

			0.00816649, // xp
			0.00532961, // x0
			0.00829682, // xt
			-0.00659988, // yp
			-0.0084614,  // y0
			+0.00373802, // yt

			-0.0102141,  // xpb
			-0.00845715, // x0b
			-0.00587413, // xtb
			+0.00243354, // ypb
			+0.00533635, // y0b
			-0.00695015  // ytb
	};

//	// vars with phiweak = phiweak + pi, r = 0.10
//	double par_input[] = {
//			0.269,
//			0.56,
//			0.941,
//			3.11,
//
//			0.0816649, // xp
//			0.0532961, // x0
//			0.0829682, // xt
//			-0.0659988, // yp
//			-0.084614,  // y0
//			+0.0373802, // yt
//
//			-0.102141,  // xpb
//			-0.0845715, // x0b
//			-0.0587413, // xtb
//			+0.0243354, // ypb
//			+0.0533635, // y0b
//			-0.0695015  // ytb
//	};

	RooRealVar ap ("ap","ap", par_input[0],0,0.5);
	RooRealVar apa("apa","apa", par_input[1],0,1);
	RooRealVar a0 ("a0","a0", par_input[2],0.8,1);
	RooRealVar a0a("a0a","a0a",0);
	RooFormulaVar at("at","sqrt(1-ap*ap-a0*a0)",RooArgSet(ap,a0));
	RooRealVar ata("ata","ata", par_input[3],2,4);

	RooRealVar xp ("xp","xp",par_input[4],-0.2,0.2);
	RooRealVar x0 ("x0","x0",par_input[5],-0.2,0.2);
	RooRealVar xt ("xt","xt",par_input[6],-0.2,0.2);

	RooRealVar yp ("yp","yp",par_input[7],-0.2,0.2);
	RooRealVar y0 ("y0","y0",par_input[8],-0.2,0.2);
	RooRealVar yt ("yt","yt",par_input[9],-0.2,0.2);

	RooRealVar xpb("xpb","xpb",par_input[10],-0.2,0.2);
	RooRealVar x0b("x0b","x0b",par_input[11],-0.2,0.2);
	RooRealVar xtb("xtb","xtb",par_input[12],-0.2,0.2);

	RooRealVar ypb("ypb","ypb",par_input[13],-0.2,0.2);
	RooRealVar y0b("y0b","y0b",par_input[14],-0.2,0.2);
	RooRealVar ytb("ytb","ytb",par_input[15],-0.2,0.2);

	DtCPPDF mixing_pdf_a("mixing_pdf_a", "mixing_pdf_a", true, perfect_tagging_,
			*thetat_,
			*thetab_,
			*phit_,

			ap,
			apa,
			a0,
			ata,
			xp,
			x0,
			xt,
			yp,
			y0,
			yt,

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

	DtCPPDF mixing_pdf_ab("mixing_pdf_ab", "mixing_pdf_ab", true, perfect_tagging_,
			*thetat_,
			*thetab_,
			*phit_,

			ap,
			apa,
			a0,
			ata,
			xpb,
			x0b,
			xtb,
			ypb,
			y0b,
			ytb,

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

	DtCPPDF mixing_pdf_b("mixing_pdf_b", "mixing_pdf_b", false, perfect_tagging_,
			*thetat_,
			*thetab_,
			*phit_,

			ap,
			apa,
			a0,
			ata,
			xp,
			x0,
			xt,
			yp,
			y0,
			yt,

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

	DtCPPDF mixing_pdf_bb("mixing_pdf_bb", "mixing_pdf_bb", false, perfect_tagging_,
			*thetat_,
			*thetab_,
			*phit_,

			ap,
			apa,
			a0,
			ata,
			xpb,
			x0b,
			xtb,
			ypb,
			y0b,
			ytb,

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


	dt_->setRange(-15, 15);

	RooArgSet varsToGenerate (*dt_, *thetat_, *thetab_, *phit_);
	RooArgSet varsToAdd(
			*tagwtag_,
			*expmc_,
			*expno_,
			*shcosthb_,
			*benergy_,
			*mbc_,
			*vrntrk_,
			*vrzerr_,
			*vrchi2_);
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

    const int num_fav = num_events*frac_fav/2.0;
    const int num_sup = num_events*frac_sup/2.0;

    printf("INFO: Creating GenSpecs...\n");
	RooAbsPdf::GenSpec* genSpec_a = mixing_pdf_a.prepareMultiGen(varsToGenerate, RooFit::NumEvents(num_fav));
    printf("INFO: 1/4 GenSpecs ready\n");
	RooAbsPdf::GenSpec* genSpec_ab = mixing_pdf_ab.prepareMultiGen(varsToGenerate, RooFit::NumEvents(num_fav));
    printf("INFO: 2/4 GenSpecs ready\n");
	RooAbsPdf::GenSpec* genSpec_b = mixing_pdf_b.prepareMultiGen(varsToGenerate, RooFit::NumEvents(num_sup));
    printf("INFO: 3/4 GenSpecs ready\n");
	RooAbsPdf::GenSpec* genSpec_bb = mixing_pdf_bb.prepareMultiGen(varsToGenerate, RooFit::NumEvents(num_sup));
    printf("INFO: 4/4 GenSpecs ready\n");

	RooAbsData::setDefaultStorageType(RooAbsData::Tree);
	RooDataSet* dataset;
	TString filename;

	RooRealVar btagmcli("btagmcli", "btagmcli", -511, 511);
	RooRealVar brecflav("brecflav", "brecflav", -1, 1);

	for (int i = 1; i <= num_toys; i++) {
	//	RooDataSet* dataset_a = mixing_pdf_a.generate(varsToGenerate, RooFit::NumEvents(num_fav));
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

	//	RooDataSet* dataset_ab = mixing_pdf_ab.generate(varsToGenerate, RooFit::NumEvents(num_fav));
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

	//	RooDataSet* dataset_b = mixing_pdf_b.generate(varsToGenerate, RooFit::NumEvents(num_sup));
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

	//	RooDataSet* dataset_bb = mixing_pdf_bb.generate(varsToGenerate, RooFit::NumEvents(num_sup));
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

		dt_formula_ = new RooFormulaVar( "dt", "#Deltat [ps]", "(vrvtxz-vtvtxz)/(0.425*0.0299792458)", RooArgSet(*vrvtxz_, *vtvtxz_));
		vrzerr_formula_ = new RooFormulaVar( "vrzerr", "#sigma z_{rec} [cm]", "sqrt(vrerr6)", RooArgSet(*vrerr6_));
		vtzerr_formula_ = new RooFormulaVar( "vtzerr", "#sigma z_{tag} [cm]", "sqrt(vterr6)", RooArgSet(*vterr6_));

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
		TFile f(filename,"RECREATE");
		new_tree->Write();
		f.Close();

		delete dataset;
		//dataset_ = dataset;
	}



//	return dataset;

//	result_ = mixing_pdf_a.fitTo(*new_dataset, RooFit::ConditionalObservables(conditional_vars_argset_),
//			RooFit::Minimizer("Minuit2"), RooFit::Hesse(false), RooFit::Minos(false),
//			RooFit::Range("dtFitRange"), RooFit::Save(true), RooFit::NumCPU(num_CPUs_));
}

/**
 * Make a plot of a variable with both data points and pdf projection and display a pull plot
 * beneath it. Then save to a file.
 *
 * @param var Variable to be plotted
 * @param data Dataset against which to plot
 * @param pdf PDF to use for plotting and pull calculation
 * @param title [optional] Title for the y-axis
 */
void FitterCPV::PlotWithPull(const RooRealVar& var, const RooAbsData& data, const RooAbsPdf& pdf, const char* title) const {
	TString name = pdf.GetName();
	name += "_";
	name += var.GetName();
	TCanvas canvas(name, name, 500, 500);

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

	// Renormalization required for multi-CPU plots but not for single-CPU; possible RooFit bug
	if (num_CPUs_ > 1) {
		pdf.plotOn(plot,RooFit::ProjWData(conditional_vars_argset_, data, kFALSE), RooFit::NumCPU(num_CPUs_),
				RooFit::NormRange("dtFitRange"), RooFit::Normalization(1.0/data.numEntries()));
	} else {
		pdf.plotOn(plot,RooFit::ProjWData(conditional_vars_argset_, data, kFALSE), RooFit::NumCPU(num_CPUs_),
				RooFit::NormRange("dtFitRange"));
	}

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

/**
 * Create and return statbox with fit results to overlay on plots. Returns a NULL pointer
 * if no fit result exist.
 *
 * @param chi2 Reduced chi2 of the projection plot
 * @param position_top Should the box be displayed at the top or bottom of the plot
 * @param position_left Should the box be displayed at the left or right of the plot
 */
TPaveText* FitterCPV::CreateStatBox(const double chi2, const bool position_top, const bool position_left) const {
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

/**
 * Construct and return a cut string common for all for categories
 */
TString FitterCPV::GetCommonCutsString() const{
	TString common_cuts("evmcflag==1&&vrusable==1&&vtusable==1&&");
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
//	TTree* input_tree = dynamic_cast<TTree*>(input_file->Get("mixing_pdf_aData"));

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
	RooDataSet* temp_dataset = new RooDataSet("dataset", "dataset", input_tree, dataset_vars_argset_, common_cuts);

//	separate conditional vars from stuff like thetat

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
	vtzerr_ = static_cast<RooRealVar*>(dataset_->addColumn(*vtzerr_formula_));
	dt_ = static_cast<RooRealVar*>(dataset_->addColumn(*dt_formula_));

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
void FitterCPV::SetOutputDir(const char* output_dir) {
	gEnv->SetValue("Canvas.PrintDirectory", output_dir);
	printf("print dir: %s\n",gEnv->GetValue("Canvas.PrintDirectory","not found"));
	if (make_plots_) {
		output_file_ = new TFile(TString(output_dir) + "/plots.root", "RECREATE");
	}
}

