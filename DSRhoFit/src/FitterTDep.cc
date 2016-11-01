#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooCategory.h"
#include "RooGenericPdf.h"
#include "Constants.h"
#include "RooDataSet.h"
#include "RooRandom.h"
#include "RooChi2Var.h"
#include "RooSimultaneous.h"
#include "RooPlot.h"
#include "RooFitResult.h"

#include "Minuit2/Minuit2Minimizer.h"
#include "TPluginManager.h"
#include "TMath.h"
#include "TIterator.h"
#include "TLine.h"

#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"

#include "DSRhoPDFTDep.h"
#include "Fitter.h"
#include "FitterTDep.h"

FitterTDep::FitterTDep(Double_t* outer_par_input) :
		Fitter(outer_par_input) {
	for (Int_t i = 0; i < 16; i++)
		par_input[i] = outer_par_input[i];

	dt_bins = 80;
	dt = new RooRealVar("dt", "dt", -3, 3);

	/// Convenient shorthands that can be used in loops
	vars[0] = tht;
	vars[1] = thb;
	vars[2] = phit;
	vars[3] = dt;
	vars_bins[0] = tht_bins;
	vars_bins[1] = thb_bins;
	vars_bins[2] = phit_bins;
	vars_bins[3] = dt_bins;

	decType = new RooCategory("decType", "decType");
	decType->defineType("a", 1);
	decType->defineType("ab", 2);
	decType->defineType("b", 3);
	decType->defineType("bb", 4);

	/// Time-dep additional vars

	dm = new RooRealVar("dm", "dm", 0.507e12);

	xp = new RooRealVar("xp", "xp", par_input[4], -0.2, 0.2);
	x0 = new RooRealVar("x0", "x0", par_input[5], -0.2, 0.2);
	xt = new RooRealVar("xt", "xt", par_input[6], -0.2, 0.2);

	yp = new RooRealVar("yp", "yp", par_input[7], -0.2, 0.2);
	y0 = new RooRealVar("y0", "y0", par_input[8], -0.2, 0.2);
	yt = new RooRealVar("yt", "yt", par_input[9], -0.2, 0.2);

	xpb = new RooRealVar("xpb", "xpb", par_input[10], -0.2, 0.2);
	x0b = new RooRealVar("x0b", "x0b", par_input[11], -0.2, 0.2);
	xtb = new RooRealVar("xtb", "xtb", par_input[12], -0.2, 0.2);

	ypb = new RooRealVar("ypb", "ypb", par_input[13], -0.2, 0.2);
	y0b = new RooRealVar("y0b", "y0b", par_input[14], -0.2, 0.2);
	ytb = new RooRealVar("ytb", "ytb", par_input[15], -0.2, 0.2);

	pdf_a = new DSRhoPDFTDep("pdf_a", "pdf_a", "a", *tht, *thb, *phit, *dt, *ap, *apa, *a0, *ata, *xp, *x0, *xt, *yp, *y0, *yt);
	pdf_b = new DSRhoPDFTDep("pdf_b", "pdf_b", "b", *tht, *thb, *phit, *dt, *ap, *apa, *a0, *ata, *xpb, *x0b, *xtb, *ypb, *y0b, *ytb);
	pdf_ab = new DSRhoPDFTDep("pdf_ab", "pdf_ab", "ab", *tht, *thb, *phit, *dt, *ap, *apa, *a0, *ata, *xpb, *x0b, *xtb,	*ypb, *y0b, *ytb);
	pdf_bb = new DSRhoPDFTDep("pdf_bb", "pdf_bb", "bb", *tht, *thb, *phit, *dt, *ap, *apa, *a0, *ata, *xp, *x0, *xt, *yp, *y0, *yt);

	pdf = new RooSimultaneous("simPdf", "simPdf", *decType);
	pdf->addPdf(*pdf_a, "a");
	pdf->addPdf(*pdf_ab, "ab");
	pdf->addPdf(*pdf_b, "b");
	pdf->addPdf(*pdf_bb, "bb");

	parameters->add(*xp);
	parameters->add(*x0);
	parameters->add(*xt);
	parameters->add(*yp);
	parameters->add(*y0);
	parameters->add(*yt);
	parameters->add(*xpb);
	parameters->add(*x0b);
	parameters->add(*xtb);
	parameters->add(*ypb);
	parameters->add(*y0b);
	parameters->add(*ytb);

	variables->add(*dt);
	variables->add(*decType);

	/// numFitParameters holds # of NON-constant fit parameters
	numFitParameters = (parameters->selectByAttrib("Constant", kFALSE))->getSize();

}

FitterTDep::~FitterTDep() {
	delete gPluginMgr;
	delete thb;
	delete tht;
	delete phit;
	delete dt;
	delete decType;

	delete ap;
	delete apa;
	delete apr;
	delete api;
	delete a0;
	delete a0a;
	delete a0r;
	delete a0i;
	delete at;
	delete ata;
	delete atr;
	delete ati;

	delete dm;
	delete xp;
	delete x0;
	delete xt;
	delete yp;
	delete y0;
	delete yt;
	delete xpb;
	delete x0b;
	delete xtb;
	delete ypb;
	delete y0b;
	delete ytb;
	delete pdf_a;
	delete pdf_b;
	delete pdf_ab;
	delete pdf_bb;
	delete pdf;
	delete parameters;
	delete fitParameters;
}

void FitterTDep::GenerateDataSet(Int_t numEvents) {
	RooRandom::randomGenerator()->SetSeed(0);

	/// For some reason this is the generated fraction of favored decays, have to check why
	const Double_t fracFav = 0.812;
	const Double_t fracSup = 1 - fracFav;

	if (dataSet != 0)
		delete dataSet;

	RooDataSet* temp_dataSet;
	temp_dataSet = pdf_a->generate(RooArgSet(*tht, *thb, *phit, *dt), TMath::Nint(numEvents * fracFav / 2));
	decType->setLabel("a");
	temp_dataSet->addColumn(*decType);

	dataSet = (RooDataSet*) temp_dataSet->Clone();
	delete temp_dataSet;

	temp_dataSet = pdf_ab->generate(RooArgSet(*tht, *thb, *phit, *dt), TMath::Nint(numEvents * fracFav / 2));
	decType->setLabel("ab");
	temp_dataSet->addColumn(*decType);
	dataSet->append(*temp_dataSet);
	delete temp_dataSet;

	temp_dataSet = pdf_b->generate(RooArgSet(*tht, *thb, *phit, *dt), TMath::Nint(numEvents * fracSup / 2));
	decType->setLabel("b");
	temp_dataSet->addColumn(*decType);
	dataSet->append(*temp_dataSet);
	delete temp_dataSet;

	temp_dataSet = pdf_bb->generate(RooArgSet(*tht, *thb, *phit, *dt), TMath::Nint(numEvents * fracSup / 2));
	decType->setLabel("bb");
	temp_dataSet->addColumn(*decType);
	dataSet->append(*temp_dataSet);
	delete temp_dataSet;
}

void FitterTDep::CreateBinnedDataSet(const char* type) {
	tht->setBins(tht_bins);
	thb->setBins(thb_bins);
	phit->setBins(phit_bins);
	dt->setBins(dt_bins);

	TString cut = "decType==decType::";
	cut += type;

	RooRandom::randomGenerator()->SetSeed(0);
	dataSet_reduced = (RooDataSet*) dataSet->reduce(cut);
	binnedNumEntries = dataSet_reduced->numEntries();

	/// Create a binned dataSet which is needed for chi2 calculation
	dataSet_binned = new RooDataHist("dataSet_binned", "dataSet_binned", RooArgSet(*tht, *thb, *phit, *dt),
			*dataSet_reduced);

	//dataSet_binned = new RooDataHist("dataSet_binned","dataSet_binned",RooArgSet(*tht,*thb,*phit,*dt),*pdf_a);
	//RooDataHist* dataSet_binned = pdf->generateBinned(RooArgSet(var1,var2,var3),dataSet->numEntries(),kFALSE);
}

void FitterTDep::CreateReducedDataSet(const char* type) {
	TString cut = "decType==decType::";
	cut += type;
	RooRandom::randomGenerator()->SetSeed(0);
	dataSet_reduced = (RooDataSet*) dataSet->reduce(cut);
}

void FitterTDep::ComputeChi2(const char* type) {
	//if(dataSet_binned == NULL)
	CreateBinnedDataSet(type);

	numFitParameters = (parameters->selectByAttrib("Constant", kFALSE))->getSize();

	if (strcmp(type, "a") == 0)
		chi2Var = new RooChi2Var("chi2Var", "chi2Var", *pdf_a, *dataSet_binned);
	else if (strcmp(type, "b") == 0)
		chi2Var = new RooChi2Var("chi2Var", "chi2Var", *pdf_b, *dataSet_binned);
	else if (strcmp(type, "ab") == 0)
		chi2Var = new RooChi2Var("chi2Var", "chi2Var", *pdf_ab, *dataSet_binned);
	else if (strcmp(type, "bb") == 0)
		chi2Var = new RooChi2Var("chi2Var", "chi2Var", *pdf_bb, *dataSet_binned);

	//chi2Var = new RooChi2Var("chi2Var","chi2Var",*simPdf,*dataSet_binned);

	RooRealVar* ndof = new RooRealVar("ndof", "number of degrees of freedom", 0);
	RooRealVar* chi2red = new RooRealVar("chi2red", "reduced chi^2", 0);
	RooRealVar* prob = new RooRealVar("prob", "prob(chi2,ndof)", 0);

	ndof->setVal(dataSet_binned->numEntries() - numFitParameters);
	chi2red->setVal(chi2Var->getVal() / ndof->getVal());
	prob->setVal(TMath::Prob(chi2Var->getVal(), static_cast<int>(ndof->getVal())));

	delete dataSet_binned;

	printf("%s:\tchi2 = %f\n%s:\tndof = %f\n%s:\tchi2red = %f\n%s:\tprob = %f\n\n", type, chi2Var->getVal(), type,
			ndof->getVal(), type, chi2red->getVal(), type, prob->getVal());

}

Double_t FitterTDep::GetChi2(const char* type) {
	Double_t mychi2 = 0;
	Double_t n = 0;
	Double_t v = 0;

//    TH1D* h_dchi2 = new TH1D("dchi2","dchi2",100,0,100000);

	//if(dataSet_binned == NULL)
	CreateBinnedDataSet(type);

	DSRhoPDFTDep* pdf = 0;

	if (strcmp(type, "a") == 0)
		pdf = pdf_a;
	else if (strcmp(type, "b") == 0)
		pdf = pdf_b;
	else if (strcmp(type, "ab") == 0)
		pdf = pdf_ab;
	else if (strcmp(type, "bb") == 0)
		pdf = pdf_bb;

	Double_t binVolume = tht->getBinWidth(0) * thb->getBinWidth(0) * phit->getBinWidth(0) * dt->getBinWidth(0);
	Int_t numBins = dataSet_binned->numEntries();

	Int_t numVPrecise = 0;

	/// Cycle through the centers of all bins
	/// I'm getting width of the first bin, because all bins are of equal width
	for (*tht = tht->getMin() + tht->getBinWidth(0) / 2; tht->getVal() < tht->getMax();
			tht->setVal(tht->getVal() + tht->getBinWidth(0))) {
		for (*thb = thb->getMin() + thb->getBinWidth(0) / 2; thb->getVal() < thb->getMax();
				thb->setVal(thb->getVal() + thb->getBinWidth(0))) {
			for (*phit = phit->getMin() + phit->getBinWidth(0) / 2; phit->getVal() < phit->getMax();
					phit->setVal(phit->getVal() + phit->getBinWidth(0))) {
				for (*dt = dt->getMin() + dt->getBinWidth(0) / 2; dt->getVal() < dt->getMax();
						dt->setVal(dt->getVal() + dt->getBinWidth(0))) {
					/// Weight is actually the bin content
					n = dataSet_binned->weight(RooArgSet(*tht, *thb, *phit, *dt), 0);
					if (n < 1) {
						numBins--;
						continue;
					}

					v = pdf->getVal(RooArgSet(*tht, *thb, *phit, *dt)) * binVolume * binnedNumEntries;

					if (((n - v) * (n - v) / v) > 1) {
						v = GetVPrecise(pdf);
						numVPrecise++;
					}

//                    printf("%.10f\t%.10f\n",v,GetVPrecise(pdf));

					mychi2 += (n - v) * (n - v) / v;
//                    h_dchi2->Fill((n-v)*(n-v)/v);
				}
			}
		}
	}

	printf("# VPrecise called: %i\n", numVPrecise);

//    TCanvas* c1 = new TCanvas("c1","c1",800,600);
//    c1->SetLogy();
//    h_dchi2->Draw();

	delete dataSet_binned;

	//printf("binVolume = %f\n",binVolume);
	printf("%s: numEntries = %i\n", type, binnedNumEntries);
	printf("%s: numBins = %i\n", type, numBins);
	printf("%s: mychi2 = %f\n", type, mychi2);
	printf("%s: mychi2red = %f\n", type, mychi2 / numBins);
	printf("%s: prob = %.10f\n\n", type, TMath::Prob(mychi2, numBins));
	return mychi2;
}

Double_t FitterTDep::GetVPrecise(DSRhoPDFTDep* pdf) {
	/// The pdf seems to be varying quite rapidly at some places, so the approximation of constant pdf in a voxel is sometimes bad.
	/// This function gives a better approximation of a pdf value in a voxel by averaging through multiple points inside it.

	Double_t v = 0;

	Int_t tht_subbins = 3;
	Int_t thb_subbins = 3;
	Int_t phit_subbins = 3;
	Int_t dt_subbins = 3;

	Double_t binVolume = tht->getBinWidth(0) * thb->getBinWidth(0) * phit->getBinWidth(0) * dt->getBinWidth(0);

	Double_t tht_binmin = tht->getVal() - tht->getBinWidth(0) / 2;
	Double_t tht_binmax = tht->getVal() + tht->getBinWidth(0) / 2;
	Double_t thb_binmin = thb->getVal() - thb->getBinWidth(0) / 2;
	Double_t thb_binmax = thb->getVal() + thb->getBinWidth(0) / 2;
	Double_t phit_binmin = phit->getVal() - phit->getBinWidth(0) / 2;
	Double_t phit_binmax = phit->getVal() + phit->getBinWidth(0) / 2;
	Double_t dt_binmin = dt->getVal() - dt->getBinWidth(0) / 2;
	Double_t dt_binmax = dt->getVal() + dt->getBinWidth(0) / 2;

//    Int_t numPasses = 0;
//    Int_t pass = 0;

	/// Using *_binmax - 0.001 because when one is at a boundary of e.g. thb, thb->getVal() < thb_binmax is never violated, even though it should be equal.
	for (*tht = tht_binmin + tht->getBinWidth(0) / (2 * tht_subbins); tht->getVal() < tht_binmax - 0.001;
			tht->setVal(tht->getVal() + tht->getBinWidth(0) / tht_subbins)) {
		for (*thb = thb_binmin + thb->getBinWidth(0) / (2 * thb_subbins); thb->getVal() < thb_binmax - 0.001;
				thb->setVal(thb->getVal() + thb->getBinWidth(0) / thb_subbins)) {
			for (*phit = phit_binmin + phit->getBinWidth(0) / (2 * phit_subbins); phit->getVal() < phit_binmax - 0.001;
					phit->setVal(phit->getVal() + phit->getBinWidth(0) / phit_subbins)) {
				for (*dt = dt_binmin + dt->getBinWidth(0) / (2 * dt_subbins); dt->getVal() < dt_binmax - 0.001;
						dt->setVal(dt->getVal() + dt->getBinWidth(0) / dt_subbins)) {
					v += pdf->getVal(RooArgSet(*tht, *thb, *phit, *dt)) * binVolume * binnedNumEntries;
				}
			}
		}
	}

	/// These lines return the variables to their initial states, before calling GetVPrecise()
	tht->setVal(tht_binmin + tht->getBinWidth(0) / 2);
	thb->setVal(thb_binmin + thb->getBinWidth(0) / 2);
	phit->setVal(phit_binmin + phit->getBinWidth(0) / 2);
	dt->setVal(dt_binmin + dt->getBinWidth(0) / 2);

	v = v / (tht_subbins * thb_subbins * phit_subbins * dt_subbins);

	return v;
}

Double_t FitterTDep::GetVPrecise1D(const int i, DSRhoPDFTDep* pdf, RooDataSet* loc_dataset) {
	const Double_t num_subbins = 5;

	RooArgSet intSet;
	for (int j = 0; j < 4; j++)
		if (j != i)
			intSet.add(*vars[j]);

	Double_t original_var = vars[i]->getVal();
	Double_t bin_min = original_var - vars[i]->getBinWidth(0) / 2;
	Double_t bin_max = original_var + vars[i]->getBinWidth(0) / 2;
	Double_t v = 0;
	RooAbsReal* vr;

	for (*vars[i] = bin_min + tht->getBinWidth(0) / (2 * num_subbins); vars[i]->getVal() < bin_max - 0.001;
			vars[i]->setVal(vars[i]->getVal() + vars[i]->getBinWidth(0) / num_subbins)) {
		vr = pdf->createIntegral(intSet, RooArgSet(*vars[0], *vars[1], *vars[2], *vars[3]));
		v += vr->getVal();
	}

	vars[i]->setVal(original_var);
	v *= vars[i]->getBinWidth(0) * loc_dataset->numEntries() / num_subbins;

	return v;
}

Double_t FitterTDep::GetVPrecise1D(const int i, RooSimultaneous* spdf, RooDataSet* loc_dataset) {
	const Double_t num_subbins = 5;

	RooArgSet intSet;
	for (int j = 0; j < 4; j++)
		if (j != i)
			intSet.add(*vars[j]);

	Double_t original_var = vars[i]->getVal();
	Double_t bin_min = original_var - vars[i]->getBinWidth(0) / 2;
	Double_t bin_max = original_var + vars[i]->getBinWidth(0) / 2;
	Double_t v = 0;
	RooAbsReal* vr;

	for (*vars[i] = bin_min + tht->getBinWidth(0) / (2 * num_subbins); vars[i]->getVal() < bin_max - 0.001;
			vars[i]->setVal(vars[i]->getVal() + vars[i]->getBinWidth(0) / num_subbins)) {
		vr = spdf->createIntegral(intSet, RooArgSet(*vars[0], *vars[1], *vars[2], *vars[3]));
		v += vr->getVal();
	}

	vars[i]->setVal(original_var);
	v *= vars[i]->getBinWidth(0) * loc_dataset->numEntries() / num_subbins;

	return v;
}

void FitterTDep::SaveResiduals() {

	TFile* file = new TFile("plots/residuals.root", "RECREATE");
	TCanvas* c_residuals = new TCanvas("c2", "c2", 800, 600);

	Double_t n = 0;
	Double_t v = 0;

	/// Variables have to be binned to be able to call ->weight to get bin content
	for (int i = 0; i < 4; i++)
		vars[i]->setBins(vars_bins[i]);

	TString name;
	TString path;
	TH2F* h2_residual[7];
	TH1F* h1_residual_bar[7];
	TH1F* h1_pull[7];
	Double_t chi2[7] = { 0, 0, 0, 0, 0, 0, 0 };
	Int_t ndof[7] = { 0, 0, 0, 0, 0, 0, 0 };

	/// Loop for tht, thb and phit. Loop for dt_{a,ab,b,bb} follows.
	for (int i = 0; i < 3; i++) {
		name = "h2_residual_";
		name += vars[i]->GetName();
		h2_residual[i] = new TH2F(name, name, vars_bins[i], vars[i]->getMin(), vars[i]->getMax(), 50, -5, 5);
		name = "h1_residual_bar_";
		name += vars[i]->GetName();
		h1_residual_bar[i] = new TH1F(name, name, vars_bins[i], vars[i]->getMin(), vars[i]->getMax());
		name = "h1_pull_";
		name += vars[i]->GetName();
		h1_pull[i] = new TH1F(name, name, 40, -7, 7);

		/// Binned dataset with *only one* dimension is created from the whole dataset, because tht,thb,phit distributions
		/// are the same for all four decay types.
		RooDataHist* dataSet_binned_1D = new RooDataHist("dataSet_binned_1D", "dataSet_binned_1D", RooArgSet(*vars[i]),
				*dataSet);
//        RooAbsReal* vr;

		for (*vars[i] = vars[i]->getMin() + vars[i]->getBinWidth(0) / 2; vars[i]->getVal() < vars[i]->getMax();
				vars[i]->setVal(vars[i]->getVal() + vars[i]->getBinWidth(0))) {
			n = dataSet_binned_1D->weight(RooArgSet(*vars[i]), 0);
			if (n <= 1)
				continue;
			v = GetVPrecise1D(i, pdf_a, dataSet);

			h2_residual[i]->Fill(vars[i]->getVal(), (n - v) / sqrt(n));
			h1_residual_bar[i]->Fill(vars[i]->getVal(), (n - v) / sqrt(n));
			h1_pull[i]->Fill((n - v) / sqrt(n));

			chi2[i] += ((n - v) * (n - v)) / v;
			ndof[i]++;
		}

		delete dataSet_binned_1D;

		h2_residual[i]->GetXaxis()->SetTitle(vars[i]->GetName());
		h2_residual[i]->SetMarkerStyle(7);

		/// Prepare 3-sigma lines for the residual histograms
		TLine baseline(vars[i]->getMin(), 0, vars[i]->getMax(), 0);
		TLine three_sigma_up(vars[i]->getMin(), 3, vars[i]->getMax(), 3);
		TLine three_sigma_down(vars[i]->getMin(), -3, vars[i]->getMax(), -3);
		three_sigma_up.SetLineColor(2);
		three_sigma_down.SetLineColor(2);

		c_residuals->SetGrid();
		h2_residual[i]->Draw();
		baseline.Draw();
		three_sigma_up.Draw();
		three_sigma_down.Draw();
		h2_residual[i]->Write();
		path = "plots/residual_";
		path += vars[i]->GetName();
		path += ".png";
		c_residuals->SaveAs(path);
		c_residuals->SetGrid(0, 0);

		/// These residual_bar histograms are mainly for debugging purposes, may be disabled
		h1_residual_bar[i]->GetXaxis()->SetTitle(vars[i]->GetName());
		c_residuals->SetGrid();
		h1_residual_bar[i]->Draw();
		h1_residual_bar[i]->Write();
		path = "plots/residual_bar_";
		path += vars[i]->GetName();
		path += ".png";
		c_residuals->SaveAs(path);
		c_residuals->SetGrid(0, 0);

		h1_pull[i]->Fit("gaus");
		h1_pull[i]->GetXaxis()->SetTitle(vars[i]->GetName());
		h1_pull[i]->Draw();
		h1_pull[i]->Write();
		path = "plots/pull_";
		path += vars[i]->GetName();
		path += ".png";
		c_residuals->SaveAs(path);
	}

	DSRhoPDFTDep* pdf = 0;
	/// Loop for dt_{a,ab,b,bb}
	for (int i = 3; i < 7; i++) {
		const char* type;

		switch (i) {
		case 3:
			type = "a";
			pdf = pdf_a;
			break;

		case 4:
			type = "ab";
			pdf = pdf_ab;
			break;

		case 5:
			type = "b";
			pdf = pdf_b;
			break;

		case 6:
			type = "bb";
			pdf = pdf_bb;
			break;

		default:
			break;
		}

		name = "h2_residual_";
		name += dt->GetName();
		name += "_";
		name += type;
		h2_residual[i] = new TH2F(name, name, dt_bins, dt->getMin(), dt->getMax(), 50, -5, 5);

		name = "h2_residual_bar_";
		name += dt->GetName();
		name += "_";
		name += type;
		h1_residual_bar[i] = new TH1F(name, name, dt_bins, dt->getMin(), dt->getMax());

		name = "h1_pull_";
		name += dt->GetName();
		name += "_";
		name += type;
		h1_pull[i] = new TH1F(name, name, 40, -7, 7);

		CreateReducedDataSet(type);
		RooDataHist* dataSet_binned_1D = new RooDataHist("dataSet_binned_1D", "dataSet_binned_1D", RooArgSet(*dt),
				*dataSet_reduced);

		for (*dt = dt->getMin() + dt->getBinWidth(0) / 2; dt->getVal() < dt->getMax();
				dt->setVal(dt->getVal() + dt->getBinWidth(0))) {
			/// It seems the first and last bin collect some under/overflow and are therefore skipped.
			/// TODO: Should be more thoroughly investigated.
			if (dt->getVal() < dt->getMin() + dt->getBinWidth(0))
				continue;
			else if (dt->getVal() > dt->getMax() - dt->getBinWidth(0))
				continue;

			n = dataSet_binned_1D->weight(RooArgSet(*dt), 0);
			if (n <= 10)
				continue;
			v = GetVPrecise1D(3, pdf, dataSet_reduced);

			h2_residual[i]->Fill(dt->getVal(), (n - v) / sqrt(n));
			h1_residual_bar[i]->Fill(dt->getVal(), (n - v) / sqrt(n));
			h1_pull[i]->Fill((n - v) / sqrt(n));

			chi2[i] += ((n - v) * (n - v)) / v;
			ndof[i]++;
		}

		delete dataSet_reduced;
		delete dataSet_binned_1D;

		h2_residual[i]->GetXaxis()->SetTitle(dt->GetName());
		h2_residual[i]->SetMarkerStyle(7);

		/// Prepare 3-sigma lines for the residual histograms
		TLine baseline(dt->getMin(), 0, dt->getMax(), 0);
		TLine three_sigma_up(dt->getMin(), 3, dt->getMax(), 3);
		TLine three_sigma_down(dt->getMin(), -3, dt->getMax(), -3);
		three_sigma_up.SetLineColor(2);
		three_sigma_down.SetLineColor(2);

		c_residuals->SetGrid();
		h2_residual[i]->Draw();
		baseline.Draw();
		three_sigma_up.Draw();
		three_sigma_down.Draw();
		h2_residual[i]->Write();
		path = "plots/residual_";
		path += dt->GetName();
		path += "_";
		path += type;
		path += ".png";
		c_residuals->SaveAs(path);
		c_residuals->SetGrid(0, 0);

		/// These residual_bar histograms are mainly for debugging purposes, may be disabled
		h1_residual_bar[i]->GetXaxis()->SetTitle(dt->GetName());
		c_residuals->SetGrid();
		h1_residual_bar[i]->Draw();
		h1_residual_bar[i]->Write();
		path = "plots/residual_bar_";
		path += dt->GetName();
		path += "_";
		path += type;
		path += ".png";
		c_residuals->SaveAs(path);
		c_residuals->SetGrid(0, 0);

		h1_pull[i]->Fit("gaus");
		h1_pull[i]->GetXaxis()->SetTitle(dt->GetName());
		h1_pull[i]->Draw();
		h1_pull[i]->Write();
		path = "plots/pull_";
		path += dt->GetName();
		path += "_";
		path += type;
		path += ".png";
		c_residuals->SaveAs(path);
		c_residuals->SetGrid(0, 0);
	}

	file->Close();
	delete c_residuals;

	///This is outside of the preceding loop because it would be intersparsed by different messages
	for (int i = 0; i < 3; i++)
		printf("%s\tchi2: %.2f\tndof: %i\tchi2red: %.3f\tprob: %f\n", vars[i]->GetName(), chi2[i], ndof[i],
				chi2[i] / ndof[i], TMath::Prob(chi2[i], ndof[i]));

	for (int i = 3; i < 7; i++)
		printf("%s_%i\tchi2: %.2f\tndof: %i\tchi2red: %.3f\tprob: %f\n", dt->GetName(), i, chi2[i], ndof[i],
				chi2[i] / ndof[i], TMath::Prob(chi2[i], ndof[i]));

}

void FitterTDep::SaveParameters(char* file) {
	DSRhoPDFTDep* pdf;
	DSRhoPDFTDep* pdfs[4] = { pdf_a, pdf_ab, pdf_b, pdf_bb };
	/// The next 2 lines enable getting category items' names and therefore reduced datasets in a loop
	const RooArgSet* args = dataSet->get();
	const RooCategory* cat = (RooCategory*) args->find("decType");
	RooDataSet* datacut;
	RooPlot* frame = 0;

	const Int_t numParameters = 58;
	Double_t* parameters = new Double_t[numParameters];

	const Int_t org_pdf_type = pdf_a->getType();

	/// Getting 1D chi^2 for all 4 decay types
	for (int i = 1; i <= 4; i++) {
		frame = dt->frame();
		TString type = (char*) cat->lookupType(i)->GetName();
		TString cut = "decType==decType::" + type;
		datacut = (RooDataSet*) dataSet->reduce(*dt, cut);
		datacut->plotOn(frame, RooFit::Name("data"));

		pdf = pdfs[i - 1];
		pdf->plotOn(frame, RooFit::Project(RooArgSet(*tht, *thb, *phit)));

		parameters[i - 1] = frame->chiSquare(11);

		delete frame;
	}

	pdf_a->setType(org_pdf_type);

	FILE* pFile;
	pFile = fopen(file, "w");
	if (pFile == NULL) {
		printf("ERROR: couldn't open file %s for writing!\n", file);
		delete[] parameters;
		return;
	}

	parameters[4] = par_input[0];
	parameters[5] = ap->getVal();
	parameters[6] = ap->getError();
	parameters[7] = par_input[1];
	parameters[8] = apa->getVal();
	parameters[9] = apa->getError();
	parameters[10] = par_input[2];
	parameters[11] = a0->getVal();
	parameters[12] = a0->getError();
	parameters[13] = 0;
	parameters[14] = a0a->getVal();
	parameters[15] = a0a->getError();
	parameters[16] = sqrt(1 - par_input[0] * par_input[0] - par_input[2] * par_input[2]);
	parameters[17] = at->getVal();
	if (result == 0) {
		parameters[18] = 0;
	} else {
		parameters[18] = at->getPropagatedError(*result);
	}
	parameters[19] = par_input[3];
	parameters[20] = ata->getVal();
	parameters[21] = ata->getError();
	parameters[22] = par_input[4];
	parameters[23] = xp->getVal();
	parameters[24] = xp->getError();
	parameters[25] = par_input[5];
	parameters[26] = x0->getVal();
	parameters[27] = x0->getError();
	parameters[28] = par_input[6];
	parameters[29] = xt->getVal();
	parameters[30] = xt->getError();
	parameters[31] = par_input[7];
	parameters[32] = yp->getVal();
	parameters[33] = yp->getError();
	parameters[34] = par_input[8];
	parameters[35] = y0->getVal();
	parameters[36] = y0->getError();
	parameters[37] = par_input[9];
	parameters[38] = yt->getVal();
	parameters[39] = yt->getError();
	parameters[40] = par_input[10];
	parameters[41] = xpb->getVal();
	parameters[42] = xpb->getError();
	parameters[43] = par_input[11];
	parameters[44] = x0b->getVal();
	parameters[45] = x0b->getError();
	parameters[46] = par_input[12];
	parameters[47] = xtb->getVal();
	parameters[48] = xtb->getError();
	parameters[49] = par_input[13];
	parameters[50] = ypb->getVal();
	parameters[51] = ypb->getError();
	parameters[52] = par_input[14];
	parameters[53] = y0b->getVal();
	parameters[54] = y0b->getError();
	parameters[55] = par_input[15];
	parameters[56] = ytb->getVal();
	parameters[57] = ytb->getError();
	Int_t separators[numParameters] = { 0, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0,
			1, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 0 };

	for (Int_t i = 0; i < numParameters; i++) {
		fprintf(pFile, "%+.5f ", parameters[i]);
		if (separators[i] == 1)
			fprintf(pFile, "| ");
		else if (separators[i] == 2)
			fprintf(pFile, "|| ");
	}
	fprintf(pFile, "\n");
	fclose(pFile);
	delete[] parameters;
}

RooDataHist* FitterTDep::GetBinnedDataSet() {
	if (dataSet_binned == NULL)
		CreateBinnedDataSet("a");

	return dataSet_binned;
}

RooDataSet* FitterTDep::GetReducedDataSet() {
	if (dataSet_reduced == NULL)
		CreateBinnedDataSet("a");

	return dataSet_reduced;
}

