#include "Fitter.h"

// ROOT includes
#include "TMath.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TPaveText.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooRandom.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooHist.h"

// Local includes
#include "Constants.h"
#include "DSRhoPDF.h"
#include "tools.h"

Fitter::Fitter(Double_t* outer_par_input) {
	for (int i = 0; i < 11; i++)
		par_input[i] = outer_par_input[i];

	chi2Var = 0;
	result = 0;

	tht_bins = 60;
	thb_bins = 60;
	phit_bins = 60;

	vars_bins[0] = tht_bins;
	vars_bins[1] = thb_bins;
	vars_bins[2] = phit_bins;

	evmcflag = new RooRealVar("evmcflag", "evmcflag", -100, 100);

	tht = new RooRealVar("thetat", "#theta_{t}", 0, PI);
	thb = new RooRealVar("thetab", "#theta_{b}", kThetaBMin, kThetaBMax);
	phit = new RooRealVar("phit", "#phi_{t}", -PI, PI);

	vars[0] = tht;
	vars[1] = thb;
	vars[2] = phit;

	for (int i = 0; i < 3; i++) {
		vars[i]->setBins(vars_bins[i]);
	}

	ap = new RooRealVar("ap", "ap", par_input[0], 0, 0.5);
	apa = new RooRealVar("apa", "apa", par_input[1], 0, 1);
	apr = new RooFormulaVar("apr", "ap*cos(apa)", RooArgSet(*ap, *apa));
	api = new RooFormulaVar("api", "ap*sin(apa)", RooArgSet(*ap, *apa));
	a0 = new RooRealVar("a0", "a0", par_input[2], 0.8, 1);
	a0a = new RooRealVar("a0a", "a0a", 0);
	a0r = new RooFormulaVar("a0r", "a0*cos(a0a)", RooArgSet(*a0, *a0a));
	a0i = new RooFormulaVar("a0i", "a0*sin(a0a)", RooArgSet(*a0, *a0a));
	at = new RooFormulaVar("at", "sqrt(1-ap*ap-a0*a0)", RooArgSet(*ap, *a0));
	ata = new RooRealVar("ata", "ata", par_input[3], 2, 4);
	atr = new RooFormulaVar("atr", "at*cos(ata)", RooArgSet(*at, *ata));
	ati = new RooFormulaVar("ati", "at*sin(ata)", RooArgSet(*at, *ata));

	ap0r = new RooFormulaVar("ap0r", "ap*a0*cos(-apa+a0a)", RooArgSet(*ap, *apa, *a0, *a0a));
	a0ti = new RooFormulaVar("a0ti", "a0*at*sin(-a0a+ata)", RooArgSet(*a0, *a0a, *at, *ata));
	apti = new RooFormulaVar("apti", "ap*at*sin(-apa+ata)", RooArgSet(*ap, *apa, *at, *ata));

	varSet = new RooArgSet(*ap, *a0, *at, *ap0r, *a0ti, *apti);

	pdf = new DSRhoPDF("pdf", "pdf", *tht, *thb, *phit, *ap, *apa, *a0, *ata);

	parameters = new RooArgSet(*ap, *apa, *a0, *ata);

	variables = new RooArgList(*tht, *thb, *phit);

	/// numFitParameters holds # of NON-constant fit parameters
	numFitParameters = (parameters->selectByAttrib("Constant", kFALSE))->getSize();

//    delete dataSet;
//    dataSet = pdf->generate(RooArgSet(*vars[0],*vars[1],*vars[2]),600000);
//    dataSet->write("data/dataset_from_pdf");

}

Fitter::~Fitter() {
	delete thb;
	delete tht;
	delete phit;
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
	delete ap0r;
	delete a0ti;
	delete apti;
}

void Fitter::Fit() {
	numFitParameters = (parameters->selectByAttrib("Constant", kFALSE))->getSize();
	result = pdf->fitTo(*dataSet, RooFit::Save(), RooFit::Timer(true),
			RooFit::Minos(false), RooFit::Hesse(false), RooFit::Strategy(1), RooFit::NumCPU(4));

	const TMatrixDSym& cor = result->correlationMatrix();
	result->Print();
	cor.Print();
	//TCanvas c1;
	//result->correlationHist()->Draw("colz");
	//c1.SaveAs("corr.gif");
}

void Fitter::CreateBinnedDataSet() {
	RooRandom::randomGenerator()->SetSeed(0);
	binnedNumEntries = dataSet->numEntries();

	/// Create a binned dataSet which is needed for chi2 calculation
	dataSet_binned = new RooDataHist("dataSet_binned", "dataSet_binned", RooArgSet(*tht, *thb, *phit), *dataSet);

	//dataSet_binned = new RooDataHist("dataSet_binned","dataSet_binned",RooArgSet(*tht,*thb,*phit,*dt),*pdf_a);
	//RooDataHist* dataSet_binned = pdf->generateBinned(RooArgSet(var1,var2,var3),dataSet->numEntries(),kFALSE);
}

void Fitter::ComputeChi2() {
	//if(dataSet_binned == NULL)
	CreateBinnedDataSet();

	numFitParameters = (parameters->selectByAttrib("Constant", kFALSE))->getSize();

	chi2Var = new RooChi2Var("chi2Var", "chi2Var", *pdf, *dataSet_binned);

	RooRealVar* ndof = new RooRealVar("ndof", "number of degrees of freedom", 0);
	RooRealVar* chi2red = new RooRealVar("chi2red", "reduced chi^2", 0);
	RooRealVar* prob = new RooRealVar("prob", "prob(chi2,ndof)", 0);

	ndof->setVal(dataSet_binned->numEntries() - numFitParameters);
	chi2red->setVal(chi2Var->getVal() / ndof->getVal());
	prob->setVal(TMath::Prob(chi2Var->getVal(), static_cast<int>(ndof->getVal())));

	delete dataSet_binned;

	printf("chi2 = %f\nndof = %f\nchi2red = %f\nprob = %f\n\n", chi2Var->getVal(), ndof->getVal(), chi2red->getVal(),
			prob->getVal());

}

Double_t Fitter::GetVPrecise(DSRhoPDF* pdf) {
	/// The pdf seems to be varying quite rapidly at some places, so the approximation of constant pdf in a voxel is sometimes bad.
	/// This function gives a better approximation of a pdf value in a voxel by averaging through multiple points inside it.

	Double_t v = 0;

	Int_t tht_subbins = 3;
	Int_t thb_subbins = 3;
	Int_t phit_subbins = 3;

	Double_t binVolume = tht->getBinWidth(0) * thb->getBinWidth(0) * phit->getBinWidth(0);

	Double_t tht_binmin = tht->getVal() - tht->getBinWidth(0) / 2;
	Double_t tht_binmax = tht->getVal() + tht->getBinWidth(0) / 2;
	Double_t thb_binmin = thb->getVal() - thb->getBinWidth(0) / 2;
	Double_t thb_binmax = thb->getVal() + thb->getBinWidth(0) / 2;
	Double_t phit_binmin = phit->getVal() - phit->getBinWidth(0) / 2;
	Double_t phit_binmax = phit->getVal() + phit->getBinWidth(0) / 2;

//    Int_t numPasses = 0;
//    Int_t pass = 0;

	/// Using *_binmax - 0.001 because when one is at a boundary of e.g. thb, thb->getVal() < thb_binmax is never violated, even though it should be equal.
	for (*tht = tht_binmin + tht->getBinWidth(0) / (2 * tht_subbins); tht->getVal() < tht_binmax - 0.001;
			tht->setVal(tht->getVal() + tht->getBinWidth(0) / tht_subbins)) {
		for (*thb = thb_binmin + thb->getBinWidth(0) / (2 * thb_subbins); thb->getVal() < thb_binmax - 0.001;
				thb->setVal(thb->getVal() + thb->getBinWidth(0) / thb_subbins)) {
			for (*phit = phit_binmin + phit->getBinWidth(0) / (2 * phit_subbins); phit->getVal() < phit_binmax - 0.001;
					phit->setVal(phit->getVal() + phit->getBinWidth(0) / phit_subbins)) {
				v += pdf->getVal(RooArgSet(*tht, *thb, *phit)) * binVolume * binnedNumEntries;
			}
		}
	}

	/// These lines return the variables to their initial states, before calling GetVPrecise()
	tht->setVal(tht_binmin + tht->getBinWidth(0) / 2);
	thb->setVal(thb_binmin + thb->getBinWidth(0) / 2);
	phit->setVal(phit_binmin + phit->getBinWidth(0) / 2);

	v = v / (tht_subbins * thb_subbins * phit_subbins);

	return v;
}

Double_t Fitter::GetVPrecise1D(const int i, DSRhoPDF* pdf, RooDataSet* loc_dataset) {
	const Double_t num_subbins = 5;

	RooArgSet intSet;
	for (int j = 0; j < 3; j++)
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

void Fitter::SaveResiduals() {

	TFile* file = new TFile("plots/residuals.root", "RECREATE");
	TCanvas* c_residuals = new TCanvas("c2", "c2", 500, 500);

	Double_t n = 0;
	Double_t v = 0;

	/// Variables have to be binned to be able to call ->weight to get bin content
	for (int i = 0; i < 3; i++)
		vars[i]->setBins(vars_bins[i]);

	TString name;
	TString path;
	TH2F* h2_residual[3];
	TH1F* h1_residual_bar[3];
	TH1F* h1_pull[3];
	Double_t chi2[3] = { 0, 0, 0 };
	Int_t ndof[3] = { 0, 0, 0 };

	/// Loop for tht, thb and phit
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
			if (n == 0)
				continue;
			v = GetVPrecise1D(i, pdf, dataSet);

			h2_residual[i]->Fill(vars[i]->getVal(), (n - v) / sqrt(n));
			h1_residual_bar[i]->Fill(vars[i]->getVal(), (n - v) / sqrt(n));
			h1_pull[i]->Fill((n - v) / sqrt(n));

			chi2[i] += ((n - v) * (n - v)) / v;
			ndof[i]++;
		}

		delete dataSet_binned_1D;

		h2_residual[i]->GetXaxis()->SetTitle(vars[i]->GetTitle());
		h2_residual[i]->SetMarkerStyle(7);
		h2_residual[i]->SetTitle("");

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
		path += format;
		c_residuals->SaveAs(path);
		c_residuals->SetGrid(0, 0);

		/// These residual_bar histograms are mainly for debugging purposes, may be disabled
		h1_residual_bar[i]->GetXaxis()->SetTitle(vars[i]->GetTitle());
		h1_residual_bar[i]->SetTitle("");
		c_residuals->SetGrid();
		h1_residual_bar[i]->Draw();
		h1_residual_bar[i]->Write();
		path = "plots/residual_bar_";
		path += vars[i]->GetName();
		path += format;
		c_residuals->SaveAs(path);
		c_residuals->SetGrid(0, 0);

		h1_pull[i]->Fit("gaus");
		h1_pull[i]->GetXaxis()->SetTitle(vars[i]->GetTitle());
		h1_pull[i]->SetTitle("");
		h1_pull[i]->Draw();
		h1_pull[i]->Write();
		path = "plots/pull_";
		path += vars[i]->GetName();
		path += format;
		c_residuals->SaveAs(path);
	}

	file->Close();
	delete c_residuals;

	///This is outside of the preceding loop because it would be intersparsed by different messages
	for (int i = 0; i < 3; i++) {
		printf("%s\tchi2: %.2f\tndof: %i\tchi2red: %.3f\tprob: %f\n", vars[i]->GetName(), chi2[i], ndof[i],
				chi2[i] / ndof[i], TMath::Prob(chi2[i], ndof[i]));
	}

}

void Fitter::SaveNllPlot(RooRealVar* var) {
	const int steps = 100;
	Double_t orig_val = var->getVal();
	TCanvas c_nll("c_nll", "c_nll", 800, 600);
	TString name;
	TString basename;
	TString path;

	basename = "nll_";
	basename += var->GetName();
	name = "h1_" + basename;

	TH1F h1_nll(name, name, steps, var->getMin(), var->getMax());
	RooAbsReal* nll;
	nll = pdf->createNLL(*dataSet, RooFit::NumCPU(2));
	for (Int_t i = 0; i < steps; i++) {
		var->setVal(i * var->getMax() / steps + (var->getMax() - var->getMin()) / (2 * steps));
		printf("Computing %i/%i likelihood function.\n", i + 1, steps);
		h1_nll.Fill(var->getVal(), 2 * nll->getVal());
	}
	delete nll;
	h1_nll.GetXaxis()->SetTitle(var->GetName());
	h1_nll.SetStats(kFALSE);
	h1_nll.Draw();
	h1_nll.Write();
	path = "plots/";
	path += basename;
	path += format;
	c_nll.SaveAs(path);
	var->setVal(orig_val);
}

void Fitter::SaveNllPlot(RooRealVar* var1, RooRealVar* var2) {
	const int steps1 = 30;
	const int steps2 = 30;
	Double_t stepsize1 = (var1->getMax() - var1->getMin()) / steps1;
	Double_t stepsize2 = (var2->getMax() - var2->getMin()) / steps2;
	Double_t orig_val1 = var1->getVal();
	Double_t orig_val2 = var2->getVal();
	TCanvas c_nll("c_nll", "c_nll", 600, 600);
	TString name;
	TString basename;
	TString path;
	basename = "nll2_";
	basename += var1->GetName();
	basename += "_";
	basename += var2->GetName();
	name = "h2_" + basename;

	TH2F h2_nll(name, name, steps1, var1->getMin(), var1->getMax(), steps2, var2->getMin(), var2->getMax());
	RooAbsReal* nll;
	nll = pdf->createNLL(*dataSet, RooFit::NumCPU(2));
	for (Int_t i = 0; i < steps1; i++) {
		for (Int_t j = 0; j < steps2; j++) {
			var1->setVal(i * stepsize1 + var1->getMin() + stepsize1 / 2);
			var2->setVal(j * stepsize2 + var2->getMin() + stepsize2 / 2);
			printf("Computing %i/%i likelihood function.\n", i * steps2 + j + 1, steps1 * steps2);
			h2_nll.Fill(var1->getVal(), var2->getVal(), 2 * nll->getVal());
		}
	}

	delete nll;

	h2_nll.GetXaxis()->SetTitle(var1->GetName());
	h2_nll.GetYaxis()->SetTitle(var2->GetName());
	h2_nll.SetStats(kFALSE);
	h2_nll.SetOption("colz");
	h2_nll.Draw();
	h2_nll.Write();
	path = "plots/";
	path += basename;
	path += format;
	c_nll.SaveAs(path);
	var1->setVal(orig_val1);
	var2->setVal(orig_val2);
}

void Fitter::GetRecoveredParameters(Int_t& numParameters, Double_t** recoveredParameters) {
	RooRealVar* chi2red = new RooRealVar("chi2red", "reduced chi^2", 0);
	chi2red->setVal(chi2Var->getVal() / (dataSet_binned->numEntries() - numFitParameters));

	numParameters = 17;
	Double_t* parameters = new Double_t[numParameters];

	parameters[0] = chi2red->getVal();
	parameters[1] = ap->getVal();
	parameters[2] = ap->getError();
	parameters[3] = apa->getVal();
	parameters[4] = apa->getError();
	parameters[5] = a0->getVal();
	parameters[6] = a0->getError();
	parameters[7] = a0a->getVal();
	parameters[8] = a0a->getError();
	parameters[9] = at->getVal();
	parameters[10] = at->getPropagatedError(*result);
	parameters[11] = ata->getVal();
	parameters[12] = ata->getError();
	parameters[13] = par_input[0];
	parameters[14] = par_input[1];
	parameters[15] = par_input[2];
	parameters[16] = par_input[3];

	*recoveredParameters = parameters;
}

RooDataHist* Fitter::GetBinnedDataSet() {
	if (dataSet_binned == NULL)
		CreateBinnedDataSet();

	return dataSet_binned;
}

void Fitter::FixAllParameters() {
	RooRealVar* rooPar = 0;
	TIterator* parIter = parameters->createIterator();
	while ((rooPar = (RooRealVar*) parIter->Next()))
		rooPar->setConstant();

}

void Fitter::FixParameter(const char* par) {
	RooRealVar* rooPar = 0;
	rooPar = (RooRealVar*) parameters->find(par);
	if (rooPar != 0)
		rooPar->setConstant();
}

void Fitter::FreeParameter(const char* par) {
	RooRealVar* rooPar = 0;
	rooPar = (RooRealVar*) parameters->find(par);
	if (rooPar != 0)
		rooPar->setConstant(kFALSE);
}

void Fitter::PrintParameter(const char* par) {
	RooRealVar* rooPar = 0;
	rooPar = (RooRealVar*) parameters->find(par);
	if (rooPar != 0)
		rooPar->Print();

	else if (strcmp(par, at->GetName()) == 0)
		printf("%s = %f +- %f\n", par, at->getVal(), at->getPropagatedError(*result));
	else if (strcmp(par, a0a->GetName()) == 0)
		a0a->Print();
	else if (strcmp(par, "hp") == 0) {
		RooFormulaVar hpr("hpr", "(apr + atr)/sqrt(2)", RooArgSet(*apr, *atr));
		RooFormulaVar hpi("hpi", "(api + ati)/sqrt(2)", RooArgSet(*api, *ati));
		RooFormulaVar hp("hp", "sqrt(hpr*hpr+hpi*hpi)", RooArgSet(hpr, hpi));
		RooFormulaVar hpa("hpa", "atan2(hpi,hpr)", RooArgSet(hpr, hpi));

		printf("%s = %f +- %f\n", "hp", hp.getVal(), hp.getPropagatedError(*result));
		printf("%s = %f +- %f\n", "hpa", hpa.getVal(), hpa.getPropagatedError(*result));
	} else if (strcmp(par, "hm") == 0) {
		RooFormulaVar hmr("hmr", "(apr - atr)/sqrt(2)", RooArgSet(*apr, *atr));
		RooFormulaVar hmi("hmi", "(api - ati)/sqrt(2)", RooArgSet(*api, *ati));
		RooFormulaVar hm("hm", "sqrt(hmr*hmr+hmi*hmi)", RooArgSet(hmr, hmi));
		RooFormulaVar hma("hma", "atan2(hmi,hmr)", RooArgSet(hmr, hmi));

		printf("%s = %f +- %f\n", "hm", hm.getVal(), hm.getPropagatedError(*result));
		printf("%s = %f +- %f\n", "hma", hma.getVal(), hma.getPropagatedError(*result));
	} else
		printf("ERROR: Parameter '%s' doesn't exist! Can't print its value.\n", par);
}

void Fitter::GetHelParameters(Double_t* params) {
	RooFormulaVar hpr("hpr", "(apr + atr)/sqrt(2)", RooArgSet(*apr, *atr));
	RooFormulaVar hpi("hpi", "(api + ati)/sqrt(2)", RooArgSet(*api, *ati));
	RooFormulaVar hp("hp", "sqrt(hpr*hpr+hpi*hpi)", RooArgSet(hpr, hpi));
	RooFormulaVar hpa("hpa", "atan2(hpi,hpr)", RooArgSet(hpr, hpi));

	RooFormulaVar hmr("hmr", "(apr - atr)/sqrt(2)", RooArgSet(*apr, *atr));
	RooFormulaVar hmi("hmi", "(api - ati)/sqrt(2)", RooArgSet(*api, *ati));
	RooFormulaVar hm("hm", "sqrt(hmr*hmr+hmi*hmi)", RooArgSet(hmr, hmi));
	RooFormulaVar hma("hma", "atan2(hmi,hmr)", RooArgSet(hmr, hmi));

	params[0] = hp.getVal();
	params[1] = hp.getPropagatedError(*result);
	params[2] = hpa.getVal();
	params[3] = hpa.getPropagatedError(*result);
	params[4] = a0->getVal();
	params[5] = a0->getError();
	params[6] = hm.getVal();
	params[7] = hm.getPropagatedError(*result);
	params[8] = hma.getVal();
	params[9] = hma.getPropagatedError(*result);
}

void Fitter::SaveNllPlot(const char* par) {
	RooRealVar* rooPar = 0;
	rooPar = (RooRealVar*) parameters->find(par);
	if (rooPar != 0)
		SaveNllPlot(rooPar);
	else
		printf("ERROR: Parameter '%s' doesn't exist! Can't save Nll plot.\n", par);
}

void Fitter::SaveNllPlot(const char* par1, const char* par2) {
	RooRealVar* rooPar1 = 0;
	RooRealVar* rooPar2 = 0;
	rooPar1 = (RooRealVar*) parameters->find(par1);
	rooPar2 = (RooRealVar*) parameters->find(par2);
	if (rooPar1 != 0 && rooPar2 != 0)
		SaveNllPlot(rooPar1, rooPar2);
	if (rooPar1 == 0)
		printf("ERROR: Parameter %s doesn't exist! Can't save Nll plot.\n", par1);
	if (rooPar2 == 0)
		printf("ERROR: Parameter %s doesn't exist! Can't save Nll plot.\n", par2);
}

void Fitter::ReadDataSet(const char* file) {
//	dataSet = RooDataSet::read(file, *variables);
	TChain* chain = new TChain("h2000");
	chain->Add(file);
	dataSet = new RooDataSet("dataSet", "dataSet", chain, RooArgSet(*tht, *thb, *phit, *evmcflag), "evmcflag==1");
	dataSet->Print();
}

void Fitter::SaveParameters(char* file) {
	const Int_t numParameters = 18;
	Double_t* parameters = new Double_t[numParameters];

	FILE* pFile;
	pFile = fopen(file, "w");
	if (pFile == NULL) {
		printf("ERROR: couldn't open file %s for writing!\n", file);
		delete[] parameters;
		return;
	}

	parameters[0] = par_input[0];
	parameters[1] = ap->getVal();
	parameters[2] = ap->getError();
	parameters[3] = par_input[1];
	parameters[4] = apa->getVal();
	parameters[5] = apa->getError();
	parameters[6] = par_input[2];
	parameters[7] = a0->getVal();
	parameters[8] = a0->getError();
	parameters[9] = 0;
	parameters[10] = a0a->getVal();
	parameters[11] = a0a->getError();
	parameters[12] = sqrt(1 - par_input[0] * par_input[0] - par_input[2] * par_input[2]);
	parameters[13] = at->getVal();
	if (result == 0) {
		parameters[14] = 0;
	} else {
		parameters[14] = at->getPropagatedError(*result);
	}
	parameters[15] = par_input[3];
	parameters[16] = ata->getVal();
	parameters[17] = ata->getError();

	Int_t separators[numParameters] = { 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0 };

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

void Fitter::SaveHelParameters(char* file) {
	const Int_t numParameters = 18;
	Double_t* parameters = new Double_t[numParameters];

	FILE* pFile;
	pFile = fopen(file, "w");
	if (pFile == NULL) {
		printf("ERROR: couldn't open file %s for writing!\n", file);
		delete[] parameters;
		return;
	}

	RooFormulaVar hp("hp", "sqrt((ap*ap + at*at + 2*ap*at*cos(apa-ata))/2)", RooArgSet(*ap, *apa, *at, *ata));
	RooFormulaVar hpa("hpa", "atan((ap*sin(apa) + at*sin(ata))/(ap*cos(apa) + at*cos(ata)))", RooArgSet(*ap, *apa, *at, *ata));

	RooFormulaVar hm("hm", "sqrt((ap*ap + at*at - 2*ap*at*cos(apa-ata))/2)", RooArgSet(*ap, *apa, *at, *ata));
	RooFormulaVar hma("hma", "atan((ap*sin(apa) - at*sin(ata))/(ap*cos(apa) - at*cos(ata)))", RooArgSet(*ap, *apa, *at, *ata));

	parameters[0] = par_input[0];
	parameters[1] = hp.getVal();
	parameters[2] = 0;
	parameters[3] = par_input[1];
	parameters[4] = hpa.getVal();
	parameters[5] = 0;
	parameters[6] = par_input[2];
	parameters[7] = a0->getVal();
	parameters[8] = a0->getError();
	parameters[9] = 0;
	parameters[10] = a0a->getVal();
	parameters[11] = a0a->getError();
	parameters[12] = sqrt(1 - par_input[0] * par_input[0] - par_input[2] * par_input[2]);
	parameters[13] = hm.getVal();
	parameters[14] = 0;
	parameters[15] = par_input[3];
	parameters[16] = hma.getVal();
	parameters[17] = 0;

	if (result) {
		parameters[2] = hp.getPropagatedError(*result);
		parameters[5] = hpa.getPropagatedError(*result);
		parameters[14] = hm.getPropagatedError(*result);
		parameters[17] = hma.getPropagatedError(*result);
	}

	Int_t separators[numParameters] = { 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0 };

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

void Fitter::SaveVarPlots() {
	const TString dir = "plots/";
	TString path;
	TString name;

	path = dir + "projections.root";
	TFile* file = new TFile(path, "RECREATE");

	TCanvas* c1 = new TCanvas("c1", "c1", 500, 500);
	TPad *pad_main = new TPad("pad_main", "pad_main", 0, 0.25, 1, 1);
	TPad *pad_pull = new TPad("pad_pull", "pad_pull", 0, 0, 1, 0.25);
	pad_main->Draw();
	pad_pull->Draw();

	pad_main->SetBottomMargin(0.0);
	pad_pull->SetTopMargin(0.0);
	pad_pull->SetBottomMargin(0.35);

	RooPlot* frame = 0;

	/// This is a quick and dirty solution to be able to loop through all variables (except dt, which is treated separetely)
	const int numVars = 3;

	/// Saving simple projections on each of the variables
	for (int i = 0; i < numVars; i++) {
		pad_main->cd();
		name = "proj_";
		name += vars[i]->GetName();
		frame = vars[i]->frame();
		dataSet->plotOn(frame, RooFit::Name("data"));
		//if(i == 1)
		pdf->plotOn(frame);	//, RooFit::Project(RooArgSet(*vars[0], *vars[1], *vars[2])));
		frame->SetName(name);
		frame->SetTitle("");
		// This line makes sure the 0 is not drawn as it would overlap with the lower pad
		frame->GetYaxis()->SetRangeUser(frame->GetMinimum() + 0.001, frame->GetMaximum());
		frame->GetYaxis()->SetTitle("");
		frame->Draw();
		TPaveText* stat_box = tools::CreateStatBox(frame->chiSquare(), NULL, 1, 1);
		stat_box->Draw();

		pad_pull->cd();
		// Create a new frame to draw the pull distribution and add the distribution to the frame
		RooPlot* plot_eff_pull_ = vars[i]->frame(RooFit::Title("Pull Distribution"));
		plot_eff_pull_->SetTitle("");
		RooHist* hpull = frame->pullHist();
		hpull->SetFillColor(kGray);
		// The only working way to get rid of error bars; HIST draw option doesn't work with RooPlot
		for (int i = 0; i < hpull->GetN(); i++) {
			hpull->SetPointError(i, 0, 0, 0, 0);
		}
		plot_eff_pull_->addPlotable(hpull, "B");
		// We plot again without bars, so the points are not half covered by bars as in case of "BP" draw option.
		// We need to create and plot a clone, because the ROOT object ownership is transfered to the RooPlot by addPlotable().
		// If we just added the hpull twice, we would get a segfault.
		RooHist* hpull_clone = static_cast<RooHist*>(hpull->Clone());
		plot_eff_pull_->addPlotable(hpull_clone, "P"); // We plot again without bars, so the points are not half covered by bars

		plot_eff_pull_->GetXaxis()->SetTickLength(0.03 * pad_main->GetAbsHNDC() / pad_pull->GetAbsHNDC());
		plot_eff_pull_->GetXaxis()->SetTitle(TString(vars[i]->GetTitle()) + " [rad]");
		plot_eff_pull_->GetXaxis()->SetTitleOffset(4);
		plot_eff_pull_->GetXaxis()->SetLabelOffset(0.01 * pad_main->GetAbsHNDC() / pad_pull->GetAbsHNDC());
		plot_eff_pull_->GetYaxis()->SetRangeUser(-5, 5);
		plot_eff_pull_->GetYaxis()->SetNdivisions(505);
		plot_eff_pull_->Draw();

		frame->Write();
		path = dir + name + format;
		c1->SaveAs(path);
		delete frame;
	}

	delete c1;

	/// Create a binned pdf with the same number of events as the data, so that the 2d plots of pdf
	/// and data are the same scale
	printf("Generation of binned dataset starting...\n");
	RooDataHist* pdf_binned = pdf->generateBinned(RooArgSet(*vars[0], *vars[1], *vars[2]), dataSet->numEntries(),
			kTRUE);
	printf("Generation of binned dataset finished\n");

	TH2* h2_pdf = 0;
	TH2* h2_data = 0;
	TH2* h2_pull = 0;
	TH2* h2_errors = 0;

	TCanvas* c2 = new TCanvas("c2", "c2", 500, 500);
	c2->SetRightMargin(0.11);

	// Saving projections of both data and pdf on 2 dimensions
	// as well as pullmaps
	for (int i = 0; i < numVars; i++) {
		for (int j = i + 1; j < numVars; j++) {
			name = "proj_";
			name += vars[i]->GetName();
			name += "_";
			name += vars[j]->GetName();
			h2_pdf = static_cast<TH2*>(pdf_binned->createHistogram(name + "_pdf", *vars[i], RooFit::YVar(*vars[j])));
			h2_pdf->SetOption("colz");
			h2_pdf->SetStats(kFALSE);
			h2_pdf->SetTitle("");
			h2_pdf->GetZaxis()->SetTitle("");
			h2_pdf->SetMinimum(0);
			h2_data = static_cast<TH2*>(dataSet->createHistogram(name + "_data", *vars[i], RooFit::YVar(*vars[j])));
			h2_data->SetOption("colz");
			h2_data->SetStats(kFALSE);
			h2_data->SetTitle("");
			h2_data->GetZaxis()->SetTitle("");
			h2_data->SetMinimum(0);
			h2_data->SetMaximum(h2_pdf->GetMaximum());
			h2_pdf->Draw();
			h2_pdf->Write();
			path = dir + name + "_pdf" + format;
			c2->SaveAs(path);
			h2_data->Draw();
			h2_data->Write();
			path = dir + name + "_data" + format;
			c2->SaveAs(path);

			h2_pull = static_cast<TH2*>(h2_data->Clone("h2_pull"));
			h2_pull->Add(h2_pdf, -1);

			h2_errors = static_cast<TH2*>(h2_data->Clone("h2_pull"));
			h2_errors->Reset();
			for (int i = 1; i <= h2_pdf->GetNbinsX(); i++) {
				for (int j = 1; j <= h2_pdf->GetNbinsY(); j++) {
					h2_errors->SetBinContent(i, j, h2_data->GetBinError(i,j));
				}
			}

			h2_pull->Divide(h2_errors);

			// Whiten bins with too few data entries
			for (int i = 1; i <= h2_pdf->GetNbinsX(); i++) {
				for (int j = 1; j <= h2_pdf->GetNbinsY(); j++) {
					if (h2_data->GetBinContent(i,j) == 0) {
						h2_pull->SetBinContent(i, j, -1000);
					}
				}
			}

			h2_pull->SetMinimum(-5);
			h2_pull->SetMaximum(5);
			h2_pull->SetOption("colz");
			h2_pull->Draw();
			h2_pull->Write();
			path = dir + name + "_pullmap" + format;
			c2->SaveAs(path);

			delete h2_pdf;
			delete h2_data;
			delete h2_pull;
		}
	}

	file->Close();
	delete pdf_binned;
	delete c2;
}

