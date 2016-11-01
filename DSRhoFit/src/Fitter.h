#ifndef FITTERTRANSTINDEP_H
#define FITTERTRANSTINDEP_H

// ROOT includes
#include "RooChi2Var.h"

// Local includes
#include "DSRhoPDF.h"

class Fitter {
public:
	explicit Fitter(Double_t* outer_par_input);
	~Fitter();
	void Fit();
	void ComputeChi2();
	Double_t GetChi2() { return chi2Var->getVal(); };
	void SaveNllPlot(const char* par);
	void SaveNllPlot(const char* par1, const char* par2);
	void GetRecoveredParameters(Int_t& numParameters, Double_t** recoveredParameters);
	RooRealVar* GetTht() { return tht; };
	RooRealVar* GetThb() { return thb; };
	RooRealVar* GetPhit() {	return phit; };
	DSRhoPDF* GetPdf() { return pdf; };
	RooDataHist* GetBinnedDataSet();
	RooDataSet* GetReducedDataSet();
	RooDataSet* GetDataSet() { return dataSet; };
	void FixAllParameters();
	void FixParameter(const char* par);
	void FreeParameter(const char* par);
	void PrintParameter(const char* par);
	void GetHelParameters(Double_t* params);
	void SaveResiduals();
	void ReadDataSet(const char* file);
	void SaveParameters(char* file);
	void SaveHelParameters(char* file);
	void SaveVarPlots();

protected:
	void SaveNllPlot(RooRealVar* var);
	void SaveNllPlot(RooRealVar* var1, RooRealVar* var2);

	void CreateBinnedDataSet();
	Double_t GetVPrecise(DSRhoPDF* pdf);
	Double_t GetVPrecise1D(const int i, DSRhoPDF* pdf, RooDataSet* loc_dataset);

	RooDataSet* dataSet { NULL };
	RooDataSet* dataSet_reduced { NULL };
	RooDataHist* dataSet_binned { NULL };
	Int_t binnedNumEntries { 0 };
	Double_t par_input[11];
	Bool_t doFit { 0 };
	RooChi2Var* chi2Var;
	RooFitResult* result;
	RooArgSet* parameters;
	RooArgList* variables;

	Int_t tht_bins;
	Int_t thb_bins;
	Int_t phit_bins;
	Int_t vars_bins[3];
	RooRealVar* vars[3];

	RooRealVar* evmcflag;

	RooRealVar* tht;
	RooRealVar* thb;
	RooRealVar* phit;

	RooRealVar* ap;
	RooRealVar* apa;
	RooFormulaVar* apr;
	RooFormulaVar* api;
	RooRealVar* a0;
	RooRealVar* a0a;
	RooFormulaVar* a0r;
	RooFormulaVar* a0i;
	RooFormulaVar* at;
	RooRealVar* ata;
	RooFormulaVar* atr;
	RooFormulaVar* ati;

	RooFormulaVar* ap0r;
	RooFormulaVar* a0ti;
	RooFormulaVar* apti;

	RooArgSet* varSet;

	/// numFitParameters holds # of NON-constant fit parameters
	RooArgSet* fitParameters { NULL };
	Int_t numFitParameters;

	DSRhoPDF* pdf;

};

#endif // FITTERTRANSTINDEP_H
