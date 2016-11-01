#ifndef FITTERTRANS_H
#define FITTERTRANS_H

#include "RooSimultaneous.h"
#include "RooChi2Var.h"

#include "Fitter.h"
#include "DSRhoPDFTDep.h"

class FitterTDep: public Fitter {
public:
	explicit FitterTDep(Double_t* outer_par_input);
	~FitterTDep();
	void ComputeChi2(const char* type);
	Double_t GetChi2(const char* type);
	void SaveParameters(char* file);
	RooRealVar* GetDt() { return dt; };
	DSRhoPDFTDep* GetPdf() { return pdf_a; };
	RooDataHist* GetBinnedDataSet();
	RooDataSet* GetReducedDataSet();
	void SaveResiduals();
	void GenerateDataSet(Int_t numEvents);

protected:
	void CreateReducedDataSet(const char* type);
	void CreateBinnedDataSet(const char* type);
	Double_t GetVPrecise(DSRhoPDFTDep* pdf);
	Double_t GetVPrecise1D(const int i, DSRhoPDFTDep* pdf, RooDataSet* loc_dataset);
	Double_t GetVPrecise1D(const int i, RooSimultaneous* spdf, RooDataSet* loc_dataset);

	Double_t par_input[16];
	RooSimultaneous* pdf;

	Int_t dt_bins;
	RooRealVar* dt;

	Int_t vars_bins[4];
	RooRealVar* vars[4];

	RooCategory* decType;

	/// Time-dep additional vars

	RooRealVar* dm;

	RooRealVar* xt;
	RooRealVar* xp;
	RooRealVar* x0;

	RooRealVar* yt;
	RooRealVar* yp;
	RooRealVar* y0;

	RooRealVar* xtb;
	RooRealVar* xpb;
	RooRealVar* x0b;

	RooRealVar* ytb;
	RooRealVar* ypb;
	RooRealVar* y0b;

	DSRhoPDFTDep* pdf_a;
	DSRhoPDFTDep* pdf_b;
	DSRhoPDFTDep* pdf_ab;
	DSRhoPDFTDep* pdf_bb;

};

#endif // FITTERTRANS_H
