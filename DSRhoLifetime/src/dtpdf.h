/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

#ifndef DTPDF_H_
#define DTPDF_H_

// ROOT includes
#include "TF1.h"
#include "RooRealProxy.h"
#include "RooAbsPdf.h"

// BASF includes
#include "tatami/tatami.h"

// Local includes
#include "constants.h"

class DtPDF : public RooAbsPdf {
public:
	// Constructor without mixing parameters for lifetime fits
	DtPDF(const char *name, const char *title,
			RooAbsReal& _dt,
			RooAbsReal& _tau,
			RooAbsReal& _expmc,
			RooAbsReal& _expno,
			RooAbsReal& _costhetabz,
			RooAbsReal& _benergy,
			RooAbsReal& _mbc,
			RooAbsReal& _vrntrk,
			RooAbsReal& _vrzerr,
			RooAbsReal& _vrchi2,
			RooAbsReal& _vrndf,
			RooAbsReal& _vtntrk,
			RooAbsReal& _vtzerr,
			RooAbsReal& _vtchi2,
			RooAbsReal& _vtndf,
			RooAbsReal& _vtistagl);

	DtPDF(const char *name, const char *title, bool CKM_favored, bool perfect_tagging,
			RooAbsReal& _S,
			RooAbsReal& _A,
			RooAbsReal& _wtag,
			RooAbsReal& _dt,
			RooAbsReal& _tau,
			RooAbsReal& _dm,
			RooAbsReal& _expmc,
			RooAbsReal& _expno,
			RooAbsReal& _costhetabz,
			RooAbsReal& _benergy,
			RooAbsReal& _mbc,
			RooAbsReal& _vrntrk,
			RooAbsReal& _vrzerr,
			RooAbsReal& _vrchi2,
			RooAbsReal& _vrndf,
			RooAbsReal& _vtntrk,
			RooAbsReal& _vtzerr,
			RooAbsReal& _vtchi2,
			RooAbsReal& _vtndf,
			RooAbsReal& _vtistagl);

	DtPDF(const DtPDF& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const {return new DtPDF(*this,newname);}
	inline virtual ~DtPDF() { }

	Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
	Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:

	RooRealProxy S ;
	RooRealProxy A ;
	RooRealProxy wtag ;
	RooRealProxy dt ;
	RooRealProxy tau ;
	RooRealProxy dm ;
	RooRealProxy expmc ;
	RooRealProxy expno ;
	RooRealProxy costhetabz ;
	RooRealProxy benergy ;
	RooRealProxy mbc ;
	RooRealProxy vrntrk ;
	RooRealProxy vrzerr ;
	RooRealProxy vrchi2 ;
	RooRealProxy vrndf ;
	RooRealProxy vtntrk ;
	RooRealProxy vtzerr ;
	RooRealProxy vtchi2 ;
	RooRealProxy vtndf ;
	RooRealProxy vtistagl ;

	Double_t evaluate() const;
	int GetRBin(double r) const;
	double GetWTag(int expno, int rbin, bool mc) const;
	double GetDeltaWTag(int expno, int rbin, bool mc) const;

	bool mixing;
	bool CKM_favored;
	bool perfect_tagging;

private:
	//ClassDef(DSRhoPDF,1) // Your description goes here...
};

#endif // DTPDF_H_
