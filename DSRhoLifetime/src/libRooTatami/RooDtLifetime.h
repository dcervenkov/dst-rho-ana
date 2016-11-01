/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

#ifndef ROODTLIFETIME
#define ROODTLIFETIME

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class RooDtLifetime : public RooAbsPdf {
public:

	RooDtLifetime() : is_mc(false), dt_ll(0), dt_ul(0), addoutlier(false), alpha(0) {};
	RooDtLifetime(const char *name, const char *title,
			RooAbsReal& _dt,
			RooAbsReal& _tau_b,
			RooAbsReal& _expno,
			RooAbsReal& _costh,
			RooAbsReal& _ecms,
			RooAbsReal& _rec_vtntrk,
			RooAbsReal& _rec_vterr,
			RooAbsReal& _rec_vtchi2,
			RooAbsReal& _rec_vtndf,
			RooAbsReal& _asc_vtntrk,
			RooAbsReal& _asc_vterr,
			RooAbsReal& _asc_vtchi2,
			RooAbsReal& _asc_vtndf,
			RooAbsReal& _keeptagl,
			bool _is_mc,
			bool _addoutlier = true,
			double _alpha = 1.0);
	RooDtLifetime(const char *name, const char *title,
			RooAbsReal& _dt,
			RooAbsReal& _tau_b,
			RooAbsReal& _expno,
			RooAbsReal& _costh,
			RooAbsReal& _ecms,
			RooAbsReal& _rec_vtntrk,
			RooAbsReal& _rec_vterr,
			RooAbsReal& _rec_vtchi2,
			RooAbsReal& _rec_vtndf,
			RooAbsReal& _asc_vtntrk,
			RooAbsReal& _asc_vterr,
			RooAbsReal& _asc_vtchi2,
			RooAbsReal& _asc_vtndf,
			RooAbsReal& _keeptagl,
			bool _is_mc,
			double _dt_ll,
			double _dt_ul,
			bool _addoutlier = true,
			double _alpha = 1.0);

	RooDtLifetime(const RooDtLifetime& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooDtLifetime(*this,newname); }
	inline virtual ~RooDtLifetime() { }

	Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
	Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:

	RooRealProxy dt;
	RooRealProxy tau_b;
	RooRealProxy expno;
	RooRealProxy costh;
	RooRealProxy ecms;
	RooRealProxy rec_vtntrk;
	RooRealProxy rec_vterr;
	RooRealProxy rec_vtchi2;
	RooRealProxy rec_vtndf;
	RooRealProxy asc_vtntrk;
	RooRealProxy asc_vterr;
	RooRealProxy asc_vtchi2;
	RooRealProxy asc_vtndf;
	RooRealProxy keeptagl;

	const bool is_mc;
	const double dt_ll;
	const double dt_ul;
	const bool addoutlier;
	const double alpha;

	Double_t evaluate() const;

private:

	ClassDef(RooDtLifetime,1) // Your description goes here...
};

#endif
