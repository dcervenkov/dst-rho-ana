/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

#ifndef DTBKG_H_
#define DTBKG_H_

// ROOT includes
#include "TF1.h"
#include "RooRealProxy.h"
#include "RooAbsPdf.h"

// Local includes
#include "constants.h"

class DtBKG : public RooAbsPdf {
public:
    // Constructor without mixing parameters for lifetime fits
    DtBKG(const char *name, const char *title,
            RooAbsReal& _dt,
            RooAbsReal& _vrerr6,
            RooAbsReal& _vterr6,
            RooAbsReal& _tau,
            RooAbsReal& _f_delta,
            RooAbsReal& _mu_delta,
            RooAbsReal& _mu_lifetime,
            RooAbsReal& _f_tail,
            RooAbsReal& _S_main,
            RooAbsReal& _S_tail,
            RooAbsReal& _f_outlier,
            RooAbsReal& _S_outlier
         );

    DtBKG(const DtBKG& other, const char* name=0) ;
    virtual TObject* clone(const char* newname) const {return new DtBKG(*this,newname);}
    inline virtual ~DtBKG() { }

    Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
    Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:

    RooRealProxy dt ;
    RooRealProxy vrerr6 ;
    RooRealProxy vterr6 ;
    RooRealProxy tau ;
    RooRealProxy f_delta ;
    RooRealProxy mu_delta ;
    RooRealProxy mu_lifetime ;
    RooRealProxy f_tail ;
    RooRealProxy S_main ;
    RooRealProxy S_tail ;
    RooRealProxy f_outlier ;
    RooRealProxy S_outlier ;

    Double_t evaluate() const;

private:
    //ClassDef(DSRhoPDF,1) // Your description goes here...
};

#endif // DTBKG_H_
