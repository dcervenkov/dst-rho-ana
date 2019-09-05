/**
 *  @file    angularpdf.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2018-06-19
 *
 *  @brief Time-independent angular PDF for B -> VV decays
 *
 */

#pragma once

// ROOT includes
#include "RooRealProxy.h"

// Local includes
#include "efficiency.h"

class AngularPDF : public RooAbsPdf {
   public:
    AngularPDF(){};
    AngularPDF(const char* name, const char* title, bool _B_bar, int _efficiency_model, std::vector<std::string> _efficiency_files,
               RooAbsReal& _tht,
               RooAbsReal& _thb,
               RooAbsReal& _phit,
               RooAbsReal& _ap,
               RooAbsReal& _apa,
               RooAbsReal& _a0,
               RooAbsReal& _ata);
    AngularPDF(const AngularPDF& other, const char* name = 0);
    virtual TObject* clone(const char* newname) const { return new AngularPDF(*this, newname); }
    // inline virtual ~AngularPDF() {}

    Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
                                const char* rangeName = 0) const;
    Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;

   protected:
    RooRealProxy tht;
    RooRealProxy thb;
    RooRealProxy phit;
    RooRealProxy ap;
    RooRealProxy apa;
    RooRealProxy a0;
    RooRealProxy ata;

    bool B_bar;
    int efficiency_model;
    Efficiency eff;

    // Angular terms from PDF x Efficiency
    Double_t f1(const double* vars);
    Double_t f2(const double* vars);
    Double_t f3(const double* vars);
    Double_t f4(const double* vars);
    Double_t f5(const double* vars);
    Double_t f6(const double* vars);

    // Numerical integrals of the above
    Double_t int_tht_thb_phit[6];

    Double_t evaluate() const;

   private:
    // ClassDef(AngularPDF,1) // Your description goes here...
};
