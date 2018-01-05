/**
 *  @file    dtscfpdf.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2017-12-25
 *
 *  @brief Class that holds the full angular time-dependent self-cross-feed PDF
 *
 */

#ifndef DTSCFPDF_H_
#define DTSCFPDF_H_

// ROOT includes
#include "TF1.h"
#include "RooRealProxy.h"
#include "RooAbsPdf.h"

// BASF includes
#include "tatami/tatami.h"

// Local includes
#include "constants.h"
#include "efficiency.h"

class DtSCFPDF : public RooAbsPdf {
public:
    // Constructor without mixing and CPV parameters for lifetime fits
    DtSCFPDF(const char *name, const char *title,
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

    DtSCFPDF(const char *name, const char *title, bool B_bar, bool CKM_favored, bool perfect_tagging, const int efficiency_model,
            RooAbsReal& _tht,
            RooAbsReal& _thb,
            RooAbsReal& _phit,
            RooAbsReal& _ap,
            RooAbsReal& _apa,
            RooAbsReal& _a0,
            RooAbsReal& _ata,
            RooAbsReal& _xp,
            RooAbsReal& _x0,
            RooAbsReal& _xt,
            RooAbsReal& _yp,
            RooAbsReal& _y0,
            RooAbsReal& _yt,

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

    DtSCFPDF(const DtSCFPDF& other, const char* name=0) ;
    virtual TObject* clone(const char* newname) const {return new DtSCFPDF(*this,newname);}
    inline virtual ~DtSCFPDF() { }

    Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
    Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:

    RooRealProxy tht ;
    RooRealProxy thb ;
    RooRealProxy phit ;
    RooRealProxy ap ;
    RooRealProxy apa ;
    RooRealProxy a0 ;
    RooRealProxy ata ;
    RooRealProxy xp ;
    RooRealProxy x0 ;
    RooRealProxy xt ;
    RooRealProxy yp ;
    RooRealProxy y0 ;
    RooRealProxy yt ;

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
    bool IsTimeIntegrated(int code) const;
    int GetRBin(double r) const;
    double GetWTag(int expno, int rbin, bool mc) const;
    double GetDeltaWTag(int expno, int rbin, bool mc) const;
    void CalculateAmplitudeTerms(double& Ap2, double& A02, double& At2,
                                 double& Ap0r, double& A0ti, double& Apti,
                                 const double& constant, const double& cos, 
                                 const double& sin) const;

    int efficiency_model;
    bool mixing;
    bool B_bar;
    bool CKM_favored;
    bool perfect_tagging;

    Efficiency eff;

    // Angular terms from PDF x Efficiency
    Double_t f1 (const double * vars);
    Double_t f2 (const double * vars);
    Double_t f3 (const double * vars);
    Double_t f4 (const double * vars);
    Double_t f5 (const double * vars);
    Double_t f6 (const double * vars);

    // Numerical integrals of the above
    Double_t int_tht_thb_phit[6] ;

private:
    //ClassDef(DSRhoPDF,1) // Your description goes here...
};

#endif // DTSCFPDF_H_
