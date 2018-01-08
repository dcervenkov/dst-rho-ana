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
#include "RooBifurGauss.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooBifurGauss.h"
#include "RooPolynomial.h"
#include "RooGenericPdf.h"
#include "RooExponential.h"
#include "RooFormulaVar.h"

// BASF includes
#include "tatami/tatami.h"

// Local includes
#include "constants.h"

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

    DtSCFPDF(const char *name, const char *title, bool B_bar, bool CKM_favored, bool perfect_tagging,
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

    // Self-cross-feed phit model
    RooRealVar* scf_phit_poly_p2_;
    RooRealVar* scf_phit_f_;
    RooPolynomial* scf_phit_poly_;
    RooRealVar* scf_phit_offset_;
    RooFormulaVar* scf_phit_phit_;
    RooGenericPdf* scf_phit_cos_;
    RooAddPdf* scf_phit_model_;

    // Self-cross-feed thetat model
    RooRealVar* scf_thetat_f_;
    RooFormulaVar* scf_thetat_thetat_;
    RooGenericPdf* scf_thetat_model_;

    // Self-cross-feed thetab model
    RooRealVar* scf_thetab_gaus_mu_;
    RooRealVar* scf_thetab_gaus_sigma_l_;
    RooRealVar* scf_thetab_gaus_sigma_r_;
    RooBifurGauss* scf_thetab_gaus_;
    RooRealVar* scf_thetab_exp_alpha_;
    RooExponential* scf_thetab_exp_;
    RooRealVar* scf_thetab_f_;
    RooAddPdf* scf_thetab_model_;

    Double_t evaluate() const;
    bool IsTimeIntegrated(int code) const;
    int GetRBin(double r) const;
    double GetWTag(int expno, int rbin, bool mc) const;
    double GetDeltaWTag(int expno, int rbin, bool mc) const;
    void CalculateAmplitudeTerms(double& Ap2, double& A02, double& At2,
                                 const double& constant, const double& cos, 
                                 const double& sin) const;

    bool mixing;
    bool B_bar;
    bool CKM_favored;
    bool perfect_tagging;


private:
    //ClassDef(DSRhoPDF,1) // Your description goes here...
};

#endif // DTSCFPDF_H_
