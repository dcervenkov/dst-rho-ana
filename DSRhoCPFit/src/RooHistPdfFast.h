/**
 *  @file    RooHistPdfFast.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2019-03-28
 *
 *  @brief Wrapper around RooHistPdf to make it *much* faster
 *
 *  When calling analyticalIntegral() of RooHistPdf, it calls sum() of the
 *  underlying RooDataHist *every* time, instead of having it cached. This
 *  makes it tremedously slow. This class is a workaround to make it as fast as
 *  it should be.
 *
 */

#pragma once

// ROOT includes
#include "RooHistPdf.h"
#include "RooRealProxy.h"

class RooHistPdfFast : public RooHistPdf {
    // using RooHistPdf::RooHistPdf;
   public:
    RooHistPdfFast(const char* name, const char* title, const RooArgSet& vars,
                   const RooDataHist& dhist, Int_t intOrder = 0);
    RooHistPdfFast(const RooHistPdfFast& other, const char* name);

    // Double_t evaluate() const;
    // Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,
    //                             const char* rangeName = 0) const;
    Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;
    // Double_t analyticalIntegralWN(Int_t code, const RooArgSet* normSet,
    //                               const char* rangeName = 0) const;
    // Double_t getValV(const RooArgSet* nset = 0) const;

    TObject* clone(const char* newname) const { return new RooHistPdfFast(*this, newname); }

   private:
    mutable Double_t cached_integral = 0;
};
