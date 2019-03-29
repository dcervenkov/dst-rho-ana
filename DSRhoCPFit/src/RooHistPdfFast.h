/**
 *  @file    RooHistPdfFast.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2019-03-28
 *
 *  @brief Wrapper around RooHistPdf to make it *much* faster
 *
 *  When calling analyticalIntegral() of RooHistPdf, it calls sum() of the
 *  underlying RooDataHist *every* time, instead of having it cached. This
 *  makes it tremendously slow when evaluating if the RooHistPdf can be cached (analyticalIntegral()
 *  is called for every event when doing that, even though it is always true). This class is a
 *  simple workaround to make it as fast as it should be.
 *
 */

#pragma once

// Standard includes
#include <string>

// ROOT includes
#include "RooHistPdf.h"

class RooHistPdfFast : public RooHistPdf {
   public:
    RooHistPdfFast(const char* name, const char* title, const RooArgSet& vars,
                   const RooDataHist& dhist, Int_t intOrder = 0)
        : RooHistPdf(name, title, vars, dhist, intOrder){};

    RooHistPdfFast(const RooHistPdfFast& other, const char* name)
        : RooHistPdf(other, name),
          cached_integral(other.cached_integral),
          cached_code(other.cached_code),
          cached_rangeName(other.cached_rangeName){};

    Double_t analyticalIntegral(Int_t code, const char* rangeName = 0) const;
    TObject* clone(const char* newname) const { return new RooHistPdfFast(*this, newname); }

   private:
    mutable Double_t cached_integral = 0;
    mutable Int_t cached_code = 0;
    mutable std::string cached_rangeName;
};
