/**
 *  @file    RooHistPdfFast.cc
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

#include "RooHistPdfFast.h"

// Standard includes
#include <string>

Double_t RooHistPdfFast::analyticalIntegral(Int_t code, const char* rangeName) const {
    // This is the first initialization of the cached integral; the integral
    // should be 0 only before the first initialization.
    if (cached_integral == 0) {
        cached_integral = RooHistPdf::analyticalIntegral(code, rangeName);
        cached_code = code;
        // We can't assign a NULL pointer to a std::string, therefore we assign
        // a null-terminated string instead
        cached_rangeName = rangeName ? rangeName : "";
    // If another integral/range than the currently cached one is requested,
    // recompute it.
    } else if (code != cached_code || cached_rangeName != (rangeName ? rangeName : "")) {
        cached_integral = RooHistPdf::analyticalIntegral(code, rangeName);
        cached_code = code;
        cached_rangeName = rangeName ? rangeName : "";
    }
    return cached_integral;
}
