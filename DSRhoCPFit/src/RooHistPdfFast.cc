#include "RooHistPdfFast.h"

// Local includes
#include "log.h"

RooHistPdfFast::RooHistPdfFast(const char* name, const char* title, const RooArgSet& vars,
                               const RooDataHist& dhist, Int_t intOrder)
    : RooHistPdf(name, title, vars, dhist, intOrder) {
    Log::print(Log::debug, "In RooHistPdfFast constructor\n");
}

RooHistPdfFast::RooHistPdfFast(const RooHistPdfFast& other, const char* name) :
  RooHistPdf(other,name),
  cached_integral(other.cached_integral)
{
    Log::print(Log::debug, "In RooHistPdfFast copy constructor\n");
}

Double_t RooHistPdfFast::analyticalIntegral(Int_t code, const char* rangeName) const {
    if (cached_integral == 0) {
        cached_integral = RooHistPdf::analyticalIntegral(code);
    }
    // assert(code == 14);
    // Log::print(Log::debug, "analyticalIntegralCode = %i\n", code);
    // Log::print(Log::debug, "analyticalIntegral = %f\n", cached_integral);
    // return RooHistPdf::analyticalIntegral(code);
    return cached_integral;
}

// Double_t RooHistPdfFast::analyticalIntegralWN(Int_t code, const RooArgSet* normSet,
//                                               const char* rangeName) const {
//     Log::print(Log::debug, "In analyticalIntegralWN()\n");
//     return 1;
// }

// Double_t RooHistPdfFast::evaluate() const {
//     Log::print(Log::debug, "In evaluate()\n");
//     return 1;
// }

// Double_t RooHistPdfFast::getValV(const RooArgSet* nset) const {
//     Log::print(Log::debug, "In getValV()\n");
//     return 1;
// }
