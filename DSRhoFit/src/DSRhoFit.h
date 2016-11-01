#include "Fitter.h"
#include "FitterTDep.h"

int ProcessTrans(FitterTDep* fitter, Int_t doFit, Int_t doPlot);
int ProcessTrans(Fitter* fitter, Int_t doFit, Int_t doPlot);
void ConvertBetweenHelAndTrans(Double_t* par_input);
Double_t Round(Double_t number, Int_t digits);
