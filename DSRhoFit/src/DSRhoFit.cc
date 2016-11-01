#include "DSRhoFit.h"

// ROOT includes
#include "TStyle.h"
#include "TApplication.h"
#include "TStopwatch.h"
#include "TComplex.h"

// Local includes
#include "Constants.h"
#include "ASSERT.h"
#include "Fitter.h"
#include "FitterTDep.h"
#include "tools.h"
#include "colors.h"

#define DEBUG
//#define VERBOSE
//#define GRAPHIC

char* inputFile = 0;
char* outputFile = 0;

int main(int argc, char* argv[]) {
#ifdef GRAPHIC
	TApplication* rootapp = new TApplication("DSRhoSens",&argc,argv); //For graphic-output apps only
	/**
	 * When using TApplication it removes arguments it "handles" from
	 * the argument array. E.g. -b, -x, -q, --help, <dir>, <file>, etc.
	 * For more info read TApplication's GetOptions function help.
	 * The solution is to use rootapp->Argc() and rootapp->Argv(i).
	 * The next few lines are for compatibility of GRAPHIC vs. non-GRAPHIC -
	 * they recreate the original argc and argv when using GRAPHIC.
	 **/
	argc = rootapp->Argc();
	for(int i = 0; i < argc; i++) {
		argv[i] = rootapp->Argv(i);
	}
#endif

	if (argc != 21 && argc != 9) {
		printf("ERROR: Wrong number of arguments.\n");
		printf(
				"Usage: DSRhoFit inputFile outputFile ap apa a0 ata xp x0 xt yp y0 yt xpb x0b xtb ypb y0b ytb doFit doPlot\n");
		return 85;
	}

	tools::SetupPlotStyle();
	colors::setColors();
	gStyle->SetOptFit(1);

	TStopwatch timer;
	timer.Start();

	/// This is so I have to change only the next block if I change the
	/// ordering, etc. of arguments
	inputFile = argv[1];
	outputFile = argv[2];
	Int_t numPars = argc - 5;
	Double_t par_input[numPars];
	for (Int_t i = 0; i < numPars; i++)
		par_input[i] = atof(argv[i + 3]);

	Int_t doFit = atoi(argv[argc - 2]);
	Int_t doPlot = atoi(argv[argc - 1]);

	if (doFit == 3) {
		ConvertBetweenHelAndTrans(par_input);
//	} else if (doFit == 4) {
//		FitterTDep* fitter = new FitterTDep(par_input);
//		fitter->GenerateDataSet(100000);
//		fitter->GetDataSet()->write(inputFile);
//	} else if (doFit == 1) {
//		//ConvertBetweenHelAndTrans(par_input);
//		FitterTDep* fitter = new FitterTDep(par_input);
//		fitter->ReadDataSet(inputFile);
//		ProcessTrans(fitter, doFit, doPlot);
	} else if (doFit == 2) {
		//ConvertBetweenHelAndTrans(par_input);
		Fitter* fitter = new Fitter(par_input);
		fitter->ReadDataSet(inputFile);
		ProcessTrans(fitter, doFit, doPlot);
	} else if (doFit == 0) {
		//ConvertBetweenHelAndTrans(par_input);
		Fitter* fitter = new Fitter(par_input);
		fitter->ReadDataSet(inputFile);
		ProcessTrans(fitter, doFit, doPlot);
	}

	timer.Stop();
	timer.Print();

#ifdef GRAPHIC
	printf("\nProgram execution has finished.\n");
	rootapp->Run();
#endif
	return 0;
}

int ProcessTrans(FitterTDep* fitter, Int_t doFit, Int_t doPlot) {
	if (doFit) {
		fitter->FixAllParameters();
		fitter->FreeParameter("ap");
		fitter->FreeParameter("apa");
		fitter->FreeParameter("a0");
		fitter->FreeParameter("ata");
//        fitter->FreeParameter("xp");
//        fitter->FreeParameter("x0");
//        fitter->FreeParameter("xt");
//        fitter->FreeParameter("yp");
//        fitter->FreeParameter("y0");
//        fitter->FreeParameter("yt");
//        fitter->FreeParameter("xpb");
//        fitter->FreeParameter("x0b");
//        fitter->FreeParameter("xtb");
//        fitter->FreeParameter("ypb");
//        fitter->FreeParameter("y0b");
//        fitter->FreeParameter("ytb");
		fitter->Fit();
		fitter->SaveParameters(outputFile);
	}

	if (doPlot == kTRUE) {
		//SaveChi2Maps(fitter->GetBinnedDataSet(),dataSet->numEntries(),fitter->GetPdf(),*(fitter->GetTht()),*(fitter->GetThb()),*(fitter->GetPhit()));
		//Double_t mychi2 = fitter->SaveChi2Maps("a");
		//fitter->SaveResiduals();
		//fitter->SaveNllPlot("yt");
		//printf("mychi2 from SaveChi2Maps = %f\n",mychi2);
//		SavePlots(fitter->GetDataSet(), fitter->GetPdf(), *(fitter->GetTht()), *(fitter->GetThb()),
//				*(fitter->GetPhit()), *(fitter->GetDt()));
	}

	return 0;
}

int ProcessTrans(Fitter* fitter, Int_t doFit, Int_t doPlot) {
	if (doFit) {
		fitter->FixAllParameters();
		fitter->FreeParameter("ap");
		fitter->FreeParameter("apa");
		fitter->FreeParameter("a0");
		fitter->FreeParameter("ata");
		fitter->Fit();
		fitter->SaveParameters(outputFile);
//		fitter->SaveHelParameters(outputFile);
	}

	if (doPlot == kTRUE) {
		fitter->SaveVarPlots();
//		Double_t mychi2 = fitter->SaveChi2Maps();
//		fitter->SaveResiduals();
		//fitter->SaveNllPlot("yt");
//		printf("mychi2 from SaveChi2Maps = %f\n",mychi2);
	}

	return 0;
}

void ConvertBetweenHelAndTrans(Double_t* par_input) {
	/// The variables are named as if converting parameters from helicity to transversity
	/// but the transformation is symmetric and can be used to convert from transversity
	/// to helicity as well.

	TComplex hp(par_input[0], par_input[1], true);
	TComplex h0(par_input[2], 0, true);
	TComplex hm(sqrt(1 - hp.Rho2() - h0.Rho2()), par_input[3], true);

	TComplex hRhop(par_input[5], par_input[8], true);
	TComplex hRho0(par_input[6], par_input[9], true);
	TComplex hRhom(par_input[7], par_input[10], true);

	TComplex hps = hRhop * hp;
	TComplex h0s = hRho0 * h0;
	TComplex hms = hRhom * hm;

	TComplex ap = (hp + hm) / sqrt(2);
	TComplex a0 = h0;
	TComplex at = (hp - hm) / sqrt(2);

	TComplex aps = (hps + hms) / sqrt(2);
	TComplex a0s = h0s;
	TComplex ats = (hps - hms) / sqrt(2);

	TComplex tRhop = aps / ap;
	TComplex tRho0 = a0s / a0;
	TComplex tRhot = ats / at;

	printf("original hel:\t");
	printf(
			"hp = %.3f\thpa = %.2f\th0 = %.3f\thma = %.2f\tphiw = %.4f\trp = %.3f\tr0 = %.3f\trm = %.3f\tsp = %.3f\ts0 = %.3f\tsm = %.3f\n",
			par_input[0], par_input[1], par_input[2], par_input[3], par_input[4], par_input[5], par_input[6],
			par_input[7], par_input[8], par_input[9], par_input[10]);

	par_input[0] = ap.Rho();
	par_input[1] = ap.Theta();
	par_input[2] = a0.Rho();
	par_input[3] = at.Theta();
	par_input[5] = tRhop.Rho();
	par_input[6] = tRho0.Rho();
	par_input[7] = tRhot.Rho();
	par_input[8] = tRhop.Theta();
	par_input[9] = tRho0.Theta();
	par_input[10] = tRhot.Theta();

	printf("converted trans:");
	printf(
			"ap = %.3f\tapa = %.2f\ta0 = %.3f\tata = %.2f\tphiw = %.4f\trp = %.3f\tr0 = %.3f\trt = %.3f\tsp = %.3f\ts0 = %.3f\tst = %.3f\n",
			par_input[0], par_input[1], par_input[2], par_input[3], par_input[4], par_input[5], par_input[6],
			par_input[7], par_input[8], par_input[9], par_input[10]);

}

Double_t Round(Double_t number, Int_t digits) {
	number = number * pow(10, digits);
	if (fmod(number, 1) > 0.5)
		number = ceil(number);
	else
		number = floor(number);

	return number / pow(10, digits);
}
