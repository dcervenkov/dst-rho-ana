/**
 *  @file    dsrhoyield.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2015-07-16
 *
 *  @brief Main file
 *
 */

#include "dsrhoyield.h"

// Standard includes
#include <iostream>

// ROOT includes
#include "TApplication.h"
#include "TCanvas.h"
#include "TEnv.h"

// Local includes
#include "fitter.h"
#include "tools.h"
#include "colors.h"

//#define GRAPHICS

int main(int argc, char* argv[]) {
#ifdef GRAPHICS
	TApplication* rootapp = new TApplication("DSRhoYield", &argc, argv);
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
#endif // GRAPHICS

	if (argc != 4){
		printf("ERROR: Wrong number of arguments.\n");
		printf("Usage: DSRhoPeek TRAINING-DIR DATA-DIR OUTPUT-DIR\n");
		return 2;
	}

	const char* trainingDir = argv[1];
	const char* dataDir = argv[2];
	const char* outputDir = argv[3];

	tools::SetupPlotStyle();
	colors::setColors();

	gEnv->SetValue("Canvas.PrintDirectory",outputDir);
	printf("print dir: %s\n",gEnv->GetValue("Canvas.PrintDirectory","not found"));

	// The various parameters are first fixed from training data (6 MC streams)
	TChain* training_data = tools::ReadDataFromDir(trainingDir);
	Fitter fitter(outputDir);

	fitter.Setup(Components::signal);
	fitter.FitTo(training_data);
	fitter.FixShape(Components::signal);
	fitter.WriteFitResults();
	fitter.Plot();

	fitter.Setup(Components::crossfeed);
	fitter.FitTo(training_data);
	fitter.FixShape(Components::crossfeed);
	fitter.WriteFitResults();
	fitter.Plot();

	fitter.Setup(Components::signal_plus_crossfeed);
	fitter.FitTo(training_data);
	fitter.FixShape(Components::signal_plus_crossfeed);
	fitter.WriteFitResults();
	fitter.Plot();

	fitter.Setup(Components::background);
	fitter.FitTo(training_data);
	fitter.FixShape(Components::background);
	fitter.WriteFitResults();
	fitter.Plot();

	// The final fit is done on a different dataset
	TChain* data = tools::ReadDataFromDir(dataDir);
	fitter.Setup(Components::all);
	fitter.FitTo(data);
	fitter.WriteFitResults();
	fitter.Plot();

//	fitter.SPlotFull(data);
//	fitter.SPlotSB(data);
//	fitter.SPlotSC(data);
//	fitter.SPlotCB(data);
//
//	fitter.Setup(Components::all);
//	printf("Correlation - all: %f\n", fitter.GetCorrelation(training_data, fitter.de_, fitter.thetab_, true));
//	fitter.Setup(Components::signal);
//	printf("Correlation - signal: %f\n", fitter.GetCorrelation(training_data, fitter.de_, fitter.thetab_, true));
//	fitter.Setup(Components::crossfeed);
//	printf("Correlation - crossfeed: %f\n", fitter.GetCorrelation(training_data, fitter.de_, fitter.thetab_, true));
//	fitter.Setup(Components::background);
//	printf("Correlation - background: %f\n", fitter.GetCorrelation(training_data, fitter.de_, fitter.thetab_, true));


#ifdef GRAPHICS
	// Write the file to disk so the ROOT file is complete, but don't close the TApp windows
	fitter.CloseOutput();
	rootapp->Run();
#endif // GRAPHICS
	return 0;
}







