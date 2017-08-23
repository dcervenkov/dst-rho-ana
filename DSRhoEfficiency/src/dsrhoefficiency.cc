/**
 *  @file    dsrhoefficiency.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2015-11-24
 *
 *  @brief Main file
 *
 */

#include "dsrhoefficiency.h"

// Standard includes
#include <stdio.h>

// Local includes
#include "fitter.h"
#include "tools.h"
#include "colors.h"
#include "constants.h"
#include "tools.h"

int main(int argc, char* argv[]) {
	if (argc != 4){
		printf("ERROR: Wrong number of arguments.\n");
		printf("Usage: DSRhoEfficiency EVTGEN-FILE GSIM-FILE OUTPUT-DIR\n");
		return 2;
	}

	const char* evtgen_filepath = argv[1];
	const char* gsim_filepath = argv[2];
	const char* output_dir = argv[3];

	tools::SetupPlotStyle();
	colors::setColors();

	Fitter fitter(evtgen_filepath, gsim_filepath, output_dir);

	fitter.thetat_.setBins(40);
	fitter.thetab_.setBins(40);
	fitter.phit_.setBins(40);

	fitter.PlotVar(fitter.dt_);

	fitter.PlotVar(fitter.thetat_);
	fitter.FitEfficiency(fitter.thetat_);
	fitter.PlotEfficiency(fitter.thetat_, true, false);

	fitter.PlotVar(fitter.thetab_);
	fitter.FitEfficiency(fitter.thetab_);
	fitter.PlotEfficiency(fitter.thetab_, true, false);

	fitter.PlotVar(fitter.phit_);
	fitter.FitEfficiency(fitter.phit_);
	fitter.PlotEfficiency(fitter.phit_, true, false);

	fitter.PlotVars2D(fitter.thetab_, fitter.thetat_);
	fitter.PlotVars2D(fitter.thetat_, fitter.phit_);
	fitter.PlotVars2D(fitter.thetab_, fitter.phit_);

	fitter.PlotVars2D(fitter.dt_, fitter.thetat_);
	fitter.PlotVars2D(fitter.dt_, fitter.thetab_);
	fitter.PlotVars2D(fitter.dt_, fitter.phit_);

	fitter.PlotEfficiency2D(fitter.thetab_, fitter.thetat_);
	fitter.PlotEfficiency2D(fitter.thetat_, fitter.phit_);
	fitter.PlotEfficiency2D(fitter.thetab_, fitter.phit_);

	fitter.PlotEfficiency2D(fitter.dt_, fitter.thetat_);
	fitter.PlotEfficiency2D(fitter.dt_, fitter.thetab_);
	fitter.PlotEfficiency2D(fitter.dt_, fitter.phit_);

	return 0;
}







