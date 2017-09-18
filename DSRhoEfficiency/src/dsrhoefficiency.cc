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
#include "colors.h"
#include "constants.h"
#include "fitter.h"
#include "tools.h"

int main(int argc, char* argv[]) {
    if (argc != 5) {
        printf("ERROR: Wrong number of arguments.\n");
        printf("Usage: DSRhoEfficiency MODEL-NUM EVTGEN-FILE GSIM-FILE OUTPUT-DIR\n");
        return 2;
    }

    const int model_num = strtol(argv[1], NULL, 10);
    const char* evtgen_filepath = argv[2];
    const char* gsim_filepath = argv[3];
    const char* output_dir = argv[4];

    tools::SetupPlotStyle();
    colors::setColors();

    Fitter fitter(evtgen_filepath, gsim_filepath, output_dir);

    fitter.thetat_.setBins(50);
    fitter.thetab_.setBins(50);
    fitter.phit_.setBins(50);

    // fitter.Process1DKDEEfficiency();
    fitter.ProcessKDEEfficiency();
    // fitter.ProcessBinnedEfficiency();

    // fitter.SetEfficiencyModel(model_num);

    // fitter.PlotVar(fitter.thetat_);
    // fitter.FitEfficiency(fitter.thetat_);
    // fitter.PlotEfficiency(fitter.thetat_, true, false);

    // fitter.PlotVar(fitter.thetab_);
    // fitter.FitEfficiency(fitter.thetab_);
    // fitter.PlotEfficiency(fitter.thetab_, true, false);

    // fitter.PlotVar(fitter.phit_);
    // fitter.FitEfficiency(fitter.phit_);
    // fitter.PlotEfficiency(fitter.phit_, true, false);

    // fitter.PlotVars2D(fitter.thetat_, fitter.thetab_);
    // fitter.PlotVars2D(fitter.thetat_, fitter.phit_);
    // fitter.PlotVars2D(fitter.thetab_, fitter.phit_);

    // fitter.PlotEfficiency2D(fitter.thetat_, fitter.thetab_);
    // fitter.PlotEfficiency2D(fitter.thetat_, fitter.phit_);
    // fitter.PlotEfficiency2D(fitter.thetab_, fitter.phit_);

    return 0;
}
