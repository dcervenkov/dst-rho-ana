/**
 *  @file    dsrhobackground.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2017-12-15
 *
 *  @brief Main file
 *
 */

#include "dsrhobackground.h"

// Standard includes
#include <stdio.h>

// Local includes
#include "colors.h"
#include "constants.h"
#include "fitterbkg.h"
#include "tools.h"

int main(int argc, char* argv[]) {
    if (argc != 3) {
        printf("ERROR: Wrong number of arguments.\n");
        printf(
            "Usage: %s [OPTION]... INPUT-FILE OUTPUT-DIR\n",
            argv[0]);
        return 2;
    }

    /// This is so I have to change only the next block if I change the
    /// ordering, etc. of arguments
    const char* file_path = argv[1];
    const char* results_path = argv[2];

    tools::SetupPlotStyle();
    colors::setColors();

    FitterBKG fitter;
    fitter.ReadInFile(file_path);
    fitter.SetPlotDir(results_path);

    fitter.Fit(&fitter.scf_dt_model_, fitter.dataset_);
    fitter.PlotWithPull(fitter.dt_, fitter.dataset_, &fitter.scf_dt_model_);

    fitter.Fit(&fitter.scf_phit_model_, fitter.dataset_);
    fitter.PlotWithPull(fitter.phit_, fitter.dataset_, &fitter.scf_phit_model_);

    fitter.Fit(&fitter.scf_thetat_model_, fitter.dataset_);
    fitter.PlotWithPull(fitter.thetat_, fitter.dataset_, &fitter.scf_thetat_model_);

    fitter.Fit(&fitter.scf_thetab_model_, fitter.dataset_);
    fitter.PlotWithPull(fitter.thetab_, fitter.dataset_, &fitter.scf_thetab_model_);

    return 0;
}
