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
#include <array>
#include <getopt.h>
#include <stdio.h>

// Local includes
#include "colors.h"
#include "constants.h"
#include "fitterbkg.h"
#include "tools.h"

int main(int argc, char* argv[]) {
    // The {} causes the struct's members to be initialized to 0. Without it
    // they would have unspecified values
    fitter_options options = {};

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
    fitter.PlotVar(*(fitter.phit_));

    return 0;
}
