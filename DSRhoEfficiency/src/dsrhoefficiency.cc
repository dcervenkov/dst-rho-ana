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
#include <getopt.h>
#include <stdio.h>
#include <array>
#include <sstream>

// Local includes
#include "colors.h"
#include "constants.h"
#include "fitter.h"
#include "tools.h"

int main(int argc, char* argv[]) {
    char** optionless_argv = NULL;
    fitter_options options;
    const int optionless_argc = ProcessCmdLineOptions(argc, argv, optionless_argv, options);

    if (optionless_argc != 4) {
        printf("ERROR: Wrong number of arguments.\n");
        printf("Usage: %s [OPTION]... EVTGEN-FILE GSIM-FILE PLOT-DIR\n", optionless_argv[0]);
        return 2;
    }

    const char* evtgen_filepath = optionless_argv[1];
    const char* gsim_filepath = optionless_argv[2];
    const char* output_dir = optionless_argv[3];

    tools::SetupPlotStyle();
    colors::setColors();

    Fitter fitter(evtgen_filepath, gsim_filepath, output_dir);

    fitter.thetat_.setBins(50);
    fitter.thetab_.setBins(50);
    fitter.phit_.setBins(50);

    // If no efficiency model is specified, use the KDE one
    if (!options.efficiency_model) {
        options.efficiency_model = 5;
    }

    if (options.efficiency_model == 5) {

        if (!options.efficiency_file) {
            options.efficiency_file = "efficiency";
        }
        std::array<double, 6> bin_kde_pars = {50, 50, 50, 0.2, 0.2, 0.4};
        std::array<double, 6> ada_kde_pars = {50, 50, 50, 0.1, 0.1, 0.2};
        if (options.bin_kde_pars) {
            bin_kde_pars = DecodeStringPars(options.bin_kde_pars);
        }
        if (options.ada_kde_pars) {
            ada_kde_pars = DecodeStringPars(options.ada_kde_pars);
        }
        fitter.ProcessKDEEfficiency(options.efficiency_file, bin_kde_pars, ada_kde_pars,
                                    options.mirror_margin);

    } else if (options.efficiency_model == 6) {

        if (!options.efficiency_file) {
            options.efficiency_file = "efficiency.root";
        }
        std::array<double, 6> bin_kde_pars = {50, 50, 50, 0.2, 0.2, 0.4};
        std::array<double, 6> ada_kde_pars = {50, 50, 50, 0.1, 0.1, 0.2};
        if (options.bin_kde_pars) {
            bin_kde_pars = DecodeStringPars(options.bin_kde_pars);
        }
        if (options.ada_kde_pars) {
            ada_kde_pars = DecodeStringPars(options.ada_kde_pars);
        }
        fitter.ProcessKDEEfficiency2(options.efficiency_file, bin_kde_pars, ada_kde_pars,
                                    options.mirror_margin);

    } else if (options.efficiency_model == 7) {

        if (!options.efficiency_file) {
            options.efficiency_file = "efficiency.root";
        }
        fitter.ProcessNormalizedEfficiency(options.efficiency_file);

    } else {

        fitter.SetEfficiencyModel(options.efficiency_model);

        fitter.PlotVar(fitter.thetat_);
        fitter.FitEfficiency(fitter.thetat_);
        fitter.PlotEfficiency(fitter.thetat_, true, false);

        fitter.PlotVar(fitter.thetab_);
        fitter.FitEfficiency(fitter.thetab_);
        fitter.PlotEfficiency(fitter.thetab_, true, false);

        fitter.PlotVar(fitter.phit_);
        fitter.FitEfficiency(fitter.phit_);
        fitter.PlotEfficiency(fitter.phit_, true, false);

        fitter.PlotVars2D(fitter.thetat_, fitter.thetab_);
        fitter.PlotVars2D(fitter.thetat_, fitter.phit_);
        fitter.PlotVars2D(fitter.thetab_, fitter.phit_);

        fitter.PlotEfficiency2D(fitter.thetat_, fitter.thetab_);
        fitter.PlotEfficiency2D(fitter.thetat_, fitter.phit_);
        fitter.PlotEfficiency2D(fitter.thetab_, fitter.phit_);
    }

    return 0;
}

std::array<double, 6> DecodeStringPars(const char* string) {
    std::vector<double> vect;
    std::stringstream ss(string);

    double i;

    while (ss >> i) {
        vect.push_back(i);

        if (ss.peek() == ',' || ss.peek() == ' ') {
            ss.ignore();
        }
    }

    std::array<double, 6> arr;
    for (int i = 0; i < 6; i++) {
        arr[i] = vect[i];
    }
    return arr;
}

/*
 * Parses command line input and extracts switches and options from it, e.g.,
 * -h or --help. Then it acts accordingly, e.g., displaying help or setting
 * variables in an option struct. It also returns optionless_argv and
 * optionless_argc (return value) for easy integration with existing code.
 *
 * CAVEAT:
 * In order to pass negative numbers as arguments, one has to use the POSIX
 * "--" end of options indicator.
 *
 * @param argc Standard argc
 * @param argv Standard argv
 * @param optionless_argv Pointer where to write the new argv with processed switches removed
 * @param options Struct which holds the variables acted upon by switches
 */
int ProcessCmdLineOptions(const int argc, char* const argv[], char**& optionless_argv,
                          fitter_options& options) {
    int c;
    struct option long_options[] = {{"ada-kde-pars", required_argument, 0, 'a'},
                                    {"bin-kde-pars", required_argument, 0, 'b'},
                                    {"efficiency-model", required_argument, 0, 'e'},
                                    {"efficiency-file", required_argument, 0, 'f'},
                                    {"mirror-margin", required_argument, 0, 'm'},
                                    {"help", no_argument, 0, 'h'},
                                    {NULL, no_argument, NULL, 0}};
    int option_index = 0;
    while ((c = getopt_long(argc, argv, "a:b:e:f:m:h", long_options, &option_index)) != -1) {
        switch (c) {
            case 0:
                printf("option %s", long_options[option_index].name);
                if (optarg) printf(" with arg %s", optarg);
                printf("\n");
                break;
            case 'a':
                options.ada_kde_pars = optarg;
                break;
            case 'b':
                options.bin_kde_pars = optarg;
                break;
            case 'e':
                options.efficiency_model = atoi(optarg);
                break;
            case 'f':
                options.efficiency_file = optarg;
                break;
            case 'm':
                options.mirror_margin = atof(optarg);
                break;
            case 'h':
                printf("Usage: %s [OPTION]... INPUT-FILE OUTPUT_DIR\n\n", argv[0]);
                printf(
                    "Mandatory arguments to long options are mandatory for short options too.\n");
                printf(
                    "-a, --ada-kde-pars=\"BINS_X,BINS_Y,BINS_Z,WIDTH_X,WIDTH_Y,WIDTH_Z\" "
                    "parameters for adaptive KDE\n");
                printf(
                    "-b, --bin-kde-pars=\"BINS_X,BINS_Y,BINS_Z,WIDTH_X,WIDTH_Y,WIDTH_Z\" "
                    "parameters for binned KDE\n");
                printf(
                    "-f, --efficiency-file=FILE_NAME  filename to which to save KDE efficiency "
                    "map\n");
                printf("-e, --efficiency-model=MODEL_NUM efficiency model to use\n");
                printf(
                    "-m, --mirror-margin=MARGIN       fraction of phasespace to mirror to each "
                    "side in each dimension\n");
                printf("-h, --help                       display this text and exit\n");
                exit(0);
                break;
            default:
                printf("?? getopt returned character code 0%o ??\n", c);
        }
    }

    // Create a char** that will become the new argv, with the options removed
    const int optionless_argc = argc - optind + 1;
    optionless_argv = new char*[optionless_argc];
    // We want to keep the program name argument
    optionless_argv[0] = argv[0];
    for (int i = 1; i < optionless_argc; i++) {
        optionless_argv[i] = argv[i - 1 + optind];
    }

    return optionless_argc;
}