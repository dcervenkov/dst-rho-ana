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
#include <getopt.h>
#include <stdio.h>
#include <vector>

#include "AdaptiveKernelDensity.hh"

// Local includes
#include "colors.h"
#include "constants.h"
#include "fitterbkg.h"
#include "log.h"
#include "tools.h"

int main(int argc, char* argv[]) {
    Log::setLogLevel(Log::debug);
    char** optionless_argv = NULL;
    fitter_options options = {};
    const int optionless_argc = ProcessCmdLineOptions(argc, argv, optionless_argv, options);

    if (optionless_argc < 2) {
        printf("ERROR: Wrong number of arguments.\n");
        printf("Usage: %s [OPTION]... INPUT-FILE(S)...\n", optionless_argv[0]);
        return 2;
    }

    std::vector<const char*> file_names;
    for (int i = 1; i < optionless_argc; i++) {
        file_names.push_back(optionless_argv[i]);
    }

    tools::SetupPlotStyle();
    colors::setColors();

    FitterBKG fitter;
    fitter.ReadInFile(file_names);
    if (options.plot_dir_set) {
        fitter.SetPlotDir(options.plot_dir);
    }

    if (options.num_CPUs_set) {
        fitter.SetNumCPUs(options.num_CPUs);
    }

    if (options.KDE && options.histo) {
        Log::print(Log::error, "Both '--kde' and '--histo' options specified. Use only one!\n");
        return 3;
    }

    if (options.KDE) {
        fitter.thetat_.setBins(50);
        fitter.thetab_.setBins(50);
        fitter.phit_.setBins(50);
        AdaptiveKernelDensity kde = fitter.FitKDE(fitter.dataset_);
        fitter.PlotKDE(kde);
    } else if (options.histo) {
        fitter.CreateHistoPDF(fitter.dataset_);
    } else {
    
        // We must plot after each fit otherwise the results printed on the plot
        // would be of the last fit
        fitter.Fit(&fitter.scf_phit_model_, fitter.dataset_);
        if (options.plot_dir_set) {
            fitter.PlotWithPull(fitter.phit_, *fitter.dataset_, fitter.scf_phit_model_);
        }

        fitter.Fit(&fitter.scf_thetat_model_, fitter.dataset_);
        if (options.plot_dir_set) {
            fitter.PlotWithPull(fitter.thetat_, *fitter.dataset_, fitter.scf_thetat_model_);
        }

        fitter.Fit(&fitter.scf_thetab_model_, fitter.dataset_);
        if (options.plot_dir_set) {
            fitter.PlotWithPull(fitter.thetab_, *fitter.dataset_, fitter.scf_thetab_model_);
        }

        fitter.Fit(&fitter.scf_dt_model_, fitter.dataset_);
        if (options.plot_dir_set) {
            fitter.PlotWithPull(fitter.dt_, *fitter.dataset_, fitter.scf_dt_model_);
        }

        fitter.Fit(&fitter.scf_dt_model_, fitter.dataset_a_);
        if (options.plot_dir_set) {
            fitter.PlotWithPull(fitter.dt_, *fitter.dataset_a_, fitter.scf_dt_model_);
        }

        fitter.Fit(&fitter.scf_dt_model_, fitter.dataset_ab_);
        if (options.plot_dir_set) {
            fitter.PlotWithPull(fitter.dt_, *fitter.dataset_ab_, fitter.scf_dt_model_);
        }

        fitter.Fit(&fitter.scf_dt_model_, fitter.dataset_b_);
        if (options.plot_dir_set) {
            fitter.PlotWithPull(fitter.dt_, *fitter.dataset_b_, fitter.scf_dt_model_);
        }

        fitter.Fit(&fitter.scf_dt_model_, fitter.dataset_bb_);
        if (options.plot_dir_set) {
            fitter.PlotWithPull(fitter.dt_, *fitter.dataset_bb_, fitter.scf_dt_model_);
        }
    }

    return 0;
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
    struct option long_options[] = {{"cpus", required_argument, 0, 'c'},
                                    {"kde", no_argument, 0, 'k'},
                                    {"histo", no_argument, 0, 's'},
                                    {"plot-dir", required_argument, 0, 'p'},
                                    {"help", no_argument, 0, 'h'},
                                    {NULL, no_argument, NULL, 0}};
    int option_index = 0;
    while ((c = getopt_long(argc, argv, "c:p:ksh", long_options, &option_index)) != -1) {
        switch (c) {
            case 0:
                printf("option %s", long_options[option_index].name);
                if (optarg) printf(" with arg %s", optarg);
                printf("\n");
                break;
            case 'c':
                options.num_CPUs = atoi(optarg);
                options.num_CPUs_set = true;
                break;
            case 'k':
                options.KDE = true;
                break;
            case 's':
                options.histo = true;
                break;
            case 'p':
                options.plot_dir = optarg;
                options.plot_dir_set = true;
                break;
            case 'h':
                printf("Usage: %s [OPTION]... INPUT-FILES\n\n", argv[0]);
                printf("Mandatory arguments to long options are mandatory for short options too.\n");
                printf("-c, --cpus=NUM_CPUS       number of CPU cores to use for fitting and plotting\n");
                printf("-k, --kde                 use kernel density estimation\n");
                printf("-s, --histo               create a histogram PDF\n");
                printf("-h, --help                display this text and exit\n");
                printf("-p, --plot-dir=PLOT_DIR   create lifetime/mixing plots\n");
                exit(0);
                break;
            default:
                printf("?? getopt returned character code 0%o ??\n", c);
                exit(1);
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
