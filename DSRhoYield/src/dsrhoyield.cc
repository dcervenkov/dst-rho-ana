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
#include <getopt.h>
#include <iostream>

// ROOT includes
#include "TApplication.h"
#include "TCanvas.h"
#include "TEnv.h"

// Local includes
#include "colors.h"
#include "fitter.h"
#include "gitversion.h"
#include "log.h"
#include "tools.h"

int main(int argc, char* argv[]) {
    if (argc < 4) {
        printf("ERROR: Wrong number of arguments.\n");
        printf("Usage: DSRhoPeek [OPTION]... TRAINING-DIR DATA-DIR OUTPUT-DIR\n");
        return 2;
    }

    char** optionless_argv = nullptr;
    fit_options options;
    ProcessCmdLineOptions(argc, argv, optionless_argv, options);

    const char* trainingDir = optionless_argv[1];
    const char* dataDir = optionless_argv[2];
    const char* outputDir = optionless_argv[3];

    Log::setLogLevel(Log::debug);

    tools::SetPlotDir(outputDir);
    tools::SetupPlotStyle();
    colors::setColors();

    // The various parameters are first fixed from training data (6 MC streams)
    TChain* training_data = tools::ReadDataFromDir(trainingDir);
    Fitter fitter(outputDir);
    fitter.SetNumCPUs(options.cpus);
    fitter.SetMC(options.mc);

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

    fitter.PrintResults();

    //	fitter.SPlotFull(data);
    //	fitter.SPlotSB(data);
    //	fitter.SPlotSC(data);
    //	fitter.SPlotCB(data);
    //
	//	fitter.Setup(Components::all);
	//	Log::print(Log::info, "Correlation - all: %f\n", fitter.GetCorrelation(training_data, fitter.de_, fitter.thetab_, true));
	//	fitter.Setup(Components::signal);
	//	Log::print(Log::info, "Correlation - signal: %f\n", fitter.GetCorrelation(training_data, fitter.de_, fitter.thetab_, true));
	//	fitter.Setup(Components::crossfeed);
	//	Log::print(Log::info, "Correlation - crossfeed: %f\n", fitter.GetCorrelation(training_data, fitter.de_, fitter.thetab_, true));
	//	fitter.Setup(Components::background);
	//	Log::print(Log::info, "Correlation - background: %f\n", fitter.GetCorrelation(training_data, fitter.de_, fitter.thetab_, true));

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
                          fit_options& options) {
    int c;
    struct option long_options[] = {{"cpus", required_argument, 0, 'c'},
                                    {"MC", no_argument, 0, 'm'},
                                    {"version", no_argument, 0, 'v'},
                                    {"help", no_argument, 0, 'h'},
                                    {nullptr, no_argument, nullptr, 0}};
    int option_index = 0;
    while ((c = getopt_long(argc, argv, "c:mvh", long_options,
                            &option_index)) != -1) {
        switch (c) {
            case 0:
                printf("option %s", long_options[option_index].name);
                if (optarg) printf(" with arg %s", optarg);
                printf("\n");
                break;
            case 'c':
                options.cpus = atoi(optarg);
                break;
            case 'm':
                options.mc = true;
                break;
            case 'v':
                printf("Version: %s\n", gitversion);
                exit(0);
            case 'h':
                printf("Usage: %s [OPTION]... TRAINING-DIR DATA-DIR OUTPUT-DIR\n\n", argv[0]);
                printf("Mandatory arguments to long options are mandatory for short options too.\n");
                printf("-c, --cpus=NUM_CPUS              number of CPU cores to use for fitting and plotting\n");
                printf("-h, --help                       display this text and exit\n");
                printf("-m, --MC                         final fit is on MC (default: data)\n");
                printf("-o, --output                     basename of the output files\n");
                printf("-p, --plot-dir=PLOT-DIR          directory for plots\n");
                printf("-v, --version                	 display version and exit\n");
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
