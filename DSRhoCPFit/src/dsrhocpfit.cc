/**
 *  @file    dsrholifetime.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2016-04-05
 *
 *  @brief Main file
 *
 */

#include "dsrhocpfit.h"

// Standard includes
#include <getopt.h>
#include <stdio.h>

// Local includes
#include "colors.h"
#include "constants.h"
#include "fittercpv.h"
#include "tools.h"

int main(int argc, char* argv[]) {
    char** optionless_argv = NULL;
    // The {} causes the struct's members to be initialized to 0. Without it
    // they would have unspecified values
    fitter_options options = {};
    const int optionless_argc = ProcessCmdLineOptions(argc, argv, optionless_argv, options);

    if (optionless_argc != 19) {
        printf("ERROR: Wrong number of arguments.\n");
        printf(
            "Usage: %s [OPTION]... -- INPUT-FILE OUTPUT_DIR "
            "AP APA A0 ATA XP X0 XT YP Y0 YT XPB X0B XTB YPB Y0B YTB\n",
            optionless_argv[0]);
        return 2;
    }

    /// This is so I have to change only the next block if I change the
    /// ordering, etc. of arguments
    const char* file_path = optionless_argv[1];
    const char* output_dir = optionless_argv[2];
    const int numPars = optionless_argc - 3;
    double par_input[numPars];
    for (Int_t i = 0; i < numPars; i++) {
        par_input[i] = atof(optionless_argv[i + 3]);
    }

    tools::SetupPlotStyle();
    colors::setColors();

    FitterCPV fitter(par_input);

    if (options.num_CPUs_set) fitter.SetNumCPUs(options.num_CPUs);
    if (options.make_plots_set) fitter.SetMakePlots(options.make_plots);
    if (options.do_mixing_fit_set) fitter.SetDoMixingFit(options.do_mixing_fit);
    if (options.perfect_tagging_set) fitter.SetPerfectTagging(options.perfect_tagging);
    if (options.fix_set) {
        if (fitter.FixParameters(options.fix)) {
            return 1;
        }
    }

    if (options.num_events_set) {
        fitter.ReadInFile(file_path, options.num_events);
    } else {
        fitter.ReadInFile(file_path);
    }
    fitter.SetOutputDir(output_dir);

    fitter.Test();
    //  fitter.GenerateToys(10000, 10);
    //  fitter.Test();

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
    struct option long_options[] = {
        {"cpus", required_argument, 0, 'c'},
        {"events", required_argument, 0, 'e'},
        {"fix", required_argument, 0, 'x'},
        {"mixing", no_argument, 0, 'm'},
        {"perfecttag", no_argument, 0, 't'},
        {"plot", no_argument, 0, 'p'},
        {"help", no_argument, 0, 'h'},
        {NULL, no_argument, NULL, 0}};
    int option_index = 0;
    while ((c = getopt_long(argc, argv, "c:e:x:plmth", long_options, &option_index)) != -1) {
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
            case 'e':
                options.num_events = atoi(optarg);
                options.num_events_set = true;
                break;
            case 'x':
                options.fix = optarg;
                options.fix_set = true;
                break;
            case 'p':
                options.make_plots = true;
                options.make_plots_set = true;
                break;
            case 'm':
                options.do_mixing_fit = true;
                options.do_mixing_fit_set = true;
                break;
            case 't':
                options.perfect_tagging = true;
                options.perfect_tagging_set = true;
                break;
            case 'h':
                printf("Usage: %s [OPTION]... INPUT-FILE OUTPUT_DIR\n\n", argv[0]);
                printf("Mandatory arguments to long options are mandatory for short options too.\n");
                printf("-c, --cpus=NUM_CPUS     number of CPU cores to use for fitting and plotting\n");
                printf("-e, --events=NUM_EVENTS number of events to be imported from the input file\n");
                printf("-x, --fix=ARG1,ARG2,... fix specified argument(s) to input values in the fit\n");
                printf("-h, --help              display this text and exit\n");
                printf("-m, --mixing            make a mixing fit\n");
                printf("-t, --perfecttag        use MC info to get perfect tagging\n");
                printf("-p, --plot              create lifetime/mixing plots\n");
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
