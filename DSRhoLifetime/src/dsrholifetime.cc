/**
 *  @file    dsrholifetime.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2016-04-05
 *
 *  @brief Main file
 *
 */

#include "dsrholifetime.h"

// Standard includes
#include <getopt.h>
#include <stdio.h>

// Local includes
#include "colors.h"
#include "config.h"
#include "constants.h"
#include "fitterlifetime.h"
#include "gitversion.h"
#include "tools.h"

int main(int argc, char* argv[]) {
    char** optionless_argv = nullptr;
    // The {} causes the struct's members to be initialized to 0. Without it
    // they would have unspecified values
    fitter_options options = {};
    const int optionless_argc = ProcessCmdLineOptions(argc, argv, optionless_argv, options);

    if (optionless_argc < 3) {
        printf("ERROR: Wrong number of arguments.\n");
        printf("Usage: %s [OPTION]... RESULTS_FILE INPUT-FILE(S)...\n", optionless_argv[0]);
        return 2;
    }

    const char* results_file = optionless_argv[1];

    std::vector<const char*> file_names;
    for (int i = 2; i < optionless_argc; i++) {
        file_names.push_back(optionless_argv[i]);
    }

    tools::SetupPlotStyle();
    colors::setColors();

    Config config;
    config.ReadInJSONFile(options.config);

    FitterLifetime fitter(config.json);

    if (options.num_CPUs_set) fitter.SetNumCPUs(options.num_CPUs);
    if (options.plot_dir_set) {
        tools::SetPlotDir(options.plot_dir);
        fitter.SetMakePlots(true);
    }
    if (options.do_lifetime_fit_set) fitter.SetDoLifetimeFit(options.do_lifetime_fit);
    if (options.do_mixing_fit_set) fitter.SetDoMixingFit(options.do_mixing_fit);
    if (options.perfect_tagging_set) fitter.SetPerfectTagging(options.perfect_tagging);
    if (options.channel_set) fitter.SetChannel(options.channel);
    if (options.physical_pdf_set) fitter.SetUsePhysicalPdf(options.physical_pdf);
    if (options.components_set) fitter.SetComponents(options.components);
    if (options.num_events_set) {
        fitter.ReadInFile(file_names, options.num_events);
    } else {
        fitter.ReadInFile(file_names);
    }

    fitter.Process();

    fitter.SaveTXTResults(results_file);

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
        {"config", required_argument, 0, 'g'},
        {"channel", required_argument, 0, 'a'},
        {"cpus", required_argument, 0, 'c'},
        {"components", required_argument, 0, 'o'},
        {"events", required_argument, 0, 'e'},
        {"plot-dir", required_argument, 0, 'p'},
        {"lifetime", no_argument, 0, 'l'},
        {"mixing", no_argument, 0, 'm'},
        {"perfecttag", no_argument, 0, 't'},
        {"physics", no_argument, 0, 'y'},
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {nullptr, no_argument, nullptr, 0}};

    int option_index = 0;
    while ((c = getopt_long(argc, argv, "a:c:e:g:o:plmtyvh", long_options, &option_index)) != -1) {
        switch (c) {
            case 0:
                printf("option %s", long_options[option_index].name);
                if (optarg) printf(" with arg %s", optarg);
                printf("\n");
                break;
            case 'a':
                options.channel = optarg;
                options.channel_set = true;
                break;
            case 'c':
                options.num_CPUs = atoi(optarg);
                options.num_CPUs_set = true;
                break;
            case 'e':
                options.num_events = atoi(optarg);
                options.num_events_set = true;
                break;
            case 'g':
                options.config = optarg;
                options.config_set = true;
                break;
            case 'o':
                options.components = optarg;
                options.components_set = true;
                break;
            case 'p':
                options.plot_dir = optarg;
                options.plot_dir_set = true;
                break;
            case 'l':
                options.do_lifetime_fit = true;
                options.do_lifetime_fit_set = true;
                break;
            case 'm':
                options.do_mixing_fit = true;
                options.do_mixing_fit_set = true;
                break;
            case 't':
                options.perfect_tagging = true;
                options.perfect_tagging_set = true;
                break;
            case 'y':
                options.physical_pdf = true;
                options.physical_pdf_set = true;
                break;
            case 'v':
                printf("Version: %s\n", gitversion);
                exit(0);
                break;
            case 'h':
                printf("Usage: %s [OPTION]... RESULTS_FILE INPUT-FILE(S)...\n\n", argv[0]);
                printf(
                    "Mandatory arguments to long options are mandatory for short options too.\n");
                printf("-a, --channel=CHANNEL        channel whose SCF and BKG pars to load\n");
                printf("-c, --cpus=NUM_CPUS          number of CPU cores to use for fitting and plotting\n");
                printf("-e, --events=NUM_EVENTS      number of events to be imported from the input file\n");
                printf("-g, --config=CONFIG-FILE     read in configuration from the specified file\n");
                printf("-h, --help                   display this text and exit\n");
                printf("-l, --lifetime               make a lifetime fit\n");
                printf("-m, --mixing                 make a mixing fit\n");
                printf("-o, --components=COMPONENTS  components which to read and fit (CR, CRSCF, all)\n");
                printf("-p, --plot-dir=DIR           create lifetime/mixing plots and save to DIR\n");
                printf("-t, --perfecttag             use MC info to get perfect tagging\n");
                printf("-y, --physics                use physics-based PDF for SCF and BKG\n");
                printf("-v, --version                print program version and exit\n");
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
