/**
 *  @file    spherical_harmonics_approx.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2019-04-06
 *
 *  @brief Main file
 *
 */

#include "spherical_harmonics_approx.h"

// Standard includes
#include <array>
#include <cstdio>
#include <cstring>
#include <getopt.h>
#include <ostream>
#include <sstream>
#include <stdio.h>

// Boost includes
#include <boost/algorithm/string/predicate.hpp>

// ROOT includes
#include "TSystem.h"

// Local includes
#include "colors.h"
#include "constants.h"
#include "fittershe.h"
#include "gitversion.h"
#include "log.h"
#include "teebuf.h"
#include "tools.h"


int main(int argc, char* argv[]) {
    // This block sets up an alternative ostream that can be used to duplicate
    // cout to both the screen and a stringstream. We use this to save a copy of
    // the log to the resulting ROOT file.
    std::ostringstream sout;
    teebuf sbuf(sout.rdbuf(), std::cout.rdbuf());
    std::ostream out(&sbuf);
    std::streambuf* orig_cout_streambuf = std::cout.rdbuf();

    Log::setLogLevel(Log::debug);

    char** optionless_argv = nullptr;
    // The {} causes the struct's members to be initialized to 0. Without it
    // they would have unspecified values
    fitter_options options = {};
    const int optionless_argc = ProcessCmdLineOptions(argc, argv, optionless_argv, options);

    if (optionless_argc < 3) {
        printf("ERROR: Not enough arguments.\n");
        printf("Usage: %s [OPTION]... OUTPUT_FILE INPUT-FILE(s)\n",
               optionless_argv[0]);
        return 2;
    }

    if (boost::algorithm::ends_with(gitversion, "-dirty")) {
        Log::print(Log::warning, "Using version from dirty Git worktree\n");
    }

    const char* results_path = optionless_argv[1];

    std::vector<const char*> file_names;
    for (int i = 2; i < optionless_argc; i++) {
        file_names.push_back(optionless_argv[i]);
    }

    if (options.save_log_set && options.save_log) {
        std::cout.rdbuf(out.rdbuf());
    }

    tools::SetupPlotStyle();
    colors::setColors();

    FitterSHE fitter;

    fitter.SetOutputFile(results_path);

    if (options.num_events_set) {
        fitter.ReadInFile(file_names, options.num_events);
    } else {
        fitter.ReadInFile(file_names);
    }

    if (options.plot_dir_set) fitter.SetPlotDir(options.plot_dir);


    fitter.AnalyzeDataset("t", -1, 20, 50, 2, 1, 5, 5, "fast");


    fitter.LogCLIArguments(argc, argv);
    fitter.LogEnvironmentMetadata();
    for (auto file_name : file_names) {
        fitter.LogText("input_file_name", file_name);
        fitter.LogFileCRC("input_file_crc", file_name);
    }

    if (options.save_log_set && options.save_log) {
        // We have to restore cout's buffer to the original, otherwise we would
        // get a segfault as our object goes out of scope sooner than cout
        std::cout.rdbuf(orig_cout_streambuf);
        Log::print(Log::info, "Saving log to ROOT file\n");
        fitter.LogText("log", sout.str().c_str());
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
    struct option long_options[] = {
        {"events", required_argument, 0, 'n'},
        {"log", no_argument, 0, 'l'},
        {"plot-dir", required_argument, 0, 'p'},
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {nullptr, no_argument, nullptr, 0}};
    int option_index = 0;
    while ((c = getopt_long(argc, argv, "n:p:lvh", long_options, &option_index)) != -1) {
        switch (c) {
            case 0:
                printf("option %s", long_options[option_index].name);
                if (optarg) printf(" with arg %s", optarg);
                printf("\n");
                break;
            case 'n':
                options.num_events = atoi(optarg);
                options.num_events_set = true;
                break;
            case 'p':
                options.plot_dir = optarg;
                options.plot_dir_set = true;
                break;
            case 'l':
                options.save_log = true;
                options.save_log_set = true;
                break;
            case 'v':
                printf("Version: %s\n", gitversion);
                exit(0);
            case 'h':
                printf("Usage: %s [OPTION]... RESULTS-FILE INPUT-FILES\n\n", argv[0]);
                printf("Mandatory arguments to long options are mandatory for short options too.\n");
                printf("-h, --help                       display this text and exit\n");
                printf("-l, --log                        save copy of log to results file\n");
                printf("-n, --events=NUM-EVENTS          number of events to be imported from the input file\n");
                printf("-p, --plot-dir=PLOT-DIR          create lifetime/mixing plots\n");
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
