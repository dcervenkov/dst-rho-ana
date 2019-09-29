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
#include <array>
#include <cstdio>
#include <cstring>
#include <ostream>
#include <sstream>

// Boost includes
#include <boost/algorithm/string/predicate.hpp>

// ROOT includes
#include "TFile.h"
#include "TSystem.h"

// Local includes
#include "colors.h"
#include "config.h"
#include "constants.h"
#include "fittercpv.h"
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
    nlohmann::json cli_config;
    ProcessCmdLineOptions(argc, argv, optionless_argv, cli_config);

    Config config;
    if (cli_config.contains("configFile")) {
        config.ReadInJSONFile(cli_config["configFile"]);
    }

    config.Update(cli_config);
    config.FillMissingDefaults();
    if (!config.IsValid()) {
        Log::print(Log::error, "Config is not valid!\n");
        exit(1);
    }
    config.RemoveExcludedChannels();

    if (config.ShouldSaveLog()) {
        Log::print(Log::debug, "Manipulating cout buffer to save copy of stdout\n");
        std::cout.rdbuf(out.rdbuf());
    }

    if (boost::algorithm::ends_with(gitversion, "-dirty")) {
        Log::print(Log::warning, "Using version from dirty Git worktree\n");
    }

    tools::SetupPlotStyle();
    colors::setColors();

    FitterCPV fitter(config.json);

    // if (options.generator_level_set) fitter.SetGeneratorLevel(options.generator_level);
    // if (options.plot_dir_set) fitter.SetPlotDir(options.plot_dir);

    // fitter.TestEfficiency();
    // fitter.PlotEfficiency();

    std::string output_filename = config.GetOutputFilename();
    output_filename += ".root";
    TFile* output_file = new TFile(output_filename.c_str(), "RECREATE");
    fitter.SetOutputFile(output_file);

    fitter.Fit(config.json);

    tools::LogCLIArguments(output_file, argc, argv);
    tools::LogEnvironmentMetadata(output_file);
    tools::LogText(output_file, "config", config.GetPrettyString());

    if (fitter.ResultExists()) {
        fitter.LogResults();
        tools::LogText(output_file, "pull_table", fitter.CreatePullTableString().c_str());
        tools::LogText(output_file, "pull_table_asym", fitter.CreatePullTableString(true).c_str());
        tools::LogText(output_file, "latex_pull_table",
                       fitter.CreateLatexPullTableString().c_str());
        tools::LogText(output_file, "latex_pull_table_asym",
                       fitter.CreateLatexPullTableString(true).c_str());
        tools::SaveTextToFile(config.GetOutputFilename(), fitter.CreateResultsString());
    }

    if (config.json.contains("plotDir")) {
        tools::SetPlotDir(config.json["plotDir"].get<std::string>().c_str());
        fitter.CreatePlots(config.json);
    }

    if (config.ShouldSaveLog()) {
        // We have to restore cout's buffer to the original, otherwise we would
        // get a segfault as our object goes out of scope sooner than cout
        std::cout.rdbuf(orig_cout_streambuf);
        Log::print(Log::info, "Saving log to ROOT file\n");
        tools::LogText(output_file, "log", sout.str().c_str());
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
                          nlohmann::json& config) {
    int c;
    struct option long_options[] = {{"cpus", required_argument, 0, 'c'},
                                    {"config", required_argument, 0, 'g'},
                                    {"components", required_argument, 0, 'e'},
                                    {"events", required_argument, 0, 'n'},
                                    {"exclude-channels", required_argument, 0, 'x'},
                                    {"MC", required_argument, 0, 'm'},
                                    {"fix", required_argument, 0, 'f'},
                                    {"plot-dir", required_argument, 0, 'p'},
                                    {"output", required_argument, 0, 'o'},
                                    {"log", no_argument, 0, 'l'},
                                    {"time-independent", no_argument, 0, 'i'},
                                    {"perfect-tag", no_argument, 0, 't'},
                                    {"generator-level", no_argument, 0, 'r'},
                                    {"version", no_argument, 0, 'v'},
                                    {"help", no_argument, 0, 'h'},
                                    {nullptr, no_argument, nullptr, 0}};
    int option_index = 0;
    while ((c = getopt_long(argc, argv, "c:g:e:n:x:m:f:p:o:litrvh", long_options,
                            &option_index)) != -1) {
        switch (c) {
            case 0:
                printf("option %s", long_options[option_index].name);
                if (optarg) printf(" with arg %s", optarg);
                printf("\n");
                break;
            case 'c':
                config["numCPUs"] = atoi(optarg);
                break;
            case 'g':
                config["configFile"] = optarg;
                break;
            case 'e':
                config["components"] = optarg;
                break;
            case 'n':
                config["events"] = atoi(optarg);
                break;
            case 'x':
                config["excludeChannels"] = optarg;
                break;
            case 'm':
                config["MC"] = bool(atoi(optarg));
                break;
            case 'f':
                config["fixedParameters"] = optarg;
                break;
            case 'p':
                config["plotDir"] = optarg;
                break;
            case 'o':
                config["output"] = optarg;
                break;
            case 'l':
                config["saveLog"] = true;
                break;
            case 'i':
                config["timeIndependent"] = true;
                break;
            case 't':
                config["perfectTagging"] = true;
                break;
            case 'r':
                config["generatorLevel"] = true;
                break;
            case 'v':
                printf("Version: %s\n", gitversion);
                exit(0);
            case 'h':
                printf("Usage: %s [OPTION]... RESULTS-FILE INPUT-FILES\n\n", argv[0]);
                printf("Mandatory arguments to long options are mandatory for short options too.\n");
                printf("-c, --cpus=NUM_CPUS              number of CPU cores to use for fitting and plotting\n");
                printf("-e, --components=CR|CRSCF|all    do a specified fit type\n");
                printf("-f, --fix=ARG1,ARG2,...          fix specified argument(s) to input values in the fit;\n");
                printf("                                 additional short-hand ARGs are: all, xy, trans and nota0\n");
                printf("-g, --config=CONFIG-FILE         read in configuration from the specified file\n");
                printf("-h, --help                       display this text and exit\n");
                printf("-i, --time-independent           make a time-independent fit\n");
                printf("-l, --log                        save copy of log to results file\n");
                printf("-m, --MC=0|1                     whether to fit MC or data\n");
                printf("-o, --output                     basename of the output files\n");
                printf("-p, --plot-dir=PLOT-DIR          create lifetime/mixing plots\n");
                printf("-r, --generator-level            do a generator level fit\n");
                printf("-t, --perfect-tag                use MC info to get perfect tagging\n");
                printf("-v, --version                	 display version and exit\n");
                printf("-x, --exclude-channels         	 exclude channels from fit\n");
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
