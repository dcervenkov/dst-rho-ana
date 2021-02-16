/**
 *  @file    dsrhosidebands.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2020-05-31
 *
 *  @brief Main file
 *
 */

#include "dsrhosidebands.h"

// Standard includes
#include <fstream>
#include <getopt.h>
#include <stdio.h>
#include <ostream>
#include <vector>

// Boost includes
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>

// ROOT includes
#include "TFile.h"

// Local includes
#include "colors.h"
#include "config.h"
#include "constants.h"
#include "fitterdelta.h"
#include "gitversion.h"
#include "log.h"
#include "teebuf.h"
#include "tools.h"

int main(int argc, char* argv[]) {
    // // This block sets up an alternative ostream that can be used to duplicate
    // // cout to both the screen and a stringstream. We use this to save a copy of
    // // the log to the resulting ROOT file.
    // std::ostringstream sout;
    // teebuf sbuf(sout.rdbuf(), std::cout.rdbuf());
    // std::ostream out(&sbuf);
    // std::streambuf* orig_cout_streambuf = std::cout.rdbuf();

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

    if (config.json.contains("plotDir")) {
        tools::SetPlotDir(config.json["plotDir"].get<std::string>().c_str());
    }

    tools::SetupPlotStyle();
    colors::setColors();

    FitterDelta fitter(config.json);


    if (boost::algorithm::ends_with(gitversion, "-dirty")) {
        Log::print(Log::warning, "Using version from dirty Git worktree\n");
    }


    std::string output_filename = config.GetOutputFilename();
    // output_filename += ".root";
    tools::CreateDirsIfNecessary(output_filename);
    // TFile* output_file = new TFile(output_filename.c_str(), "RECREATE");
    // fitter.SetOutputFile(output_file);

    nlohmann::json json_results;
    if (config.json.contains("randomize")) {
        Log::LogLine(Log::debug) << "Generating randomized results";
        fitter.PrintCovarianceMatrix();
        for (int i = 0; i < config.json["randomize"].get<int>(); i++) {
            json_results = fitter.GetJSONResults("bkg_thetab_corr_", true);

            // Save formatted results to file
            if (json_results.size()) {
                std::string output_filename_rand = output_filename;
                tools::RemoveSubstring(output_filename_rand, ".json");
                output_filename_rand += "_rnd_";
                output_filename_rand += std::to_string(i);
                output_filename_rand += ".json";
                std::cout << json_results.dump(2) << std::endl;
                std::ofstream ofs(output_filename_rand);
                ofs << std::setw(2) << json_results << std::endl;
            }
        }
    } else {
        json_results = fitter.GetJSONResults("bkg_thetab_corr_", false);

        // Save formatted results to file
        if (json_results.size()) {
            std::cout << json_results.dump(2) << std::endl;
            std::ofstream ofs(output_filename);
            ofs << std::setw(2) << json_results << std::endl;
        }
    }

    // if (config.ShouldSaveLog()) {
    //     // We have to restore cout's buffer to the original, otherwise we would
    //     // get a segfault as our object goes out of scope sooner than cout
    //     std::cout.rdbuf(orig_cout_streambuf);
    //     Log::print(Log::info, "Saving log to ROOT file\n");
    //     tools::LogText(output_file, "log", sout.str().c_str());
    // }

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
                                    {"plot-dir", required_argument, 0, 'p'},
                                    {"output", required_argument, 0, 'o'},
                                    {"log", no_argument, 0, 'l'},
                                    {"randomize", required_argument, 0, 'r'},
                                    {"version", no_argument, 0, 'v'},
                                    {"help", no_argument, 0, 'h'},
                                    {nullptr, no_argument, nullptr, 0}};
    int option_index = 0;
    while ((c = getopt_long(argc, argv, "c:g:p:o:r:vlh", long_options, &option_index)) != -1) {
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
            case 'l':
                config["saveLog"] = true;
                break;
            case 'o':
                config["output"] = optarg;
                break;
            case 'p':
                config["plotDir"] = optarg;
                break;
            case 'r':
                config["randomize"] = atoi(optarg);
                break;
            case 'v':
                printf("Version: %s\n", gitversion);
                exit(0);
                break;
            case 'h':
                printf("Usage: %s [OPTION]... RESULTS-FILE INPUT-FILES\n\n", argv[0]);
                printf("Mandatory arguments to long options are mandatory for short options too.\n");
                printf("-c, --cpus=NUM_CPUS       number of CPU cores to use for fitting and plotting\n");
                printf("-g, --config=FILE         specify a file from which to read config\n");
                printf("-h, --help                display this text and exit\n");
                printf("-l, --log                 save copy of log to results file\n");
                printf("-o, --output=BASENAME     basename of the output files\n");
                printf("-p, --plot-dir=PLOT_DIR   create lifetime/mixing plots\n");
                printf("-r, --randomize=NUMBER    generate NUMBER randomized results\n");
                printf("-v, --version             print program version and exit\n");
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
