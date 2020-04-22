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

// Boost includes
#include <boost/filesystem.hpp>

// ROOT includes
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"

// Meerkat includes
#include "AdaptiveKernelDensity.hh"

// Local includes
#include "colors.h"
#include "constants.h"
#include "fitterbkg.h"
#include "gitversion.h"
#include "log.h"
#include "tools.h"

int main(int argc, char* argv[]) {
    Log::setLogLevel(Log::debug);
    char** optionless_argv = nullptr;
    fitter_options options = {};
    const int optionless_argc = ProcessCmdLineOptions(argc, argv, optionless_argv, options);

    if (optionless_argc < 3) {
        printf("ERROR: Wrong number of arguments.\n");
        printf("Usage: %s [OPTION]... RESULTS-FILE INPUT-FILE(S)...\n", optionless_argv[0]);
        return 2;
    }

    std::string results_file = optionless_argv[1];
    std::vector<const char*> file_names;
    if (boost::filesystem::is_regular_file(results_file)) {
        Log::LogLine(Log::warning) << "Rewriting results file '" << results_file << "'";
    }
    for (int i = 2; i < optionless_argc; i++) {
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

    const RooArgSet observables(fitter.thetab_, fitter.thetat_, fitter.phit_, fitter.dt_);
    nlohmann::json json_results;

    if (options.KDE) {
        fitter.thetat_.setBins(50);
        fitter.thetab_.setBins(50);
        fitter.phit_.setBins(50);
        AdaptiveKernelDensity kde = fitter.CreateKDEPDF(fitter.dataset_, results_file);
        fitter.PlotKDE(kde);
    }
    if (options.histo) {
        fitter.CreateHistoPDF(fitter.dataset_, results_file);
    }
    if (options.physics) {
        if (options.notail) {
            fitter.SetNoTailPDF();
        }
        if (options.nodelta) {
            fitter.SetNoDeltaPDF();
        }
        nlohmann::json local_results =
            fitter.FitDt(fitter.physics_dt_model_, "phys_", options.plot_dir_set);
        json_results = tools::MergeJSON(json_results, local_results);
    }
    if (options.empirical) {
        nlohmann::json local_results =
            fitter.FitDt(&fitter.dt_model_, "emp_", options.plot_dir_set);
        json_results = tools::MergeJSON(json_results, local_results);
    }
    if (options.angular) {
        nlohmann::json local_results = fitter.FitAngular(options.plot_dir_set);
        json_results = tools::MergeJSON(json_results, local_results);
    }

    // Save formatted results to file
    if (json_results.size()) {
        std::cout << json_results.dump(2) << std::endl;
        std::ofstream ofs(results_file);
        ofs << std::setw(2) << json_results << std::endl;
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
                                    {"angular", no_argument, 0, 'a'},
                                    {"empirical", no_argument, 0, 'e'},
                                    {"plot-dir", required_argument, 0, 'p'},
                                    {"physics", no_argument, 0, 'y'},
                                    {"nodelta", no_argument, 0, 'd'},
                                    {"notail", no_argument, 0, 't'},
                                    {"version", no_argument, 0, 'v'},
                                    {"help", no_argument, 0, 'h'},
                                    {nullptr, no_argument, nullptr, 0}};
    int option_index = 0;
    while ((c = getopt_long(argc, argv, "c:p:kaesydtvh", long_options, &option_index)) != -1) {
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
            case 'a':
                options.angular = true;
                break;
            case 'e':
                options.empirical = true;
                break;
            case 'y':
                options.physics = true;
                break;
            case 'd':
                options.nodelta = true;
                break;
            case 't':
                options.notail = true;
                break;
            case 'v':
                printf("Version: %s\n", gitversion);
                exit(0);
                break;
            case 'h':
                printf("Usage: %s [OPTION]... INPUT-FILES\n\n", argv[0]);
                printf("Mandatory arguments to long options are mandatory for short options too.\n");
                printf("-c, --cpus=NUM_CPUS       number of CPU cores to use for fitting and plotting\n");
                printf("-k, --kde                 use kernel density estimation\n");
                printf("-s, --histo               create a histogram PDF\n");
                printf("-a, --angular             fit angular PDFs\n");
                printf("-e, --empirical           fit empirical dt PDFs\n");
                printf("-y, --physics             fit physics-based dt PDF\n");
                printf("-d, --nodelta             disable certain parts of physics-based PDF\n");
                printf("-t, --notail              disable certain parts of physics-based PDF\n");
                printf("-h, --help                display this text and exit\n");
                printf("-p, --plot-dir=PLOT_DIR   create lifetime/mixing plots\n");
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
