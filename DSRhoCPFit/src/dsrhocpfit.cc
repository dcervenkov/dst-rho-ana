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

    /// This is so I have to change only the next block if I change the
    /// ordering, etc. of arguments
    const char* results_path = optionless_argv[1];

    // TODO: Move par_input to option parameter string
    std::array<double, 16> par_input;

    std::vector<const char*> file_names;
    for (int i = 2; i < optionless_argc; i++) {
        file_names.push_back(optionless_argv[i]);
    }

    par_input = constants::par_input;

    if (options.save_log_set && options.save_log) {
        std::cout.rdbuf(out.rdbuf());
    }

    if (options.scf_kde_file_set && options.scf_histo_file_set) {
        Log::print(Log::error, "Both '--scf-kde' and '--scf-histo' set. Use only one!\n");
        return 3;
    }

    tools::SetupPlotStyle();
    colors::setColors();

    FitterCPV fitter;

    std::string output_filename(results_path);
    output_filename += ".root";
    fitter.SetOutputFile(output_filename.c_str());

    if (options.generator_level_set) fitter.SetGeneratorLevel(options.generator_level);

    fitter.InitVars(par_input);

    if (options.config_file_set) {
        rapidjson::Document config = FitterCPV::ReadJSONConfig(options.config_file);
        fitter.ApplyJSONConfig(config);
    }

    if (options.perfect_tagging_set) fitter.SetPerfectTagging(options.perfect_tagging);

    if (options.num_events_set) {
        fitter.ReadInFile(file_names, options.num_events);
    } else {
        fitter.ReadInFile(file_names);
    }

    if (options.num_CPUs_set) fitter.SetNumCPUs(options.num_CPUs);
    if (options.efficiency_model_set) {
        fitter.SetEfficiencyModel(options.efficiency_model);
    } else {
        fitter.SetEfficiencyModel(1);
    }
    if (options.efficiency_files_set) {
        fitter.SetEfficiencyFiles(options.efficiency_files);
    }
    if (options.plot_dir_set) fitter.SetPlotDir(options.plot_dir);

    if (options.scf_kde_file_set) fitter.SetSCFKDE(options.scf_kde_file);
    if (options.scf_histo_file_set) fitter.SetSCFHisto(options.scf_histo_file);

    if (options.do_mixing_fit_set) fitter.SetDoMixingFit(options.do_mixing_fit);
    if (options.do_time_independent_fit_set) {
        fitter.SetDoTimeIndependentFit(options.do_time_independent_fit);
    } else {
        fitter.SetDoTimeIndependentFit(false);
    }
    if (options.fix_set) {
        if (fitter.FixParameters(options.fix)) {
            return 1;
        }
    }

    if (!options.fit_set) {
        options.fit = (char*)"all";
    }

    // fitter.TestEfficiency();
    // fitter.PlotEfficiency();

    if (std::strcmp(options.fit, "CR") == 0) {
        if (fitter.GetDoTimeIndependentFit()) {
            fitter.FitAngularCR();
        } else {
            fitter.FitCR();
        }
    } else if (std::strcmp(options.fit, "CRSCF") == 0) {
        if (fitter.GetDoTimeIndependentFit()) {
            fitter.FitAngularCRSCF();
        } else {
            fitter.FitCRSCF();
        }
    } else if (std::strcmp(options.fit, "all") == 0) {
        if (fitter.GetDoTimeIndependentFit()) {
            fitter.FitAngularAll();
        } else {
            fitter.FitAll();
        }
    }
    // fitter.GenerateToys(10000, 10);

    fitter.LogCLIArguments(argc, argv);
    fitter.LogEnvironmentMetadata();
    for (auto file_name : file_names) {
        fitter.LogText("input_file_name", file_name);
        fitter.LogFileCRC("input_file_crc", file_name);
    }
    fitter.LogText("efficiency_model", std::to_string(fitter.GetEfficiencyModel()).c_str());
    for (auto efficiency_file : fitter.GetEfficiencyFiles()) {
        fitter.LogText("efficiency_file_name", efficiency_file);
        fitter.LogFileCRC("efficiency_file_crc", efficiency_file);
    }
    if (fitter.ResultExists()) {
        fitter.LogResults();
        fitter.LogText("pull_table", fitter.CreatePullTableString().c_str());
        fitter.LogText("pull_table_asym", fitter.CreatePullTableString(true).c_str());
        fitter.LogText("latex_pull_table", fitter.CreateLatexPullTableString().c_str());
        fitter.LogText("latex_pull_table_asym", fitter.CreateLatexPullTableString(true).c_str());

        fitter.SaveTXTResults(results_path);
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
        {"cpus", required_argument, 0, 'c'},
        {"config", required_argument, 0, 'g'},
        {"efficiency-file", required_argument, 0, 'y'},
        {"efficiency-model", required_argument, 0, 'e'},
        {"events", required_argument, 0, 'n'},
        {"fit", required_argument, 0, 'f'},
        {"fix", required_argument, 0, 'x'},
        {"generator-level", no_argument, 0, 'r'},
        {"log", no_argument, 0, 'l'},
        {"mixing", no_argument, 0, 'm'},
        {"time-independent", no_argument, 0, 'i'},
        {"perfect-tag", no_argument, 0, 't'},
        {"plot-dir", required_argument, 0, 'p'},
        {"scf-kde", required_argument, 0, 'k'},
        {"scf-histo", required_argument, 0, 's'},
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {nullptr, no_argument, nullptr, 0}};
    int option_index = 0;
    while ((c = getopt_long(argc, argv, "c:g:y:e:n:f:x:p:k:s:lmitvh", long_options, &option_index)) != -1) {
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
            case 'g':
                options.config_file = optarg;
                options.config_file_set = true;
                break;
            case 'e':
                options.efficiency_model = atoi(optarg);
                options.efficiency_model_set = true;
                break;
            case 'y':
                options.efficiency_files.push_back(optarg);
                options.efficiency_files_set = true;
                break;
            case 'n':
                options.num_events = atoi(optarg);
                options.num_events_set = true;
                break;
            case 'f':
                options.fit = optarg;
                options.fit_set = true;
                break;
            case 'x':
                options.fix = optarg;
                options.fix_set = true;
                break;
            case 'p':
                options.plot_dir = optarg;
                options.plot_dir_set = true;
                break;
            case 'l':
                options.save_log = true;
                options.save_log_set = true;
                break;
            case 'm':
                options.do_mixing_fit = true;
                options.do_mixing_fit_set = true;
                break;
            case 'i':
                options.do_time_independent_fit = true;
                options.do_time_independent_fit_set = true;
                break;
            case 't':
                options.perfect_tagging = true;
                options.perfect_tagging_set = true;
                break;
            case 'r':
                options.generator_level = true;
                options.generator_level_set = true;
                break;
            case 'k':
                options.scf_kde_file = optarg;
                options.scf_kde_file_set = true;
                break;
            case 'v':
                printf("Version: %s\n", gitversion);
                exit(0);
            case 's':
                options.scf_histo_file = optarg;
                options.scf_histo_file_set = true;
                break;
            case 'h':
                printf("Usage: %s [OPTION]... RESULTS-FILE INPUT-FILES\n\n", argv[0]);
                printf("Mandatory arguments to long options are mandatory for short options too.\n");
                printf("-c, --cpus=NUM_CPUS              number of CPU cores to use for fitting and plotting\n");
                printf("-e, --efficiency-model=MODEL-NUM number of the efficiency model to be used\n");
                printf("-f, --fit=CR|CRSCF|all           do a specified fit type\n");
                printf("-g, --config=CONFIG-FILE         read in configuration from the specified file\n");
                printf("-h, --help                       display this text and exit\n");
                printf("-i, --time-independent           make a time-independent fit\n");
                printf("-k, --scf-kde                    use SCF KDE from file\n");
                printf("-l, --log                        save copy of log to results file\n");
                printf("-m, --mixing                     make a mixing fit\n");
                printf("-n, --events=NUM-EVENTS          number of events to be imported from the input file\n");
                printf("-p, --plot-dir=PLOT-DIR          create lifetime/mixing plots\n");
                printf("-r, --generator-level            do a generator level fit\n");
                printf("-s, --scf-histo                  use histo SCF from file\n");
                printf("-t, --perfect-tag                use MC info to get perfect tagging\n");
                printf("-v, --version                	 display version and exit\n");
                printf("-x, --fix=ARG1,ARG2,...          fix specified argument(s) to input values in the fit;\n");
                printf("                                 additional short-hand ARGs are: all, xy, trans and nota0\n");
                printf("-y, --efficiency-file=ROOT-FILE  file from which to read in efficiency histogram\n");
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
