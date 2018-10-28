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

// ROOT includes
#include "TSystem.h"

// Local includes
#include "colors.h"
#include "constants.h"
#include "fittercpv.h"
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

    char** optionless_argv = NULL;
    // The {} causes the struct's members to be initialized to 0. Without it
    // they would have unspecified values
    fitter_options options = {};
    const int optionless_argc = ProcessCmdLineOptions(argc, argv, optionless_argv, options);

    if (optionless_argc != 19 && optionless_argc != 3) {
        printf("ERROR: Wrong number of arguments.\n");
        printf(
            "Usage: %s [OPTION]... -- INPUT-FILE OUTPUT_DIR "
            "[AP APA A0 ATA XP X0 XT YP Y0 YT XPB X0B XTB YPB Y0B YTB]\n",
            optionless_argv[0]);
        return 2;
    }

    /// This is so I have to change only the next block if I change the
    /// ordering, etc. of arguments
    const char* file_path = optionless_argv[1];
    const char* results_path = optionless_argv[2];
    std::array<double, 16> par_input;

    if (optionless_argc == 19) {
        for (size_t i = 0; i < par_input.size(); i++) {
            par_input[i] = atof(optionless_argv[i + 3]);
        }
    } else {
        par_input = constants::par_input;
    }

    if (options.save_log_set && options.save_log) {
        std::cout.rdbuf(out.rdbuf());
    }

    tools::SetupPlotStyle();
    colors::setColors();

    FitterCPV fitter(par_input);

    std::string output_filename(results_path);
    output_filename += ".root";
    fitter.SetOutputFile(output_filename.c_str());

    if (options.perfect_tagging_set) fitter.SetPerfectTagging(options.perfect_tagging);

    if (options.num_events_set) {
        fitter.ReadInFile(file_path, options.num_events);
    } else {
        fitter.ReadInFile(file_path);
    }

    // Config from file is applied first, so that it can be overriden by CLI
    // switches. However, ReadInFile() has to be called first, as some parts of
    // a config need the data to be present.
    if (options.config_file_set) {
        rapidjson::Document config = FitterCPV::ReadJSONConfig(options.config_file);
        fitter.ApplyJSONConfig(config);
    }

    if (options.num_CPUs_set) fitter.SetNumCPUs(options.num_CPUs);
    if (options.efficiency_model_set) {
        fitter.SetEfficiencyModel(options.efficiency_model);
    } else {
        fitter.SetEfficiencyModel(1);
    }
    if (options.plot_dir_set) fitter.SetPlotDir(options.plot_dir);
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
            // printf("ERROR: Time independent CRSCF fit not implemented!\n");
            // return 2;
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
    fitter.LogText("input_file_name", file_path);
    fitter.LogFileCRC("input_file_crc", file_path);
    fitter.LogText("efficiency_model", std::to_string(fitter.GetEfficiencyModel()).c_str());
    // TODO: This definitely shouldnt't be hardcoded
    fitter.LogFileCRC("meerkat_efficiency_crc", "efficiency");
    fitter.LogFileCRC("histo_efficiency_crc", "efficiency.root");
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
        {"efficiency-model", required_argument, 0, 'e'},
        {"events", required_argument, 0, 'n'},
        {"fit", required_argument, 0, 'f'},
        {"fix", required_argument, 0, 'x'},
        {"log", no_argument, 0, 'l'},
        {"mixing", no_argument, 0, 'm'},
        {"time-independent", no_argument, 0, 'i'},
        {"perfect-tag", no_argument, 0, 't'},
        {"plot-dir", required_argument, 0, 'p'},
        {"help", no_argument, 0, 'h'},
        {NULL, no_argument, NULL, 0}};
    int option_index = 0;
    while ((c = getopt_long(argc, argv, "c:g:e:n:f:x:p:lmith", long_options, &option_index)) != -1) {
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
            case 'h':
                printf("Usage: %s [OPTION]... INPUT-FILE OUTPUT_DIR\n\n", argv[0]);
                printf("Mandatory arguments to long options are mandatory for short options too.\n");
                printf("-c, --cpus=NUM_CPUS              number of CPU cores to use for fitting and plotting\n");
                printf("-e, --efficiency-model=MODEL_NUM number of the efficiency model to be used\n");
                printf("-f, --fit=CR|CRSCF|all           do a specified fit type\n");
                printf("-g, --config=CONFIG_FILE         read in configuration from the specified file");
                printf("-h, --help                       display this text and exit\n");
                printf("-i, --time-independent           make a time-independent fit\n");
                printf("-l, --log                        save copy of log to results file\n");
                printf("-m, --mixing                     make a mixing fit\n");
                printf("-n, --events=NUM_EVENTS          number of events to be imported from the input file\n");
                printf("-p, --plot-dir=PLOT_DIR          create lifetime/mixing plots\n");
                printf("-t, --perfect-tag                use MC info to get perfect tagging\n");
                printf("-x, --fix=ARG1,ARG2,...          fix specified argument(s) to input values in the fit\n");
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
