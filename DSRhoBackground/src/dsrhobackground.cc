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
    std::string JSON_formatted_results;

    if (options.KDE) {
        fitter.thetat_.setBins(50);
        fitter.thetab_.setBins(50);
        fitter.phit_.setBins(50);
        AdaptiveKernelDensity kde = fitter.CreateKDEPDF(fitter.dataset_, results_file);
        fitter.PlotKDE(kde);
    } else if (options.histo) {
        fitter.CreateHistoPDF(fitter.dataset_, results_file);
    } else if (options.physics) {
        fitter.Fit(fitter.bkg_physics_dt_model_, fitter.dataset_);
        JSON_formatted_results += tools::FormatResultsJSON(
            "Physics-based dt Parameters", fitter.bkg_physics_dt_model_, observables);
        if (options.plot_dir_set) {
            fitter.PlotWithPull(fitter.dt_, *fitter.dataset_, *fitter.bkg_physics_dt_model_);
        }
    } else {
        // We must plot after each fit otherwise the results printed on the plot
        // would be of the last fit
        fitter.Fit(&fitter.bkg_phit_model_, fitter.dataset_);
        if (options.plot_dir_set) {
            fitter.PlotWithPull(fitter.phit_, *fitter.dataset_, fitter.bkg_phit_model_);
        }

        fitter.Fit(&fitter.bkg_thetat_model_, fitter.dataset_);
        if (options.plot_dir_set) {
            fitter.PlotWithPull(fitter.thetat_, *fitter.dataset_, fitter.bkg_thetat_model_);
        }

        std::vector<const RooAbsPdf*> angular_pdfs = {
            &fitter.bkg_thetab_model_, &fitter.bkg_thetat_model_, &fitter.bkg_phit_model_};
        fitter.Fit(&fitter.bkg_thetab_model_, fitter.dataset_);
        JSON_formatted_results += tools::FormatResultsJSON("Angular PDFs' Parameters", angular_pdfs, observables);
        if (options.plot_dir_set) {
            fitter.PlotWithPull(fitter.thetab_, *fitter.dataset_, fitter.bkg_thetab_model_);
        }

        fitter.Fit(&fitter.bkg_dt_model_, fitter.dataset_);
        JSON_formatted_results += tools::FormatResultsJSON("dt Parameters", &fitter.bkg_dt_model_, observables);
        if (options.plot_dir_set) {
            fitter.PlotWithPull(fitter.dt_, *fitter.dataset_, fitter.bkg_dt_model_);
        }

        fitter.Fit(&fitter.bkg_dt_model_, fitter.dataset_a_);
        JSON_formatted_results += tools::FormatResultsJSON("dt a Parameters", &fitter.bkg_dt_model_, observables);
        if (options.plot_dir_set) {
            fitter.PlotWithPull(fitter.dt_, *fitter.dataset_a_, fitter.bkg_dt_model_);
        }

        fitter.Fit(&fitter.bkg_dt_model_, fitter.dataset_ab_);
        JSON_formatted_results += tools::FormatResultsJSON("dt ab Parameters", &fitter.bkg_dt_model_, observables);
        if (options.plot_dir_set) {
            fitter.PlotWithPull(fitter.dt_, *fitter.dataset_ab_, fitter.bkg_dt_model_);
        }

        fitter.Fit(&fitter.bkg_dt_model_, fitter.dataset_b_);
        JSON_formatted_results += tools::FormatResultsJSON("dt b Parameters", &fitter.bkg_dt_model_, observables);
        if (options.plot_dir_set) {
            fitter.PlotWithPull(fitter.dt_, *fitter.dataset_b_, fitter.bkg_dt_model_);
        }

        fitter.Fit(&fitter.bkg_dt_model_, fitter.dataset_bb_);
        JSON_formatted_results += tools::FormatResultsJSON("dt bb Parameters", &fitter.bkg_dt_model_, observables);
        if (options.plot_dir_set) {
            fitter.PlotWithPull(fitter.dt_, *fitter.dataset_bb_, fitter.bkg_dt_model_);
        }

        RooDataSet* dataset_cf = static_cast<RooDataSet*>(
                fitter.dataset_->emptyClone("dataset_cf", "dataset_cf"));
        dataset_cf->append(*fitter.dataset_a_);
        dataset_cf->append(*fitter.dataset_ab_);
        fitter.Fit(&fitter.bkg_dt_model_, dataset_cf);
        JSON_formatted_results += tools::FormatResultsJSON("dt cf Parameters", &fitter.bkg_dt_model_, observables);
        if (options.plot_dir_set) {
            fitter.PlotWithPull(fitter.dt_, *dataset_cf, fitter.bkg_dt_model_);
        }

        RooDataSet* dataset_dcs = static_cast<RooDataSet*>(
                fitter.dataset_->emptyClone("dataset_dcs", "dataset_dcs"));
        dataset_dcs->append(*fitter.dataset_b_);
        dataset_dcs->append(*fitter.dataset_bb_);
        fitter.Fit(&fitter.bkg_dt_model_, dataset_dcs);
        JSON_formatted_results += tools::FormatResultsJSON("dt dcs Parameters", &fitter.bkg_dt_model_, observables);
        if (options.plot_dir_set) {
            fitter.PlotWithPull(fitter.dt_, *dataset_dcs, fitter.bkg_dt_model_);
        }

        if (options.plot_dir_set) {
            fitter.thetat_.setBins(50);
            fitter.thetab_.setBins(50);
            fitter.phit_.setBins(50);

            RooProdPdf* bkg_pdf =
                new RooProdPdf("bkg_pdf", "bkg_pdf",
                                RooArgList(fitter.bkg_thetat_model_, fitter.bkg_thetab_model_,
                                            fitter.bkg_phit_model_));

            RooDataHist* pdf_hist = bkg_pdf->generateBinned(RooArgSet(fitter.thetat_,
            fitter.thetab_, fitter.phit_), fitter.dataset_->numEntries(), true);

            RooDataHist data_hist("data_hist", "data_hist",
                            RooArgSet(fitter.thetat_, fitter.thetab_, fitter.phit_),
                            *fitter.dataset_);

            tools::PlotVars2D(fitter.thetat_, fitter.thetab_, data_hist, *pdf_hist, "",
            constants::format); tools::PlotVars2D(fitter.thetat_, fitter.phit_, data_hist,
            *pdf_hist, "", constants::format); tools::PlotVars2D(fitter.thetab_,
            fitter.phit_, data_hist, *pdf_hist, "", constants::format);

            tools::PlotPull2D(fitter.thetat_, fitter.thetab_, data_hist, *pdf_hist, "",
            constants::format, true); tools::PlotPull2D(fitter.thetat_, fitter.phit_,
            data_hist, *pdf_hist, "", constants::format, true);
            tools::PlotPull2D(fitter.thetab_, fitter.phit_, data_hist, *pdf_hist, "",
            constants::format, true);
        }

    }
    if (!JSON_formatted_results.empty()) {
        std::cout << JSON_formatted_results;
        tools::SaveTextToFile(results_file, JSON_formatted_results);
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
                                    {"physics", no_argument, 0, 'y'},
                                    {"version", no_argument, 0, 'v'},
                                    {"help", no_argument, 0, 'h'},
                                    {nullptr, no_argument, nullptr, 0}};
    int option_index = 0;
    while ((c = getopt_long(argc, argv, "c:p:ksyvh", long_options, &option_index)) != -1) {
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
            case 'y':
                options.physics = true;
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
                printf("-y, --physics             fit physics-based dt PDF\n");
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
