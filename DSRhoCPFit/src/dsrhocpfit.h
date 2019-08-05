/**
 *  @file    dsrhocpfit.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2016-06-24
 *
 *  @brief Main header
 *
 */

#ifndef DSRHOCPFIT_H_
#define DSRHOCPFIT_H_

// Standard includes
#include <string>
#include <vector>

struct fitter_options {
    int num_CPUs;
    bool num_CPUs_set;
    bool do_mixing_fit;
    bool do_mixing_fit_set;
    bool do_time_independent_fit;
    bool do_time_independent_fit_set;
    char* fix;
    bool fix_set;
    char* plot_dir;
    bool plot_dir_set;
    int num_events;
    bool num_events_set;
    bool perfect_tagging;
    bool perfect_tagging_set;
    std::vector<std::string> efficiency_files;
    bool efficiency_files_set;
    int efficiency_model;
    bool efficiency_model_set;
    bool generator_level;
    bool generator_level_set;
    char* fit;
    bool fit_set;
    char* config_file;
    bool config_file_set;
    bool save_log;
    bool save_log_set;
    char* scf_kde_file;
    bool scf_kde_file_set;
    char* scf_histo_file;
    bool scf_histo_file_set;
};

int ProcessCmdLineOptions(const int argc, char* const argv[], char**& optionless_argv,
                          fitter_options& options);

#endif /* DSRHOCPFIT_H_ */
