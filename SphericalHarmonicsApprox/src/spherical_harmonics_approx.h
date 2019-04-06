/**
 *  @file    spherical_harmonics_approx.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2019-04-06
 *
 *  @brief Main header
 *
 */

#pragma once

struct fitter_options {
    char* plot_dir;
    bool plot_dir_set;
    int num_events;
    bool num_events_set;
    bool save_log;
    bool save_log_set;
};

int ProcessCmdLineOptions(const int argc, char* const argv[], char**& optionless_argv,
                          fitter_options& options);
