/**
 *  @file    dsrhobackground.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2017-12-15
 *
 *  @brief Main header
 *
 */

#ifndef DSRHOBACKGROUND_H_
#define DSRHOBACKGROUND_H_

struct fitter_options {
    int num_CPUs;
    bool num_CPUs_set;
    bool do_mixing_fit;
    bool do_mixing_fit_set;
    char* fix;
    bool fix_set;
    char* plot_dir;
    bool plot_dir_set;
    int num_events;
    bool num_events_set;
    bool perfect_tagging;
    bool perfect_tagging_set;
    int efficiency_model;
    bool efficiency_model_set;
};

int ProcessCmdLineOptions(const int argc, char* const argv[], char**& optionless_argv,
                          fitter_options& options);

#endif /* DSRHOBACKGROUND_H_ */