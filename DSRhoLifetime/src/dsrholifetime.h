/**
 *  @file    dsrholifetime.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2016-04-05
 *
 *  @brief Main header
 *
 */

#ifndef DSRHOLIFETIME_H_
#define DSRHOLIFETIME_H_

struct fitter_options {
    int num_CPUs;
    bool num_CPUs_set;
    bool do_lifetime_fit;
    bool do_lifetime_fit_set;
    bool do_mixing_fit;
    bool do_mixing_fit_set;
    const char* plot_dir;
    bool plot_dir_set;
    const char* config;
    bool config_set;
    const char* channel;
    bool channel_set;
    const char* components;
    bool components_set;
    int num_events;
    bool num_events_set;
    bool perfect_tagging;
    bool perfect_tagging_set;
};

int ProcessCmdLineOptions(const int argc, char* const argv[], char**& optionless_argv,
                          fitter_options& options);

#endif /* DSRHOLIFETIME_H_ */
