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

struct fitter_options {
    int num_CPUs;
    bool num_CPUs_set;
    bool do_mixing_fit;
    bool do_mixing_fit_set;
    char* fix;
    bool fix_set;
    bool make_plots;
    bool make_plots_set;
    int num_events;
    bool num_events_set;
    bool perfect_tagging;
    bool perfect_tagging_set;
};

int ProcessCmdLineOptions(const int argc, char* const argv[], char**& optionless_argv,
                          fitter_options& options);

#endif /* DSRHOCPFIT_H_ */
