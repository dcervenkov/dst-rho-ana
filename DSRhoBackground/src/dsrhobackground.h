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
    bool KDE;
    bool histo;
    bool physics;
    int num_CPUs;
    bool num_CPUs_set;
    char* plot_dir;
    bool plot_dir_set;
};

int ProcessCmdLineOptions(const int argc, char* const argv[], char**& optionless_argv,
                          fitter_options& options);

#endif /* DSRHOBACKGROUND_H_ */
