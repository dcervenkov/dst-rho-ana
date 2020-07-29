/**
 *  @file    dsrhoefficiency.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2015-11-24
 *
 *  @brief Main header
 *
 */

#ifndef DSRHOEFFICIENCY_H_
#define DSRHOEFFICIENCY_H_

// Standard includes
#include <array>

struct fitter_options {
    const char* ada_kde_pars;
    const char* bin_kde_pars;
    int efficiency_model;
    int random_models;
    const char* efficiency_file;
    double mirror_margin;
};

std::array<double, 6> DecodeStringPars(const char* string);

int ProcessCmdLineOptions(const int argc, char* const argv[], char**& optionless_argv,
                          fitter_options& options);

#endif /* DSRHOYIELD_H_ */
