/**
 *  @file    dsrhoyield.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2015-07-27
 *
 *  @brief Main header
 *
 */

#pragma once

struct fit_options {
    bool mc = false;
    int cpus = 1;
    int rbin;
    bool rbin_set = false;
};

int ProcessCmdLineOptions(const int argc, char* const argv[], char**& optionless_argv,
                          fit_options& options);
