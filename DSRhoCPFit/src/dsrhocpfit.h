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

// Local includes
#include "nlohmann/json.hpp"

int ProcessCmdLineOptions(const int argc, char* const argv[], char**& optionless_argv,
                          nlohmann::json& config);

#endif /* DSRHOCPFIT_H_ */
