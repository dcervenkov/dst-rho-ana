/**
 *  @file    dsrhocpfit.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2016-06-24
 *
 *  @brief Main header
 *
 */

#pragma once

// Local includes
#include "nlohmann/json.hpp"

int ProcessCmdLineOptions(const int argc, char* const argv[], char**& optionless_argv,
                          nlohmann::json& config);
