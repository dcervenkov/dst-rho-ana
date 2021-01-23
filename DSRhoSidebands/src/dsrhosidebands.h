/**
 *  @file    dsrhosidebands.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2020-05-31
 *
 *  @brief Main header
 *
 */

#pragma once

#include "nlohmann/json.hpp"

int ProcessCmdLineOptions(const int argc, char* const argv[], char**& optionless_argv,
                          nlohmann::json& config);
