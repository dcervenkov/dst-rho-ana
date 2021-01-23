/**
 *  @file    config.cc
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2019-08-28
 *
 *  @brief Class holding and handling configuration of the fitter, PDFs, etc.
 *
 */

#include "config.h"

// Standard includes
#include <fstream>

// Local includes
#include "constants.h"
#include "log.h"
#include "tools.h"

Config::Config(/* args */) {}

Config::~Config() {}

/**
 * Read a JSON config from a specified filename.
 */
void Config::ReadInJSONFile(const char* filename) {
    std::ifstream filestream(filename);
    if (!filestream.good()) {
        Log::print(Log::error, "Specified config file '%s' doesn't exist!\n", filename);
        exit(5);
    }
    Log::print(Log::debug, "Reading JSON config file '%s'\n", filename);
    filestream >> json;
}

/**
 * Read a JSON config from a specified filename.
 */
void Config::ReadInJSONFile(const std::string filename) { ReadInJSONFile(filename.c_str()); }

/**
 * Update the current config by appending/overwriting values from a JSON
 * object.
 */
void Config::Update(const nlohmann::json new_config) {
    Log::print(Log::debug, "Updating config...\n");
    json.update(new_config);
}

/**
 * Get a pretty-formatted string of the config.
 *
 * @return std::string Pretty-formatted string
 */
std::string Config::GetPrettyString() const { return json.dump(4); }

/**
 * Fill mandatory missing values with sensible defaults
 */
void Config::FillMissingDefaults() {
    if (!json.contains("numCPUs")) {
        json["numCPUs"] = 1;
    }
}

/**
 * Check whether the current config is valid.
 */
bool Config::IsValid() const {
    bool is_valid = true;
    // output check
    if (!json.contains("output")) {
        Log::LogLine(Log::error) << "No output file specified";
        is_valid = false;
    }

    // inputFiles check
    auto input_files = json.find("inputFiles");
    if (input_files == json.end()) {
        Log::LogLine(Log::error) << "Input files not specified";
        is_valid = false;
    } else {
        if (!input_files.value().contains("signal")) {
            Log::LogLine(Log::error) << "Signal window input files not specified";
            is_valid = false;
        }
        if (!input_files.value().contains("sidebands")) {
            Log::LogLine(Log::error) << "Sidebands input files not specified";
            is_valid = false;
        }
    }

    return is_valid;
}
