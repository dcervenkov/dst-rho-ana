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
 * Remove channels that were excluded by 'excludedChannels' key from the JSON
 */
void Config::RemoveExcludedChannels() {
    if (json.contains("excludeChannels")) {
        std::vector<std::string> excluded_channels =
            tools::SplitString(json["excludeChannels"], ',');
        for (auto excluded_channel : excluded_channels) {
            Log::LogLine(Log::info) << "Excluding channel " << excluded_channel;
            json["channels"].erase(excluded_channel);
        }
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

    // components check
    auto components = json.find("components");
    if (components == json.end()) {
        Log::LogLine(Log::error) << "Fit components not specified";
        is_valid = false;
    } else {
        if (components.value() != "CR" && components.value() != "CRSCF" &&
            components.value() != "all") {
            Log::LogLine(Log::error) << "Invalid fit components " << components.value();
            is_valid = false;
        }
    }

    // MC check
    if (!json.contains("MC")) {
        Log::LogLine(Log::error) << "Not specified whether to fit MC or data";
        is_valid = false;
    }

    if (json.contains("excludeChannels")) {
        std::vector<std::string> excluded_channels =
            tools::SplitString(json["excludeChannels"], ',');
        for (auto excluded_channel : excluded_channels) {
            if (!json["channels"].contains(excluded_channel)) {
                Log::LogLine(Log::error) << "Channel '" << excluded_channel
                                         << "' which is to be excluded is not present in config";
                is_valid = false;
            }
        }
    }

    for (auto& chan : json["channels"].items()) {
        const char* channel_name = chan.key().c_str();
        nlohmann::json channel = chan.value();

        // efficiencyModel check
        auto efficiency_model = channel.find("efficiencyModel");
        if (efficiency_model == channel.end()) {
            Log::LogLine(Log::error)
                << "No efficiency model specified for channel " << channel_name;
            is_valid = false;
        } else {
            if (efficiency_model.value() < 0 ||
                efficiency_model.value() > constants::max_efficiency_model_number) {
                Log::LogLine(Log::error)
                    << "Efficiency model " << efficiency_model.value() << " does not exist";
                is_valid = false;
            }
        }

        // efficiencyFile check
        if (efficiency_model.value() > 4 && !channel.contains("efficiencyFile")) {
            Log::LogLine(Log::error) << "Efficiency file not specified";
            is_valid = false;
        }

        // scfHisto vs scfKDE
        if (channel.contains("scfHisto") && channel.contains("scfKDE")) {
            Log::LogLine(Log::error) << "Both SCF histo and KDE specified";
            is_valid = false;
        }

        // inputFiles check
        auto input_files = channel.find("inputFiles");
        if (input_files == channel.end()) {
            Log::LogLine(Log::error) << "Input files not specified for channel " << channel_name;
            is_valid = false;
        } else {
            if (json["MC"] == true) {
                if (!input_files.value().contains("signalMC")) {
                    Log::LogLine(Log::error) << "Signal MC input files not specified";
                    is_valid = false;
                }
                if (!input_files.value().contains("genericMC") && components.value() == "all") {
                    Log::LogLine(Log::error) << "Generic MC input files not specified";
                    is_valid = false;
                }
            } else {
                if (!input_files.value().contains("data")) {
                    Log::LogLine(Log::error) << "Data input files not specified";
                    is_valid = false;
                }
            }
        }
    }

    return is_valid;
}
