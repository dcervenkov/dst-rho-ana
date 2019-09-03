/**
 *  @file    config.h
 *  @author  Daniel Cervenkov, cervenkov(at)ipnp.mff.cuni.cz
 *  @date    2019-08-28
 *
 *  @brief Class holding and handling configuration of the fitter, PDFs, etc.
 *
 */

#pragma once

// Standard includes
#include <string>

// Local includes
#include "nlohmann/json.hpp"

class Config
{
private:
public:
    Config(/* args */);
    ~Config();
    void ReadInJSONFile(const char* filename);
    void ReadInJSONFile(const std::string filename);
    void Update(const nlohmann::json new_config);
    std::string GetPrettyString() const;
    bool IsValid() const;
    bool ShouldSaveLog() const { return json.contains("saveLog"); };
    std::string GetOutputFilename() const { return json["output"].get<std::string>(); };

    nlohmann::json json;
};
