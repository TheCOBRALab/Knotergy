#pragma once

#include <string>
#include <vector>

namespace knotergy {

enum class ParamStatus {
    LoadedUserFile,
    Defaulted,
    Fallback,
    Missing,
};

struct ParamMessage {
    std::string text;
    bool warning = false;
};

struct ParamSourceInfo {
    std::string label;  // "ViennaRNA", "Pseudoknot", "Modified Bases"
    ParamStatus status{};
    std::string requested_path;  // what user passed
    std::string resolved_path;   // what actually got used
    std::string resolved_name;   // human-readable parameter set name
    std::vector<ParamMessage> messages;
};

struct ParameterLoadReport {
    ParamSourceInfo vienna;
    ParamSourceInfo pseudoknot;
    ParamSourceInfo modified_bases;
};

}  // namespace knotergy