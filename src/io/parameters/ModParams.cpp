#include "ModParams.hpp"

#include "utils/FileUtils.hpp"

#include <iostream>
namespace knotergy {

const modified_base_param* all_mod_params::get_modified_base_param(
    const std::string& modified_base) const {
    auto it = mod_param_lookup.find(modified_base);
    if (it != mod_param_lookup.end()) {
        return it->second;
    }
    return nullptr;
}

const std::string* all_mod_params::get_unmodified_base(const std::string& modified_base) const {
    auto it = mod_to_unmod_lookup.find(modified_base);
    if (it != mod_to_unmod_lookup.end()) {
        return it->second;
    }
    return nullptr;
}

void all_mod_params::build_lookup() {
    mod_param_lookup.reserve(mod_params_.size());
    mod_to_unmod_lookup.reserve(mod_params_.size());
    for (const modified_base_param& param : mod_params_) {
        mod_param_lookup[param.modified_base()]    = &param;
        mod_to_unmod_lookup[param.modified_base()] = &param.fallback_base();
    }
}

// ---------------- Load Modified Base Parameters ----------------

std::vector<modified_base_param> ModParams::load_modified_energy_parameters(
    const std::string& path) {
    if (path.empty()) {
        return {};
    }

    if (!FileUtils::file_exists(path)) {
        THROW_ERROR("Modified parameters JSON file \"" + path + "\" not found.");
    }

    string_list all_files;
    if (FileUtils::is_file(path)) {
        all_files.push_back(path);
    } else if (FileUtils::is_directory(path)) {
        all_files = FileUtils::get_files_in_dir(path);
    }

    std::vector<modified_base_param> params_list;

    for (const std::string& file : all_files) {
        if (file.size() >= 5 && file.substr(file.size() - 5) == ".json") {
            params_list.push_back(parse_modified_base_json(file));
        }
    }

    return params_list;
}

modified_base_param ModParams::parse_modified_base_json(const std::string& jsonFile) {
    std::ifstream f(jsonFile);
    if (!f.is_open()) {
        THROW_ERROR("Error: Unable to open modified base parameter file: " + jsonFile);
    }
    json        data = json::parse(f);
    const json& mod  = data.at("modified_base");

    warn_if_missing(mod, "name", jsonFile);
    warn_if_missing(mod, "unmodified", jsonFile);
    warn_if_missing(mod, "one_letter_code", jsonFile);
    warn_if_missing(mod, "fallback", jsonFile);
    warn_if_missing(mod, "pairing_partners", jsonFile);

    modified_base_param params(
        mod.value("name", ""), mod.value("unmodified", ""), mod.value("one_letter_code", ""),
        mod.value("fallback", ""), mod.value("pairing_partners", string_list{}),
        mod.value("stacking_energies", param_map{}), mod.value("stacking_enthalpies", param_map{}),
        mod.value("terminal_energies", param_map{}), mod.value("terminal_enthalpies", param_map{}),
        mod.value("mismatch_energies", param_map{}), mod.value("mismatch_enthalpies", param_map{}),
        mod.value("dangle5_energies", param_map{}), mod.value("dangle5_enthalpies", param_map{}),
        mod.value("dangle3_energies", param_map{}), mod.value("dangle3_enthalpies", param_map{}));

    return params;
}

void ModParams::warn_if_missing(const json& j, const std::string& key, const std::string& file) {
    if (!j.contains(key)) {
        std::cerr << "Warning: modified_base missing required key '" << key << "' in file " << file
                  << std::endl;
    }
}

}  // namespace knotergy