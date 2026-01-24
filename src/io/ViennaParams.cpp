#include "ViennaParams.hpp"

namespace knotergy {

void ViennaParams::load_energy_parameters(const std::string& paramFile, int dangle,
                                          const std::string& seq) {
    if (p) {
        free(p);
        p = nullptr;
    }  // free previous params if reloading

    // Find user specified parameter file
    vrna_md_set_default(&ViennaParams::md);
    ViennaParams::md.dangles = dangle;
    if (!paramFile.empty()) {
        if (FileUtils::file_exists(paramFile)) {
            int loaded = vrna_params_load(paramFile.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
            if (!loaded) {
                THROW_ERROR("Failed to load parameter file: " + paramFile);
            }
            std::cout << "Successfully loaded parameter file: " << paramFile << std::endl;
            ViennaParams::p = vrna_params(&ViennaParams::md);
            return;
        } else {
            std::cerr << "Warning: Parameter file \"" << paramFile << "\" not found." << std::endl;
        }
    } else {
        std::cerr << "No parameter file provided. ";
    }

    // Default fallback based on sequence
    if (seq.find('T') != std::string::npos) {  //  (detect DNA but currently disabled)
        std::cerr << "Defaulting to DNA parameters (Mathews 2004)." << std::endl;
        vrna_params_load_DNA_Mathews2004();
    } else {
        // Default RNA parameters (Dirks&Pierce 2009)
        const std::string default_path = "./params/common/rna_DirksPierce09.par";
        std::cerr << "Defaulting to RNA parameters (Dirks&Pierce 2009)." << std::endl;
        int loaded = vrna_params_load(default_path.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);

        // in case someone deletes the default param file 😭
        if (!loaded) {
            std::cerr << ("Failed to load parameter file: " + default_path) << std::endl;
            std::cerr << "Defaulting to RNA parameters (Turner 2004)." << std::endl;
            vrna_params_load_RNA_Turner2004();
        }
    }
    ViennaParams::p = vrna_params(&ViennaParams::md);
    return;
}

std::vector<modified_base_params> ViennaParams::load_modified_energy_parameters(
    const std::string& path) {
    if (path.empty()) {
        return {};
    }

    if (!FileUtils::file_exists(path)) {
        THROW_ERROR("Modified parameters JSON file \"" + path + "\" not found.");
    }

    std::vector<std::string> all_files;
    if (FileUtils::is_file(path)) {
        all_files.push_back(path);
    } else if (FileUtils::is_directory(path)) {
        all_files = FileUtils::get_files_in_dir(path);
    }

    std::vector<modified_base_params> params_list;

    for (const std::string& file : all_files) {
        if (file.size() >= 5 && file.substr(file.size() - 5) == ".json") {
            params_list.push_back(parse_modified_base_json(file));
        }
    }

    return params_list;
}

modified_base_params ViennaParams::parse_modified_base_json(const std::string& jsonFile) {
    std::ifstream f(jsonFile);
    if (!f.is_open()) {
        THROW_ERROR("Error: Unable to open modified base parameter file: " + jsonFile);
    }
    json data = json::parse(f);
    json mod = data["modified_base"].get<json>();

    warn_if_missing(mod, "name", jsonFile);
    warn_if_missing(mod, "unmodified", jsonFile);
    warn_if_missing(mod, "one_letter_code", jsonFile);
    warn_if_missing(mod, "fallback", jsonFile);
    warn_if_missing(mod, "pairing_partners", jsonFile);

    modified_base_params params(mod.value("name", ""), mod.value("unmodified", ""),
                                mod.value("one_letter_code", ""), mod.value("fallback", ""),
                                mod.value("pairing_partners", std::vector<std::string>{}),
                                mod.value("stacking_energies", std::map<std::string, float>{}),
                                mod.value("stacking_enthalpies", std::map<std::string, float>{}),
                                mod.value("terminal_energies", std::map<std::string, float>{}),
                                mod.value("terminal_enthalpies", std::map<std::string, float>{}),
                                mod.value("mismatch_energies", std::map<std::string, float>{}),
                                mod.value("mismatch_enthalpies", std::map<std::string, float>{}),
                                mod.value("dangle5_energies", std::map<std::string, float>{}),
                                mod.value("dangle5_enthalpies", std::map<std::string, float>{}),
                                mod.value("dangle3_energies", std::map<std::string, float>{}),
                                mod.value("dangle3_enthalpies", std::map<std::string, float>{}));

    return params;
}

void ViennaParams::warn_if_missing(const json& j, const std::string& key, const std::string& file) {
    if (!j.contains(key)) {
        std::cerr << "Warning: modified_base missing required key '" << key << "' in file " << file
                  << std::endl;
    }
}

}  // namespace knotergy
