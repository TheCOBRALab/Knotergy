#pragma once

#include <iostream>
#include "shared.hpp"

#include <fstream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

extern "C" {
// used for load_energy_parameters
#include <ViennaRNA/params/io.h>
}

namespace knotergy{
    struct modified_base_params {

        modified_base_params(const std::string& unmod, const std::string& mod,
                             const std::string& fallback,
                             const std::vector<std::string>& partners,
                             const std::map<std::string, float>& stacking,
                             const std::map<std::string, float>& enthalpies,
                             const std::map<std::string, float>& terminal_en,
                             const std::map<std::string, float>& terminal_h)
            : unmodified_base(unmod),
              modified_base(mod),
              fallback_base(fallback),
              pairing_partners(partners),
              stacking_energies(stacking),
              stacking_enthalpies(enthalpies),
              terminal_energies(terminal_en),
              terminal_enthalpies(terminal_h) {}

        const std::string unmodified_base;
        const std::string modified_base;
        const std::string fallback_base;
        const std::vector<std::string> pairing_partners;
        const std::map<std::string, float> stacking_energies;
        const std::map<std::string, float> stacking_enthalpies;
        const std::map<std::string, float> terminal_energies;
        const std::map<std::string, float> terminal_enthalpies;
    };

    class ViennaParams {
       public:
        static void load_energy_parameters(const std::string& paramFile = "", const std::string& seq = "") {
            if (!paramFile.empty()) {
                if (file_exists(paramFile)) {
                    int loaded = vrna_params_load(paramFile.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
                    if (!loaded) {
                        THROW_ERROR("Failed to load parameter file: " + paramFile);
                    }
                    std::cout << "Successfully loaded parameter file: " << paramFile << std::endl;
                    return;
                } else {
                    std::cerr << "Warning: Parameter file \"" << paramFile << "\" not found." << std::endl;
                }
            } else {
                std::cerr << "Warning: No parameter file provided. ";
            }

            // Default fallback based on sequence
            if (seq.find('T') != std::string::npos) {
                std::cerr << "Defaulting to DNA parameters (Mathews 2004)." << std::endl;
                vrna_params_load_DNA_Mathews2004();
            } else {
                // std::cerr << "Defaulting to RNA parameters (Turner 2004)." << std::endl;
                // vrna_params_load_RNA_Turner2004();

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
        }

        static modified_base_params load_modified_energy_parameters(const std::string& jsonFile) {
            if (!file_exists(jsonFile)) {
                THROW_ERROR("Modified parameters JSON file \"" + jsonFile + "\" not found.");
            }

            std::vector<std::string> all_files;
            if (is_file(jsonFile)) {
                all_files.push_back(jsonFile);
            } else if (is_directory(jsonFile)) {
                // TODO: implement directory reading
            }

            std::ifstream f(jsonFile);
            json data = json::parse(f);
            json mod = data["modified_base"].get<json>();

            modified_base_params params(
                mod["unmodified"],
                mod["one_letter_code"],
                mod["fallback"],
                mod["pairing_partners"].get<std::vector<std::string>>(),
                mod["stacking_energies"].get<std::map<std::string, float>>(),
                mod["stacking_enthalpies"].get<std::map<std::string, float>>(),
                mod["terminal_energies"].get<std::map<std::string, float>>(),
                mod["terminal_enthalpies"].get<std::map<std::string, float>>()
            );

            return params;
        }
    };
}