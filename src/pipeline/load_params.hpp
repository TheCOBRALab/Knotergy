#pragma once

#include <iostream>
#include <filesystem>
#include "shared.hpp"

extern "C" {
// used for load_energy_parameters
#include <ViennaRNA/params/io.h>
}

namespace knotergy{
    class ViennaParams {
       public:
        static void load_energy_parameters(const std::string& paramFile = "", const std::string& seq = "") {
            if (!paramFile.empty()) {
                if (std::filesystem::exists(paramFile)) {
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
                std::cerr << "Warning: No parameter file provided." << std::endl;
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

        static void load_modified_energy_parameters();

    };
}