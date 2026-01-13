#pragma once

#include <iostream>
#include "shared.hpp"

#include <fstream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

extern "C" {
// used for load_energy_parameters
#include <ViennaRNA/params/io.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/utils/basic.h>
}

namespace knotergy{
    struct modified_base_params {

        modified_base_params(const std::string& mod_name, const std::string& unmod, const std::string& mod,
                             const std::string& fallback,
                             const std::vector<std::string>& partners,
                             const std::map<std::string, float>& stacking,
                             const std::map<std::string, float>& enthalpies,
                             const std::map<std::string, float>& terminal_en,
                             const std::map<std::string, float>& terminal_h,
                             const std::map<std::string, float>& mismatch_en,
                             const std::map<std::string, float>& mismatch_h,
                             const std::map<std::string, float>& dangle5_en,
                             const std::map<std::string, float>& dangle5_h,
                             const std::map<std::string, float>& dangle3_en,
                             const std::map<std::string, float>& dangle3_h)
            : name(mod_name),
              unmodified_base(unmod),
              modified_base(mod),
              fallback_base(fallback),
              pairing_partners(partners),
              stacking_energies(stacking),
              stacking_enthalpies(enthalpies),
              terminal_energies(terminal_en),
              terminal_enthalpies(terminal_h),
              mismatch_energies(mismatch_en),
              mismatch_enthalpies(mismatch_h),
              dangle5_energies(dangle5_en),
              dangle5_enthalpies(dangle5_h),
              dangle3_energies(dangle3_en),
              dangle3_enthalpies(dangle3_h){}

        const std::string name;
        const std::string unmodified_base;
        const std::string modified_base;
        const std::string fallback_base;
        const std::vector<std::string> pairing_partners;
        const std::map<std::string, float> stacking_energies;
        const std::map<std::string, float> stacking_enthalpies;
        const std::map<std::string, float> terminal_energies;
        const std::map<std::string, float> terminal_enthalpies;
        const std::map<std::string, float> mismatch_energies;
        const std::map<std::string, float> mismatch_enthalpies;
        const std::map<std::string, float> dangle5_energies;
        const std::map<std::string, float> dangle5_enthalpies;
        const std::map<std::string, float> dangle3_energies;
        const std::map<std::string, float> dangle3_enthalpies;
    };

    class ViennaParams {
       public:
        inline static vrna_md_t md{};    // model details
        inline static vrna_param_t* P = nullptr; // parameters

        ~ViennaParams() {
            if (P) free(P);
        }

        static void load_energy_parameters(const std::string& paramFile = "", int dangle = 2, const std::string& seq = "") { //  (detecting DNA disabled)
            vrna_md_set_default(&ViennaParams::md);
            ViennaParams::md.dangles = dangle;
            if (!paramFile.empty()) {
                if (file_exists(paramFile)) {
                    int loaded = vrna_params_load(paramFile.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
                    if (!loaded) {
                        THROW_ERROR("Failed to load parameter file: " + paramFile);
                    }
                    std::cout << "Successfully loaded parameter file: " << paramFile << std::endl;
                    ViennaParams::P = vrna_params(&ViennaParams::md);
                    return;
                } else {
                    std::cerr << "Warning: Parameter file \"" << paramFile << "\" not found." << std::endl;
                }
            } else {
                std::cerr << "No parameter file provided. ";
            }

            // Default fallback based on sequence
            if (seq.find('T') != std::string::npos) { //  (detect DNA but currently disabled)
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
            ViennaParams::P = vrna_params(&ViennaParams::md);
            return;
        }

        static std::vector<modified_base_params> load_modified_energy_parameters(const std::string& jsonFile) {
            if (jsonFile.empty()) {
                return {};
            }
            
            if (!file_exists(jsonFile)) {
                THROW_ERROR("Modified parameters JSON file \"" + jsonFile + "\" not found.");
            }

            std::vector<std::string> all_files;
            if (is_file(jsonFile)) {
                all_files.push_back(jsonFile);
            } else if (is_directory(jsonFile)) {
                all_files = list_files_in_dir(jsonFile);
            }

            std::vector<modified_base_params> params_list;

            for (const std::string& file : all_files) {
                if (file.size() >= 5 && file.substr(file.size() - 5) == ".json") {
                    params_list.push_back(parse_modified_base_json(file));
                }
            }

            return params_list;
        }


        static modified_base_params parse_modified_base_json(const std::string& jsonFile) {
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

            modified_base_params params(
                mod.value("name", ""),
                mod.value("unmodified", ""),
                mod.value("one_letter_code", ""),
                mod.value("fallback", ""),
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
                mod.value("dangle3_enthalpies", std::map<std::string, float>{})
            );

            return params;
        }

        static std::vector<std::string> get_supported_modified_bases(const std::vector<modified_base_params>& params_list) {
            std::vector<std::string> supported_bases;
            for (const auto& params : params_list) {
                supported_bases.push_back(params.modified_base);
            }
            return supported_bases;
        }

    private:
        static void warn_if_missing(const json& j, const std::string& key, const std::string& file)
        {
            if (!j.contains(key)) {
                std::cerr << "Warning: modified_base missing required key '"
                          << key << "' in file " << file << std::endl;
            }
        }
    };
}