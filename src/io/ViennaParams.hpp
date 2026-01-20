#pragma once

#include <iostream>
#include "FileUtils.hpp"

#include <fstream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

extern "C" {
// used for load_energy_parameters
#include <ViennaRNA/params/io.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/utils/basic.h>
}

namespace knotergy {

    struct modified_base_params {
        modified_base_params(const std::string& mod_name, const std::string& unmod, const std::string& mod,
                             const std::string& fallback,
                             const std::vector<std::string>& partners,
                             const std::map<std::string, float>& stacking,
                             const std::map<std::string, float>& enthalpies,
                             const std::map<std::string, float>& terminal_e,
                             const std::map<std::string, float>& terminal_h,
                             const std::map<std::string, float>& mismatch_e,
                             const std::map<std::string, float>& mismatch_h,
                             const std::map<std::string, float>& dangle5_e,
                             const std::map<std::string, float>& dangle5_h,
                             const std::map<std::string, float>& dangle3_e,
                             const std::map<std::string, float>& dangle3_h)
            : name(mod_name),
              unmodified_base(unmod),
              modified_base(mod),
              fallback_base(fallback),
              pairing_partners(partners),
              stacking_energies(stacking),
              stacking_enthalpies(enthalpies),
              terminal_energies(terminal_e),
              terminal_enthalpies(terminal_h),
              mismatch_energies(mismatch_e),
              mismatch_enthalpies(mismatch_h),
              dangle5_energies(dangle5_e),
              dangle5_enthalpies(dangle5_h),
              dangle3_energies(dangle3_e),
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
        inline static vrna_param_t* p = nullptr; // parameters
        
        ~ViennaParams() {if (p) free(p);}

        static void load_energy_parameters(const std::string& paramFile = "", int dangle = 2, const std::string& seq = "");

        [[nodiscard]] static std::vector<modified_base_params> load_modified_energy_parameters(const std::string& path);

        [[nodiscard]] static modified_base_params parse_modified_base_json(const std::string& jsonFile);
        
        private:
        static void warn_if_missing(const json& j, const std::string& key, const std::string& file);
    };

    
    


} // namespace knotergy