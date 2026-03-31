#pragma once

#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

#include "FileUtils.hpp"
#include "Report.hpp"
using json = nlohmann::json;

extern "C" {
// used for load_energy_parameters
#include <ViennaRNA/model.h>
#include <ViennaRNA/params/io.h>
#include <ViennaRNA/utils/basic.h>
}

using string_list = std::vector<std::string>;
using param_map = std::map<std::string, float>;

namespace knotergy {

struct vrna_md_param {
    ~vrna_md_param() {
        if (p) free(p);
    }

    vrna_md_t md{};     ///< ViennaRNA model details.
    vrna_param_t* p{};  ///< ViennaRNA parameters.

    // ------- Details about loading the parameters, for reporting purposes -------

    [[nodiscard]] const ParamSourceInfo& get_source_info() const { return source_info; }
    void set_source_info(const ParamSourceInfo& info) { source_info = info; }

   private:
    ParamSourceInfo source_info;
};

/**
 * @brief Parameters for modified RNA bases.
 *
 * Contains all energy parameters needed to compute thermodynamic contributions
 * of modified nucleotides in RNA secondary structures. Parameters include
 * stacking, terminal, mismatch, and dangle energies.
 */
struct modified_base_param {
    /**
     * @brief Construct modified base parameters.
     *
     * @param mod_name Name of the modified base.
     * @param unmod Unmodified base character it replaces.
     * @param mod Modified base representation.
     * @param fallback Fallback base for energy calculations.
     * @param partners Valid pairing partners.
     * @param stacking Stacking energy parameters.
     * @param enthalpies Stacking enthalpy parameters.
     * @param terminal_e Terminal mismatch energy parameters.
     * @param terminal_h Terminal mismatch enthalpy parameters.
     * @param mismatch_e Mismatch energy parameters.
     * @param mismatch_h Mismatch enthalpy parameters.
     * @param dangle5_e 5' dangle energy parameters.
     * @param dangle5_h 5' dangle enthalpy parameters.
     * @param dangle3_e 3' dangle energy parameters.
     * @param dangle3_h 3' dangle enthalpy parameters.
     */
    modified_base_param(const std::string& mod_name, const std::string& unmod,
                        const std::string& mod, const std::string& fallback,
                        const string_list& partners, const param_map& stacking,
                        const param_map& enthalpies, const param_map& terminal_e,
                        const param_map& terminal_h, const param_map& mismatch_e,
                        const param_map& mismatch_h, const param_map& dangle5_e,
                        const param_map& dangle5_h, const param_map& dangle3_e,
                        const param_map& dangle3_h)
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
          dangle3_enthalpies(dangle3_h) {}

    const std::string name;               ///< Modified base name.
    const std::string unmodified_base;    ///< Unmodified base it replaces.
    const std::string modified_base;      ///< Modified base representation.
    const std::string fallback_base;      ///< Fallback base for calculations.
    const string_list pairing_partners;   ///< Valid pairing partners.
    const param_map stacking_energies;    ///< Stacking energies.
    const param_map stacking_enthalpies;  ///< Stacking enthalpies.
    const param_map terminal_energies;    ///< Terminal mismatch energies.
    const param_map terminal_enthalpies;  ///< Terminal mismatch enthalpies.
    const param_map mismatch_energies;    ///< Mismatch energies.
    const param_map mismatch_enthalpies;  ///< Mismatch enthalpies.
    const param_map dangle5_energies;     ///< 5' dangle energies.
    const param_map dangle5_enthalpies;   ///< 5' dangle enthalpies.
    const param_map dangle3_energies;     ///< 3' dangle energies.
    const param_map dangle3_enthalpies;   ///< 3' dangle enthalpies.
};

/**
 * @brief Manages ViennaRNA energy parameters and model details.
 *
 * This class provides methods to load and access ViennaRNA energy parameters,
 * including support for custom parameter files and modified base parameters.
 */
class ViennaParams {
   public:
    /**
     * @brief Load ViennaRNA energy parameters.
     *
     * Initializes ViennaRNA with energy parameters from a file or defaults.
     *
     * @param paramFile Path to custom parameter file (empty for defaults).
     * @param dangle Dangle model to use (0, 1, 2, or 3). Default is 2.
     * @param seq Optional sequence for parameter initialization.
     */
    static vrna_md_param load_energy_parameters(const std::string& paramFile = "", int dangle = 2,
                                                const std::string& seq = "");

    /**
     * @brief Load modified base energy parameters from a file or directory.
     *
     * @param path Path to a JSON file or directory containing modified base parameter files.
     * @return Vector of modified_base_param structures.
     */
    [[nodiscard]] static std::vector<modified_base_param> load_modified_energy_parameters(
        const std::string& path);

    /**
     * @brief Parse a single modified base JSON parameter file.
     *
     * @param jsonFile Path to JSON file containing modified base parameters.
     * @return modified_base_param structure with parsed parameters.
     */
    [[nodiscard]] static modified_base_param parse_modified_base_json(const std::string& jsonFile);

   private:
    /**
     * @brief Warn if a required JSON key is missing.
     *
     * @param j JSON object to check.
     * @param key Key to look for.
     * @param file File name for error message.
     */
    static void warn_if_missing(const json& j, const std::string& key, const std::string& file);
};

}  // namespace knotergy