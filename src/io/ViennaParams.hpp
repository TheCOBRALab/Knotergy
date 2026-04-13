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

static const uint32_t BigEndianMarker = 0x01020304;

struct ParamCacheHeader {
    char magic[8] = {'V', 'R', 'N', 'A', 'P', 'R', 'M', '1'};
    std::uint32_t cache_version = 1;
    std::uint32_t param_struct_size = sizeof(vrna_param_t);
    std::uint32_t md_struct_size = sizeof(vrna_md_t);
    std::uint32_t endian_marker = BigEndianMarker;
    std::int32_t dangles = 2;
    std::uint64_t source_mtime = 0;   // 0 for built-in fallback sets
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
    modified_base_param(
        std::string mod_name,
        std::string unmod,
        std::string mod,
        std::string fallback,
        string_list partners,
        param_map stacking,
        param_map enthalpies,
        param_map terminal_e,
        param_map terminal_h,
        param_map mismatch_e,
        param_map mismatch_h,
        param_map dangle5_e,
        param_map dangle5_h,
        param_map dangle3_e,
        param_map dangle3_h)
        : name(std::move(mod_name)),
          unmodified_base(std::move(unmod)),
          modified_base(std::move(mod)),
          fallback_base(std::move(fallback)),
          pairing_partners(std::move(partners)),
          stacking_energies(std::move(stacking)),
          stacking_enthalpies(std::move(enthalpies)),
          terminal_energies(std::move(terminal_e)),
          terminal_enthalpies(std::move(terminal_h)),
          mismatch_energies(std::move(mismatch_e)),
          mismatch_enthalpies(std::move(mismatch_h)),
          dangle5_energies(std::move(dangle5_e)),
          dangle5_enthalpies(std::move(dangle5_h)),
          dangle3_energies(std::move(dangle3_e)),
          dangle3_enthalpies(std::move(dangle3_h)) {}

    const std::string name;
    const std::string unmodified_base;
    const std::string modified_base;
    const std::string fallback_base;
    const string_list pairing_partners;
    const param_map stacking_energies;
    const param_map stacking_enthalpies;
    const param_map terminal_energies;
    const param_map terminal_enthalpies;
    const param_map mismatch_energies;
    const param_map mismatch_enthalpies;
    const param_map dangle5_energies;
    const param_map dangle5_enthalpies;
    const param_map dangle3_energies;
    const param_map dangle3_enthalpies;
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