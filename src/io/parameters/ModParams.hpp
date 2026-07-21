#pragma once
#include "utils/FileUtils.hpp"

#include <nlohmann/json.hpp>

#include <cstddef>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace knotergy {

using json = nlohmann::json;

[[nodiscard]] inline const std::string& default_mod_param_path() {
    static const std::string path = FileUtils::resolve_data_path("params/modified_bases");
    return path;
}

using string_list = std::vector<std::string>;
using param_map   = std::unordered_map<std::string, float>;

/**
 * @brief Parameters for modified RNA bases.
 *
 * Contains all energy parameters needed to compute thermodynamic contributions
 * of modified nucleotides in RNA secondary structures. Parameters include
 * stacking, terminal, mismatch, and dangle energies.
 */
class modified_base_param {
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
   public:
    modified_base_param(std::string mod_name, std::string unmod, std::string mod,
                        std::string fallback, string_list partners, param_map stacking,
                        param_map enthalpies, param_map terminal_e, param_map terminal_h,
                        param_map mismatch_e, param_map mismatch_h, param_map dangle5_e,
                        param_map dangle5_h, param_map dangle3_e, param_map dangle3_h)
        : name_(std::move(mod_name)),
          unmodified_base_(std::move(unmod)),
          modified_base_(std::move(mod)),
          fallback_base_(std::move(fallback)),
          pairing_partners_(std::move(partners)),
          stacking_energies_(std::move(stacking)),
          stacking_enthalpies_(std::move(enthalpies)),
          terminal_energies_(std::move(terminal_e)),
          terminal_enthalpies_(std::move(terminal_h)),
          mismatch_energies_(std::move(mismatch_e)),
          mismatch_enthalpies_(std::move(mismatch_h)),
          dangle5_energies_(std::move(dangle5_e)),
          dangle5_enthalpies_(std::move(dangle5_h)),
          dangle3_energies_(std::move(dangle3_e)),
          dangle3_enthalpies_(std::move(dangle3_h)) {}

    [[nodiscard]] const std::string& name() const noexcept { return name_; }
    [[nodiscard]] const std::string& unmodified_base() const noexcept { return unmodified_base_; }
    [[nodiscard]] const std::string& modified_base() const noexcept { return modified_base_; }
    [[nodiscard]] const std::string& fallback_base() const noexcept { return fallback_base_; }
    [[nodiscard]] const string_list& pairing_partners() const noexcept { return pairing_partners_; }
    [[nodiscard]] const param_map& stacking_energies() const noexcept { return stacking_energies_; }
    [[nodiscard]] const param_map& stacking_enthalpies() const noexcept {
        return stacking_enthalpies_;
    }
    [[nodiscard]] const param_map& terminal_energies() const noexcept { return terminal_energies_; }
    [[nodiscard]] const param_map& terminal_enthalpies() const noexcept {
        return terminal_enthalpies_;
    }
    [[nodiscard]] const param_map& mismatch_energies() const noexcept { return mismatch_energies_; }
    [[nodiscard]] const param_map& mismatch_enthalpies() const noexcept {
        return mismatch_enthalpies_;
    }
    [[nodiscard]] const param_map& dangle5_energies() const noexcept { return dangle5_energies_; }
    [[nodiscard]] const param_map& dangle5_enthalpies() const noexcept {
        return dangle5_enthalpies_;
    }
    [[nodiscard]] const param_map& dangle3_energies() const noexcept { return dangle3_energies_; }
    [[nodiscard]] const param_map& dangle3_enthalpies() const noexcept {
        return dangle3_enthalpies_;
    }

   private:
    std::string name_;
    std::string unmodified_base_;
    std::string modified_base_;
    std::string fallback_base_;
    string_list pairing_partners_;
    param_map   stacking_energies_;
    param_map   stacking_enthalpies_;
    param_map   terminal_energies_;
    param_map   terminal_enthalpies_;
    param_map   mismatch_energies_;
    param_map   mismatch_enthalpies_;
    param_map   dangle5_energies_;
    param_map   dangle5_enthalpies_;
    param_map   dangle3_energies_;
    param_map   dangle3_enthalpies_;
};

/**
 * @brief Container for all modified base parameters and lookup functionality.
 *
 * This class holds a collection of modified_base_param objects and provides efficient
 * lookup methods to retrieve parameters based on modified base characters, as well as
 * mapping from modified bases to their corresponding unmodified bases. It also includes
 * metadata about the source of the parameters for reporting purposes.
 *
 */
struct all_mod_params {
    all_mod_params() = default;
    all_mod_params(std::vector<modified_base_param> mod_params)
        : mod_params_(std::move(mod_params)) {
        build_lookup();
    }

    [[nodiscard]] const std::vector<modified_base_param>& get_all_params() const {
        return mod_params_;
    }

    [[nodiscard]] bool        empty() const { return mod_params_.empty(); }
    [[nodiscard]] std::size_t size() const { return mod_params_.size(); }

    // Key is the modified base character, returns parameters for that modified base if it exists,
    // otherwise nullptr
    [[nodiscard]] const modified_base_param* get_modified_base_param(
        const std::string& modified_base) const;

    // Key is the modified base character, returns the unmodified base character if it exists,
    // otherwise nullptr
    [[nodiscard]] const std::string* get_unmodified_base(const std::string& modified_base) const;

   private:
    const std::vector<modified_base_param>                      mod_params_;
    std::unordered_map<std::string, const modified_base_param*> mod_param_lookup;
    std::unordered_map<std::string, const std::string*>         mod_to_unmod_lookup;

    void build_lookup();
};

class ModParams {
   public:
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