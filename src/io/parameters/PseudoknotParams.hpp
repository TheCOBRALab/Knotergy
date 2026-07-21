#pragma once

#include "io/output/Report.hpp"
#include "utils/FileUtils.hpp"
#include "utils/colors.hpp"
#include "utils/common.hpp"

#include <nlohmann/json.hpp>

#include <fstream>
#include <iostream>
#include <string>

using json = nlohmann::json;

namespace knotergy {

[[nodiscard]] inline const std::string& default_pk_param_path() {
    static const std::string path =
        FileUtils::resolve_data_path("params/pseudo/rna_pk_DirksPierce09.json");
    return path;
}

/**
 * @brief Pseudoknot energy parameters.
 *
 * Contains all penalty parameters for computing pseudoknot energies, including
 * initialization penalties, band penalties, and multipliers for different loop types.
 */
struct pk_param {
    /**
     * @brief Construct with default DirksPierce09 parameters from HotKnotsV2.
     *
     * Hard coded values based on the original HotKnotsV2 implementation
     */
    pk_param()
        : name("DirksPierce09 (Hard-coded default)"),
          pk_in_ext(-138),
          pk_in_mloop(1007),
          pk_in_pk(1500),
          band(246),
          unpaired_in_pk(6),
          cr_in_pk(96),
          pk_stack_x(0.89),
          pk_internal_x(0.74),
          pk_mloop_init(341),
          pk_mloop_bp(56),
          pk_mloop_unpaired(12) {}

    /**
     * @brief Construct with custom pseudoknot parameters.
     *
     * @param param_name Name of the parameter set.
     * @param pk_ext Pseudoknot in exterior loop penalty.
     * @param pk_multi Pseudoknot in multiloop penalty.
     * @param pk_pk Pseudoknot in pseudoloop penalty.
     * @param band_pen Band penalty.
     * @param unpaired_pk Unpaired bases in pseudoknot penalty.
     * @param cr_in_pk_pen Closed region nested in pseudoknot penalty.
     * @param pk_stack_multiplier Stacked pair spanning band multiplier.
     * @param pk_internal_multiplier Internal pair spanning band multiplier.
     * @param pk_mloop_init_pen Multiloop spanning band initialization penalty.
     * @param pk_mloop_bp_pen Base pair in multiloop spanning band penalty.
     * @param pk_mloop_unpaired_pen Unpaired bases in multiloop spanning band penalty.
     */
    pk_param(const std::string& param_name, int pk_ext, int pk_multi, int pk_pk, int band_pen,
             int unpaired_pk, int cr_in_pk_pen, double pk_stack_multiplier,
             double pk_internal_multiplier, int pk_mloop_init_pen, int pk_mloop_bp_pen,
             int pk_mloop_unpaired_pen)
        : name(param_name),
          pk_in_ext(pk_ext),
          pk_in_mloop(pk_multi),
          pk_in_pk(pk_pk),
          band(band_pen),
          unpaired_in_pk(unpaired_pk),
          cr_in_pk(cr_in_pk_pen),
          pk_stack_x(pk_stack_multiplier),
          pk_internal_x(pk_internal_multiplier),
          pk_mloop_init(pk_mloop_init_pen),
          pk_mloop_bp(pk_mloop_bp_pen),
          pk_mloop_unpaired(pk_mloop_unpaired_pen) {}

    const std::string name;               ///< Parameter set name.
    const int         pk_in_ext;          ///< Pseudoknot in exterior loop penalty.
    const int         pk_in_mloop;        ///< Pseudoknot in multiloop penalty.
    const int         pk_in_pk;           ///< Pseudoknot in pseudoloop penalty.
    const int         band;               ///< Band penalty.
    const int         unpaired_in_pk;     ///< Unpaired bases in pseudoknot penalty.
    const int         cr_in_pk;           ///< Closed region nested in pseudoknot penalty.
    const double      pk_stack_x;         ///< Stacked pair spanning band multiplier.
    const double      pk_internal_x;      ///< Internal pair spanning band multiplier.
    const int         pk_mloop_init;      ///< Multiloop spanning band initialization penalty.
    const int         pk_mloop_bp;        ///< Base pair in multiloop spanning band penalty.
    const int         pk_mloop_unpaired;  ///< Unpaired bases in multiloop spanning band penalty.

    // ------- Details about loading the parameters, for reporting purposes -------

    [[nodiscard]] const ParamSourceInfo& get_source_info() const { return source_info; }

    void set_source_info(const ParamSourceInfo& info) { source_info = info; }

   private:
    ParamSourceInfo source_info;
};

/**
 * @brief Manages pseudoknot energy parameters.
 *
 * Provides static methods to load pseudoknot parameters from JSON files
 * or use default parameters.
 */
class PseudoknotParams {
   public:
    /**
     * @brief Load pseudoknot parameters from a JSON file.
     *
     * @param paramFile Path to JSON file with pseudoknot parameters.
     * @return Loaded pk_param structure.
     * @throws DetailedException if file not found or invalid.
     */
    [[nodiscard]] static pk_param load_pk_param(
        const std::string& paramFile = default_pk_param_path()) {
        ParamSourceInfo info;
        info.label          = "Pseudoknot";
        info.requested_path = paramFile;

        if (paramFile.empty()) {
            // Prioritize using the default parameter file if it exists,
            // otherwise use hard-coded defaults.
            bool file_exists = FileUtils::file_exists(default_pk_param_path());
            if (!file_exists) {
                std::cout << WARNING
                          << " Warning: Default pseudoknot parameter file not found. Using "
                             "hard-coded defaults.\n";
            }

            pk_param pkp = file_exists ? parse_pk_json(default_pk_param_path()) : pk_param();
            info.status  = ParamStatus::Defaulted;
            pkp.set_source_info(info);  // Set the source info for reporting purposes
            return pkp;
        }

        if (!FileUtils::file_exists(paramFile)) {
            THROW_ERROR("Pseudoknot parameters JSON file \"" + paramFile + "\" not found.");
        }

        if (!FileUtils::is_file(paramFile)) {
            THROW_ERROR("Pseudoknot parameters path \"" + paramFile + "\" is not a file.");
        }

        pk_param pkp = parse_pk_json(paramFile);

        info.status        = ParamStatus::LoadedUserFile;
        info.resolved_path = paramFile;
        info.resolved_name = pkp.name;

        pkp.set_source_info(info);

        return pkp;
    }

    /**
     * @brief Parse pseudoknot parameters from a JSON file.
     *
     * @param jsonFile Path to JSON file.
     * @return Parsed pk_param structure.
     * @throws DetailedException if file cannot be opened or parsed.
     */
    [[nodiscard]] static pk_param parse_pk_json(const std::string& jsonFile) {
        std::ifstream f(jsonFile);
        if (!f.is_open()) {
            THROW_ERROR("Error: Unable to open pseudoknot parameter file: " + jsonFile);
        }

        json data = json::parse(f);

        auto it = data.find("pseudoknot_parameters");
        if (it == data.end()) {
            THROW_ERROR("Missing required key 'pseudoknot_parameters' in file " + jsonFile);
        }
        const auto& pk = *it;

        return pk_param(pk.value("name", std::string{"Nameless"}), pk.at("pk_in_ext").get<int>(),
                        pk.at("pk_in_mloop").get<int>(), pk.at("pk_in_pk").get<int>(),
                        pk.at("band").get<int>(), pk.at("unpaired_in_pk").get<int>(),
                        pk.at("cr_in_pk").get<int>(), pk.at("pk_stack_x").get<double>(),
                        pk.at("pk_internal_x").get<double>(), pk.at("pk_mloop_init").get<int>(),
                        pk.at("pk_mloop_bp").get<int>(), pk.at("pk_mloop_unpaired").get<int>());
    }
};

}  // namespace knotergy