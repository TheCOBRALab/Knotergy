#pragma once

#include "io/output/Report.hpp"
#include "utils/FileUtils.hpp"

#include <RNAstructure/params/datatable.hpp>
#include <RNAstructure/params/pseudoknot.hpp>

#include <string>
#include <vector>

namespace knotergy {

namespace rnastructure = thermorna::rnastructure;

/// Folder with one sub-folder per RNAstructure parameter set (params/RNAstructure/turner2004,
/// params/RNAstructure/m6A, params/RNAstructure/m1Ψ, ...). See its README.md.
[[nodiscard]] inline std::string default_rnastructure_param_path() {
    return FileUtils::resolve_data_path("params/RNAstructure");
}

/// The set used when none is named.
inline constexpr const char* kDefaultRNAstructureSet = "turner2004";

/**
 * @brief Folder of an RNAstructure parameter set.
 *
 * @p set is the name of a folder in default_rnastructure_param_path() ("turner2004", "dna",
 * "dmem", "m6A", "Ψ", "m1Ψ", "Ψ_full", "m1Ψ_full", or one you added), RNAstructure's name for
 * it (efn2 --alphabet: "rna", "PSU", "1", ...), or a path to a folder. A name containing '/'
 * is always a path.
 */
[[nodiscard]] std::string resolve_rnastructure_set(const std::string& set);

/**
 * @brief efn2's (RNAstructure's) nearest-neighbor parameters, with any number of modified bases.
 *
 * `data` holds every table efn2 uses (stack, helix_ends, terminal mismatches, dangles, coaxial
 * stacking, internal loops, special hairpins, miscloop scalars), in RNAstructure's units of
 * 0.1 kcal/mol. Its accessors take the same arguments as RNAstructure's datatable arrays, e.g.
 * data->tstkh[i][j][i+1][j-1] is data.tstkh(s[i], s[j], s[i+1], s[j-1]) with s = data.encode().
 * `pseudoknot` holds the constants efn2 reads for pseudoknotted structures (pseudconst.dat).
 */
struct rnastructure_param {
    rnastructure::Datatable data;
    rnastructure::PseudoknotConstants pseudoknot;

    [[nodiscard]] const ParamSourceInfo& get_source_info() const { return source_info; }
    void set_source_info(const ParamSourceInfo& info) { source_info = info; }

   private:
    ParamSourceInfo source_info;
};

class RNAstructureParams {
   public:
    /**
     * @brief Load efn2's parameters: a base set plus any number of modified-base sets.
     *
     * @param parameter_set Base set (see resolve_rnastructure_set()); empty: "turner2004".
     *        RNAstructure's data_tables folder also works: efn2's default set "rna" is used.
     * @param modifications Modified-base sets to add, e.g. {"m6A", "Ψ", "m1Ψ"}. When the base
     *        is RNAstructure's data_tables, its own sets are used ("m6A", "PSU", "1").
     * @param temperature Temperature in degrees Celsius.
     * @param keep_37C_without_enthalpies At temperatures other than 37 °C, keep the 37 °C values
     *        of sets without *.dh files (all but turner2004) instead of failing, as efn2 does.
     */
    [[nodiscard]] static rnastructure_param load(const std::string& parameter_set = "",
                                                 const std::vector<std::string>& modifications = {},
                                                 double temperature = 37.0,
                                                 bool keep_37C_without_enthalpies = false);

    /// Load with full control (explicit fallbacks, aliases such as "Ψ"). Folders are paths
    /// (resolve_rnastructure_set() turns names into folders); an empty data_directory means
    /// params/RNAstructure/turner2004. options.temperature is in Kelvin, as in efn2.
    [[nodiscard]] static rnastructure_param load(rnastructure::LoadOptions options);
};

}  // namespace knotergy
