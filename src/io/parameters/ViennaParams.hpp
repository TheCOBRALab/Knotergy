#pragma once

#include "io/output/Report.hpp"
#include "utils/FileUtils.hpp"

#include <fstream>
#include <iostream>

// used for load_energy_parameters
#include <ViennaRNA/model.hpp>
#include <ViennaRNA/params/basic.hpp>
#include <ViennaRNA/params/io.hpp>
#include <ViennaRNA/utils/basic.hpp>

namespace viennarna = thermorna::viennarna;

namespace knotergy {

[[nodiscard]] inline const std::string& cache_path() {
    bool check_exists = false;  // Doesn not need to exist and can be created if it doesn't exist
    static const std::string path =
        FileUtils::resolve_data_path("params/common/cache/", check_exists);
    return path;
}

[[nodiscard]] inline const std::string& default_param_path() {
    static const std::string path =
        FileUtils::resolve_data_path("params/common/rna_DirksPierce09.par");
    return path;
}

struct vrna_md_param {
    ~vrna_md_param() {
        if (p) free(p);
    }

    viennarna::vrna_md_t md{};     ///< ViennaRNA model details.
    viennarna::vrna_param_t* p{};  ///< ViennaRNA parameters.

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
    std::uint32_t param_struct_size = sizeof(viennarna::vrna_param_t);
    std::uint32_t md_struct_size = sizeof(viennarna::vrna_md_t);
    std::uint32_t endian_marker = BigEndianMarker;
    std::int32_t dangles = 2;
    std::uint64_t source_mtime = 0;  // 0 for built-in fallback sets
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
    [[nodiscard]] static vrna_md_param load_energy_parameters(const std::string& paramFile = "",
                                                              int dangle = 2,
                                                              const std::string& seq = "");
};

}  // namespace knotergy