#pragma once

#include "io/output/Report.hpp"
#include "utils/FileUtils.hpp"

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <new>
#include <utility>

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
        FileUtils::resolve_data_path("params/ViennaRNA/cache/", check_exists);
    return path;
}

[[nodiscard]] inline const std::string& default_param_path() {
    static const std::string path =
        FileUtils::resolve_data_path("params/ViennaRNA/rna_DirksPierce09.par");
    return path;
}

struct vrna_md_param {
    vrna_md_param() = default;

    // Copy constructor
    vrna_md_param(const vrna_md_param& other)
        : md(other.md), p(clone(other.p)), source_info(other.source_info) {
        std::memcpy(mismatchM_dangles, other.mismatchM_dangles, sizeof(mismatchM_dangles));
    }

    // Move constructor
    vrna_md_param(vrna_md_param&& other) noexcept
        : md(other.md),
          p(std::exchange(other.p, nullptr)),
          source_info(std::move(other.source_info)) {
        std::memcpy(mismatchM_dangles, other.mismatchM_dangles, sizeof(mismatchM_dangles));
    }

    // Copy assignment operator
    vrna_md_param& operator=(vrna_md_param other) noexcept {
        std::swap(md, other.md);
        std::swap(p, other.p);
        std::swap(mismatchM_dangles, other.mismatchM_dangles);
        std::swap(source_info, other.source_info);
        return *this;
    }

    ~vrna_md_param() {
        if (p) free(p);
    }

    viennarna::vrna_md_t md{};     ///< ViennaRNA model details.
    viennarna::vrna_param_t* p{};  ///< ViennaRNA parameters.

    /// ViennaRNA sets p->mismatchM to zero under dangle model 0, but the modified-base corrections
    /// of terminal mismatches in hairpin and internal loops use it as their baseline under every
    /// dangle model for modified bases.
    decltype(viennarna::vrna_param_t::mismatchM) mismatchM_dangles{};

    // Gives mismatchM even in dangle model 0, which ViennaRNA sets to zero. This is the baseline
    // for modified-base corrections of terminal mismatches in hairpin and internal loops.
    [[nodiscard]] int multiloop_mismatch(unsigned int type, int n5, int n3) const {
        return mismatchM_dangles[type][n5][n3];
    }

    // ------- Details about loading the parameters, for reporting purposes -------

    [[nodiscard]] const ParamSourceInfo& get_source_info() const { return source_info; }
    void set_source_info(const ParamSourceInfo& info) { source_info = info; }

   private:
    [[nodiscard]] static viennarna::vrna_param_t* clone(const viennarna::vrna_param_t* source) {
        if (source == nullptr) return nullptr;
        auto* copy =
            static_cast<viennarna::vrna_param_t*>(std::malloc(sizeof(viennarna::vrna_param_t)));
        if (copy == nullptr) throw std::bad_alloc();
        std::memcpy(copy, source, sizeof(viennarna::vrna_param_t));
        return copy;
    }

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
    double temperature = 37.0;
    double salt = 1.021;
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
                                                              const int dangle = 2,
                                                              const double temperature = 37.0,
                                                              const double salt = 1.021,
                                                              const bool disable_cache = false);
};

}  // namespace knotergy