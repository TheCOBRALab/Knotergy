#include "ViennaParams.hpp"

#include "io/output/colors.hpp"

#include <ViennaRNA/model.hpp>
#include <ViennaRNA/params/basic.hpp>
#include <ViennaRNA/params/io.hpp>
#include <sys/stat.h>
#include <sys/types.h>

#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <string>

namespace viennarna = thermorna::viennarna;

namespace knotergy {

namespace {
// -------------------------- Helper functions for parameter caching --------------------------

// Generate a cache file path based on the parameter file, dangle model, and sequence type
std::string make_cache_path(const std::string& paramFile, int dangle, const std::string& seq) {
    std::string file_name = FileUtils::strip_extension(paramFile);
    std::string cache_dir = cache_path();
    if (!paramFile.empty()) {
        return cache_dir + file_name + ".d" + std::to_string(dangle) + ".vrna.bin";
    }

    // Fallback cache names
    const bool is_dna = (seq.find('T') != std::string::npos);
    return is_dna ? cache_dir + "dna_Mathews2004.d" + std::to_string(dangle) + ".vrna.bin"
                  : cache_dir + "rna_default.d" + std::to_string(dangle) + ".vrna.bin";
}

// Returns true if cache was successfully loaded and is valid, false otherwise
bool load_param_cache(const std::string& cachePath, int expectedDangle,
                      std::uint64_t expectedSourceMtime, vrna_md_param& out) {
    // Read and validate cache file
    std::ifstream in(cachePath, std::ios::binary);
    if (!in) return false;

    // Read header from cache file
    ParamCacheHeader hdr{};
    in.read(reinterpret_cast<char*>(&hdr), sizeof(hdr));
    if (!in) return false;

    // Validate header fields to ensure cache compatibility
    const ParamCacheHeader expected{};
    if (std::memcmp(hdr.magic, expected.magic, sizeof(hdr.magic)) != 0) return false;
    if (hdr.cache_version != expected.cache_version) return false;
    if (hdr.param_struct_size != sizeof(viennarna::vrna_param_t)) return false;
    if (hdr.md_struct_size != sizeof(viennarna::vrna_md_t)) return false;
    if (hdr.endian_marker != BigEndianMarker) return false;
    if (hdr.dangles != expectedDangle) return false;
    if (hdr.source_mtime != expectedSourceMtime) return false;

    // Allocate memory for vrna_param_t
    viennarna::vrna_param_t* p =
        static_cast<viennarna::vrna_param_t*>(std::malloc(sizeof(viennarna::vrna_param_t)));
    if (!p) return false;

    // Load parameter into the allocated memory
    in.read(reinterpret_cast<char*>(p), sizeof(viennarna::vrna_param_t));
    if (!in) {
        std::free(p);
        return false;
    }

    // Ensure the loaded parameters have the expected dangle model
    if (p->model_details.dangles != expectedDangle) {
        std::free(p);
        return false;
    }

    // update output parameter struct
    out.md = p->model_details;
    out.p = p;
    return true;
}

void save_param_cache(const std::string& cachePath, int dangle, std::uint64_t sourceMtime,
                      const viennarna::vrna_param_t& p) {
    size_t slash = cachePath.find_last_of("/\\");
    if (slash != std::string::npos) {
        std::string dir = cachePath.substr(0, slash);
        if (!dir.empty() && !FileUtils::is_directory(dir)) {
            ::mkdir(dir.c_str(), 0755);
        }
    }

    std::ofstream out(cachePath, std::ios::binary | std::ios::trunc);
    if (!out) return;

    ParamCacheHeader hdr{};
    hdr.dangles = dangle;
    hdr.source_mtime = sourceMtime;

    out.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));
    out.write(reinterpret_cast<const char*>(&p), sizeof(viennarna::vrna_param_t));
}
}  // namespace

//------------------------- Load ViennaRNA Energy Parameters -----------------------

vrna_md_param ViennaParams::load_energy_parameters(const std::string& paramFile, int dangle,
                                                   const std::string& seq,
                                                   const bool disable_cache) {
    vrna_md_param md_param{};
    ParamSourceInfo source_info;

    source_info.label = "ViennaRNA";
    source_info.requested_path = paramFile;

    enum class SourceKind { UserFile, DefaultRNAFile, BuiltinDNA, BuiltinRNA };

    SourceKind source_kind;
    std::string cache_key;
    std::string load_path;
    std::uint64_t srcMtime = 0;

    // ----- Resolve parameter source first -----
    if (!paramFile.empty() && FileUtils::file_exists(paramFile)) {
        // If user provides a param file and it exists, use it.
        source_kind = SourceKind::UserFile;
        cache_key = paramFile;
        load_path = paramFile;
        srcMtime = FileUtils::get_file_mtime(paramFile);

        source_info.status = ParamStatus::LoadedUserFile;
        source_info.resolved_path = paramFile;
        source_info.resolved_name = FileUtils::strip_extension(paramFile);
    } else if (seq.find('T') != std::string::npos) {
        // If user never provided a param file, and the sequence contains T, use DNA parameters.
        source_kind = SourceKind::BuiltinDNA;
        cache_key = "builtin:dna_Mathews2004";

        source_info.status = ParamStatus::Fallback;
        source_info.resolved_name = "Mathews 2004 (DNA)";
    } else if (FileUtils::file_exists(default_param_path())) {
        // If no param file was provided, and the default RNA parameter file exists, use it.
        source_kind = SourceKind::DefaultRNAFile;
        cache_key = default_param_path();
        load_path = default_param_path();
        srcMtime = FileUtils::get_file_mtime(default_param_path());

        source_info.status = ParamStatus::Fallback;
        source_info.resolved_path = default_param_path();
        source_info.resolved_name = "Dirks-Pierce 2009";
    } else {
        // If the default param path was not found, fall back to built-in Turner 2004 parameters.
        // This should never fail since it's hardcoded into ViennaRNA.
        std::cout << WARNING
                  << " Default RNA parameter file not found, falling back to built-in Turner 2004 "
                     "parameters.\n";
        source_kind = SourceKind::BuiltinRNA;
        cache_key = "builtin:rna_Turner2004";

        source_info.status = ParamStatus::Fallback;
        source_info.resolved_name = "Turner 2004";
    }

    // ----- Try cache if not disabled -----
    const std::string cachePath = make_cache_path(cache_key, dangle, seq);
    if (!disable_cache && load_param_cache(cachePath, dangle, srcMtime, md_param)) {
        md_param.set_source_info(source_info);
        return md_param;
    }

    // ----- Cache miss: load actual parameter source -----
    switch (source_kind) {
        case SourceKind::UserFile:
        case SourceKind::DefaultRNAFile: {
            int loaded =
                viennarna::vrna_params_load(load_path.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
            if (!loaded) {
                THROW_ERROR("Failed to load parameter file: " + load_path);
            }
            break;
        }

        case SourceKind::BuiltinDNA:
            if (!viennarna::vrna_params_load_DNA_Mathews2004()) {
                THROW_ERROR("Failed to load built-in DNA parameters");
            }
            break;

        case SourceKind::BuiltinRNA:
            if (!viennarna::vrna_params_load_RNA_Turner2004()) {
                THROW_ERROR("Failed to load built-in RNA parameters");
            }
            break;
    }

    // ----- Initialize md after loading params -----
    viennarna::vrna_md_set_default(&md_param.md);
    md_param.md.dangles = dangle;

    md_param.p = viennarna::vrna_params(&md_param.md);
    if (!md_param.p) {
        THROW_ERROR("Failed to create scaled ViennaRNA parameter table");
    }

    md_param.set_source_info(source_info);

    // ----- Save cache with the loaded parameters -----
    if (!disable_cache) {
        save_param_cache(cachePath, dangle, srcMtime, *md_param.p);
    }

    return md_param;
}

}  // namespace knotergy
