#include "ViennaParams.hpp"

#include <ViennaRNA/model.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/params/io.h>

#include <cstdint>
#include <cstring>
#include <fstream>
#include <memory>
#include <string>
#include <cstdlib> 
#include <sys/stat.h>
#include <sys/types.h>
#include <cerrno>


namespace knotergy {

namespace {
// -------------------------- Helper functions for parameter caching --------------------------

// Generate a cache file path based on the parameter file, dangle model, and sequence type
std::string make_cache_path(const std::string& paramFile, int dangle, const std::string& seq) {
    std::string file_name = FileUtils::get_filename_no_ext(paramFile);
    std::string cache_dir = std::string(KNOTERGY_SOURCE_DIR) + "/params/common/cache/";
    if (!paramFile.empty()) {
        return cache_dir + file_name + ".d" + std::to_string(dangle) + ".vrna.bin";
    }

    // Fallback cache names
    const bool is_dna = (seq.find('T') != std::string::npos);
    return is_dna
        ? cache_dir + "dna_Mathews2004.d" + std::to_string(dangle) + ".vrna.bin"
        : cache_dir + "rna_default.d" + std::to_string(dangle) + ".vrna.bin";
}

// Returns true if cache was successfully loaded and is valid, false otherwise
bool load_param_cache(const std::string& cachePath,
                      int expectedDangle,
                      std::uint64_t expectedSourceMtime,
                      vrna_md_param& out) {
                        
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
    if (hdr.param_struct_size != sizeof(vrna_param_t)) return false;
    if (hdr.md_struct_size != sizeof(vrna_md_t)) return false;
    if (hdr.endian_marker != BigEndianMarker) return false;
    if (hdr.dangles != expectedDangle) return false;
    if (hdr.source_mtime != expectedSourceMtime) return false;
    
    // Allocate memory for vrna_param_t
    vrna_param_t *p = static_cast<vrna_param_t*>(std::malloc(sizeof(vrna_param_t)));
    if (!p) return false;
    
    // Load parameter into the allocated memory
    in.read(reinterpret_cast<char*>(p), sizeof(vrna_param_t));
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
    out.p  = p;
    return true;
}

void save_param_cache(const std::string& cachePath,
                      int dangle,
                      std::uint64_t sourceMtime,
                      const vrna_param_t& p) {
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
    out.write(reinterpret_cast<const char*>(&p), sizeof(vrna_param_t));
}
} // namespace


//------------------------- Load ViennaRNA Energy Parameters -----------------------

vrna_md_param ViennaParams::load_energy_parameters(const std::string& paramFile,
                                                   int dangle,
                                                   const std::string& seq) {
    const std::string DP_path_str =
        std::string(KNOTERGY_SOURCE_DIR) + "/params/common/rna_DirksPierce09.par";

    vrna_md_param md_param{};
    ParamSourceInfo source_info;

    source_info.label = "ViennaRNA";
    source_info.requested_path = paramFile;

    enum class SourceKind {
        UserFile,
        DefaultRNAFile,
        BuiltinDNA,
        BuiltinRNA
    };

    SourceKind source_kind;
    std::string cache_key;
    std::string load_path;
    std::uint64_t srcMtime = 0;

    // ----- Resolve parameter source first -----
    if (!paramFile.empty() && FileUtils::file_exists(paramFile)) {
        source_kind = SourceKind::UserFile;
        cache_key = paramFile;
        load_path = paramFile;
        srcMtime = FileUtils::get_file_mtime(paramFile);

        source_info.status = ParamStatus::LoadedUserFile;
        source_info.resolved_path = paramFile;
        source_info.resolved_name =FileUtils::get_filename_no_ext(paramFile);
    } else if (seq.find('T') != std::string::npos) {
        source_kind = SourceKind::BuiltinDNA;
        cache_key = "builtin:dna_Mathews2004";

        source_info.status = ParamStatus::Fallback;
        source_info.resolved_name = "Mathews 2004 (DNA)";
    } else if (FileUtils::file_exists(DP_path_str)) {
        source_kind = SourceKind::DefaultRNAFile;
        cache_key = DP_path_str;
        load_path = DP_path_str;
        srcMtime = FileUtils::get_file_mtime(DP_path_str);

        source_info.status = ParamStatus::Fallback;
        source_info.resolved_path = DP_path_str;
        source_info.resolved_name = "Dirks&Pierce 2009";
    } else {
        source_kind = SourceKind::BuiltinRNA;
        cache_key = "builtin:rna_Turner2004";

        source_info.status = ParamStatus::Fallback;
        source_info.resolved_name = "Turner 2004";
    }

    // ----- Always try cache -----
    const std::string cachePath = make_cache_path(cache_key, dangle, seq);
    if (load_param_cache(cachePath, dangle, srcMtime, md_param)) {
        md_param.set_source_info(source_info);
        return md_param;
    }

    // ----- Cache miss: load actual parameter source -----
    switch (source_kind) {
        case SourceKind::UserFile:
        case SourceKind::DefaultRNAFile: {
            int loaded =
                vrna_params_load(load_path.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
            if (!loaded) {
                THROW_ERROR("Failed to load parameter file: " + load_path);
            }
            break;
        }

        case SourceKind::BuiltinDNA:
            if (!vrna_params_load_DNA_Mathews2004()) {
                THROW_ERROR("Failed to load built-in DNA parameters");
            }
            break;

        case SourceKind::BuiltinRNA:
            if (!vrna_params_load_RNA_Turner2004()) {
                THROW_ERROR("Failed to load built-in RNA parameters");
            }
            break;
    }

    // ----- Initialize md after loading params -----
    vrna_md_set_default(&md_param.md);
    md_param.md.dangles = dangle;

    md_param.p = vrna_params(&md_param.md);
    if (!md_param.p) {
        THROW_ERROR("Failed to create scaled ViennaRNA parameter table");
    }

    md_param.set_source_info(source_info);

    // ----- Always save cache -----
    save_param_cache(cachePath, dangle, srcMtime, *md_param.p);

    return md_param;
}

std::vector<modified_base_param> ViennaParams::load_modified_energy_parameters(
    const std::string& path) {
    if (path.empty()) {
        return {};
    }

    if (!FileUtils::file_exists(path)) {
        THROW_ERROR("Modified parameters JSON file \"" + path + "\" not found.");
    }

    string_list all_files;
    if (FileUtils::is_file(path)) {
        all_files.push_back(path);
    } else if (FileUtils::is_directory(path)) {
        all_files = FileUtils::get_files_in_dir(path);
    }

    std::vector<modified_base_param> params_list;

    for (const std::string& file : all_files) {
        if (file.size() >= 5 && file.substr(file.size() - 5) == ".json") {
            params_list.push_back(parse_modified_base_json(file));
        }
    }

    return params_list;
}

modified_base_param ViennaParams::parse_modified_base_json(const std::string& jsonFile) {
    std::ifstream f(jsonFile);
    if (!f.is_open()) {
        THROW_ERROR("Error: Unable to open modified base parameter file: " + jsonFile);
    }
    json data = json::parse(f);
    const json& mod = data.at("modified_base");

    warn_if_missing(mod, "name", jsonFile);
    warn_if_missing(mod, "unmodified", jsonFile);
    warn_if_missing(mod, "one_letter_code", jsonFile);
    warn_if_missing(mod, "fallback", jsonFile);
    warn_if_missing(mod, "pairing_partners", jsonFile);

    modified_base_param params(
        mod.value("name", ""), mod.value("unmodified", ""), mod.value("one_letter_code", ""),
        mod.value("fallback", ""), mod.value("pairing_partners", string_list{}),
        mod.value("stacking_energies", param_map{}), mod.value("stacking_enthalpies", param_map{}),
        mod.value("terminal_energies", param_map{}), mod.value("terminal_enthalpies", param_map{}),
        mod.value("mismatch_energies", param_map{}), mod.value("mismatch_enthalpies", param_map{}),
        mod.value("dangle5_energies", param_map{}), mod.value("dangle5_enthalpies", param_map{}),
        mod.value("dangle3_energies", param_map{}), mod.value("dangle3_enthalpies", param_map{}));

    return params;
}

void ViennaParams::warn_if_missing(const json& j, const std::string& key, const std::string& file) {
    if (!j.contains(key)) {
        std::cerr << "Warning: modified_base missing required key '" << key << "' in file " << file
                  << std::endl;
    }
}

}  // namespace knotergy
