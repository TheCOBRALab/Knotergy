#include "ViennaParams.hpp"

namespace knotergy {

vrna_md_param ViennaParams::load_energy_parameters(const std::string& paramFile, int dangle,
                                                   const std::string& seq) {
    vrna_md_param md_param{};
    ParamSourceInfo source_info;

    source_info.label = "ViennaRNA";
    source_info.requested_path = paramFile;

    // Find user specified parameter file
    vrna_md_set_default(&md_param.md);
    md_param.md.dangles = dangle;

    // Try to load user-specified parameter file if provided
    if (!paramFile.empty()) {
        if (FileUtils::file_exists(paramFile)) {
            int loaded = vrna_params_load(paramFile.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
            if (!loaded) {
                THROW_ERROR("Failed to load parameter file: " + paramFile);
            }

            md_param.p = vrna_params(&md_param.md);
            source_info.status = ParamStatus::LoadedUserFile;
            source_info.resolved_path = paramFile;
            source_info.resolved_name = FileUtils::get_filename_no_ext(paramFile);

            if (source_info.resolved_name == "rna_DirksPierce09") {
                source_info.resolved_name += "Dirks&Pierce 2009";
            } else if (source_info.resolved_name == "rna_Turner2004") {
                source_info.resolved_name += "Turner 2004";
            } else if (source_info.resolved_name == "dna_Mathews2004") {
                source_info.resolved_name += "Mathews 2004";
            }

            md_param.set_source_info(source_info);
            md_param.p = vrna_params(&md_param.md);

            return md_param;
        }
    }

    // Fallback for DNA
    if (seq.find('T') != std::string::npos) {
        vrna_params_load_DNA_Mathews2004();
        source_info.status = ParamStatus::Fallback;
        source_info.resolved_name = "Mathews 2004 (DNA)";
    }
    // Fallback for RNA (Dirks&Pierce 2009 if available, otherwise Turner 2004)
    else if (FileUtils::file_exists("./params/common/rna_DirksPierce09.par")) {
        int loaded = vrna_params_load("./params/common/rna_DirksPierce09.par",
                                      VRNA_PARAMETER_FORMAT_DEFAULT);
        if (!loaded) {
            THROW_ERROR(
                "Failed to load default RNA parameter file: "
                "./params/common/rna_DirksPierce09.par");
        }

        source_info.status = ParamStatus::Fallback;
        source_info.resolved_path = "./params/common/rna_DirksPierce09.par";
        source_info.resolved_name = "Dirks&Pierce 2009";
    } else {
        vrna_params_load_RNA_Turner2004();
        source_info.status = ParamStatus::Fallback;
        source_info.resolved_name = "Turner 2004";
    }

    md_param.set_source_info(source_info);
    md_param.p = vrna_params(&md_param.md);

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
    json mod = data["modified_base"].get<json>();

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
