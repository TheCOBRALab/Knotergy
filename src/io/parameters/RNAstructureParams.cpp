#include "RNAstructureParams.hpp"

#include "utils/common.hpp"

#include <RNAstructure/constants.hpp>
#include <RNAstructure/error.hpp>
#include <RNAstructure/params/parameter_set.hpp>

#include <algorithm>
#include <cstddef>
#include <utility>

namespace knotergy {

namespace {

bool is_path(const std::string& set) { return set.find('/') != std::string::npos; }

// pseudconst.dat is in RNAstructure's data_tables, and in params/RNAstructure/pseudoknot.
std::string pseudoknot_folder(const std::string& data_directory) {
    if (FileUtils::is_file(FileUtils::join_path(data_directory, "pseudconst.dat"))) {
        return data_directory;
    }
    return FileUtils::join_path(default_rnastructure_param_path(), "pseudoknot");
}

rnastructure_param load_resolved(const rnastructure::LoadOptions& options, ParamSourceInfo info) {
    rnastructure_param param;
    try {
        param.data = rnastructure::Datatable::load(options);
        param.pseudoknot =
            rnastructure::PseudoknotConstants::read(pseudoknot_folder(options.data_directory));
    } catch (const rnastructure::ParameterError& error) {
        THROW_ERROR(error.what());
    }

    info.label = "RNAstructure";
    info.resolved_path = options.data_directory;
    info.resolved_name = param.data.base().name();
    for (std::size_t i = 0; i < param.data.modification_count(); ++i) {
        info.resolved_name += " + " + param.data.modification_name(i);
    }
    for (const std::string& warning : param.data.warnings()) {
        info.messages.push_back(ParamMessage{warning, true});
    }
    param.set_source_info(info);
    return param;
}

}  // namespace

std::string resolve_rnastructure_set(const std::string& set) {
    if (set.empty()) THROW_ERROR("No RNAstructure parameter set was named.");
    if (is_path(set)) {
        if (!FileUtils::is_directory(set)) {
            THROW_ERROR("RNAstructure parameter folder \"" + set + "\" not found.");
        }
        return set;
    }
    const std::string root = default_rnastructure_param_path();
    try {
        return rnastructure::find_set_folder(root, set).directory;
    } catch (const rnastructure::ParameterError& error) {
        if (FileUtils::is_directory(set)) return set;  // a folder in the working directory
        THROW_ERROR(root.empty() ? "Knotergy's params/RNAstructure folder was not found."
                                 : error.what());
    }
}

rnastructure_param RNAstructureParams::load(const std::string& parameter_set,
                                            const std::vector<std::string>& modifications,
                                            double temperature, bool keep_37C_without_enthalpies) {
    rnastructure::LoadOptions options;
    options.data_directory =
        resolve_rnastructure_set(parameter_set.empty() ? kDefaultRNAstructureSet : parameter_set);
    options.temperature = temperature + rnastructure::kCelsiusToKelvin;
    options.missing_enthalpy = keep_37C_without_enthalpies
                                   ? rnastructure::MissingEnthalpy::KeepFreeEnergy
                                   : rnastructure::MissingEnthalpy::Fail;

    // RNAstructure's data_tables hold every set: modifications it has come from there too.
    std::vector<std::string> shared;
    try {
        shared = rnastructure::list_parameter_sets(options.data_directory);
    } catch (const rnastructure::ParameterError&) {
        // reported when the tables are read
    }
    for (const std::string& name : modifications) {
        rnastructure::Modification modification;
        if (shared.size() > 1 && !is_path(name) &&
            std::find(shared.begin(), shared.end(), name) != shared.end()) {
            modification.alphabet = name;
        } else {
            modification.directory = resolve_rnastructure_set(name);
        }
        options.modifications.push_back(std::move(modification));
    }

    ParamSourceInfo info;
    info.requested_path = parameter_set;
    info.status = parameter_set.empty() ? ParamStatus::Defaulted : ParamStatus::LoadedUserFile;
    return load_resolved(options, std::move(info));
}

rnastructure_param RNAstructureParams::load(rnastructure::LoadOptions options) {
    ParamSourceInfo info;
    info.requested_path = options.data_directory;
    info.status = ParamStatus::LoadedUserFile;
    if (options.data_directory.empty()) {
        options.data_directory = resolve_rnastructure_set(kDefaultRNAstructureSet);
        info.status = ParamStatus::Defaulted;
    }
    return load_resolved(options, std::move(info));
}

}  // namespace knotergy
