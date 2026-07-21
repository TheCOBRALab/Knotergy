#pragma once

#include <energy/ComputeEnergy.hpp>
#include <gtest/gtest.h>
#include <io/input/RNAInputManager.hpp>
#include <io/parameters/ModParams.hpp>
#include <loop_tree/LoopFactory.hpp>
#include <preprocessing/ProcessedRNAEntry.hpp>
#include <preprocessing/RNAEntry.hpp>
#include <preprocessing/RNAProcessor.hpp>
#include <utils/common.hpp>

#include <string>
#include <vector>

const std::string DP_file =
    knotergy::FileUtils::resolve_data_path("params/common/rna_DirksPierce09.par");
const std::string turner_file =
    knotergy::FileUtils::resolve_data_path("params/common/rna_turner2004.par");
const std::string pkp_file =
    knotergy::FileUtils::resolve_data_path("params/pseudo/rna_pk_DirksPierce09.json");
const std::string mod_folder = knotergy::FileUtils::resolve_data_path("params/modified_bases");

inline static double get_energy(std::string sequence, std::string structure, int dangle = 2,
                                bool round = false, std::string param_file = turner_file,
                                std::string pseudoknot_param_file = pkp_file,
                                std::string mod_param_file        = mod_folder) {
    // Load parameters
    knotergy::vrna_md_param vp =
        knotergy::ViennaParams::load_energy_parameters(param_file, dangle, sequence);
    knotergy::pk_param pkp = knotergy::PseudoknotParams::load_pk_param(pseudoknot_param_file);
    std::vector<knotergy::modified_base_param> mod_params =
        knotergy::ModParams::load_modified_energy_parameters(mod_param_file);

    knotergy::all_mod_params mp{mod_params};

    // pre-process RNA entry
    knotergy::RNAEntry          rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(
        knotergy::RNAProcessor::process_rna(std::move(rna), mp));

    // build loop tree
    knotergy::LoopFactory factory(processed_rna, vp);

    // compute energy
    knotergy::ComputeEnergy energy(factory.get_root_node(), processed_rna, vp, pkp, mp, round);

    return energy.getEnergy();
}

inline static double get_energy(knotergy::RNAEntry rna, int dangle = 2, bool round = false,
                                std::string param_file            = turner_file,
                                std::string pseudoknot_param_file = pkp_file,
                                std::string mod_param_file        = mod_folder) {
    return get_energy(rna.sequence, rna.structure, dangle, round, param_file, pseudoknot_param_file,
                      mod_param_file);
}

inline static double get_energy(std::string sequence, std::string structure,
                                std::string param_file = turner_file) {
    const int  dangle = 2;
    const bool round  = false;
    return get_energy(sequence, structure, dangle, round, param_file);
}
