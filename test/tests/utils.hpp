#pragma once

#include <gtest/gtest.h>

#include <io/common.hpp>
#include <io/RNAInputManager.hpp>
#include <preprocessing/RNAEntry.hpp>
#include <preprocessing/ProcessedRNAEntry.hpp>
#include <preprocessing/RNAProcessor.hpp>
#include <loop_tree/LoopFactory.hpp>
#include <energy/ComputeEnergy.hpp>

#include <string>
#include <vector>

const std::string DP_file = std::string(KNOTERGY_SOURCE_DIR) + "/params/common/rna_DirksPierce09.par";
const std::string turner_file = std::string(KNOTERGY_SOURCE_DIR) + "/params/common/rna_turner2004.par";
const std::string pkp_file = std::string(KNOTERGY_SOURCE_DIR) + "/params/pseudo/rna_pk_DirksPierce09_HotKnotsV2.json";
const std::string mod_folder = std::string(KNOTERGY_SOURCE_DIR) + "/params/modified_bases";

inline static float get_energy(std::string sequence, std::string structure, int dangle = 2, bool round = false,
    std::string param_file = turner_file,
    std::string pseudoknot_param_file = pkp_file,
    std::string mod_param_file = mod_folder) {
    
    // Load parameters
    knotergy::vrna_md_param vp = knotergy::ViennaParams::load_energy_parameters(param_file, dangle, sequence);
    knotergy::pk_param pkp = knotergy::PseudoknotParams::load_pk_param(pseudoknot_param_file);
    std::vector<knotergy::modified_base_param> mod_params = knotergy::ViennaParams::load_modified_energy_parameters(mod_param_file);
    
    // pre-process RNA entry
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna), mod_params));

    // build loop tree
    knotergy::LoopFactory factory(processed_rna);

    // compute energy
    knotergy::ComputeEnergy energy(factory.get_root_node(), processed_rna, vp, pkp, mod_params, round);

    return energy.getEnergy();
}

inline static float get_energy(knotergy::RNAEntry rna, int dangle = 2, bool round = false,
    std::string param_file = turner_file,
    std::string pseudoknot_param_file = pkp_file,
    std::string mod_param_file = mod_folder) {
    return get_energy(rna.sequence, rna.structure, dangle, round, param_file, pseudoknot_param_file, mod_param_file);
}

inline static float get_energy(std::string sequence, std::string structure, 
    std::string param_file = turner_file) {
    const int dangle = 2;
    const bool round = false;
    return get_energy(sequence, structure, dangle, round, param_file);
}

