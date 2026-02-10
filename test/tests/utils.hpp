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


static float get_energy(std::string sequence, std::string structure, int dangle, bool round = false,
    std::string param_file = "../../params/common/rna_turner2004.par",
    std::string pseudoknot_param_file = "../../params/pseudo/rna_pk_DirksPierce09_HotKnotsV2.json",
    std::string mod_param_file = "../../params/modified/") {
    
    // Load parameters
    knotergy::ViennaParams::load_energy_parameters(param_file, dangle, sequence);
    knotergy::pk_param pseudoknot_param = knotergy::PseudoknotParams::load_pk_param(pseudoknot_param_file);
    std::vector<knotergy::modified_base_param> mod_params = knotergy::ViennaParams::load_modified_energy_parameters(mod_param_file);
    
    // pre-process RNA entry
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));

    // build loop tree
    knotergy::LoopFactory factory(processed_rna);

    // compute energy
    knotergy::ComputeEnergy energy(factory.get_root_node(), processed_rna, mod_params, dangle, round);

    return energy.getEnergy();
}

static float get_energy(knotergy::RNAEntry rna, int dangle, bool round = false,
    std::string param_file = "../../params/common/rna_turner2004.par",
    std::string pseudoknot_param_file = "../../params/pseudo/rna_pk_DirksPierce09_HotKnotsV2.json",
    std::string mod_param_file = "../../params/modified/") {
    return get_energy(rna.sequence, rna.structure, dangle, round, param_file, pseudoknot_param_file, mod_param_file);
}



