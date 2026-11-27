#include "energy/EnergyFunctions.hpp"

#include "energy/modified_bases/ModBaseUtils.hpp"
#include "energy/modified_bases/ModExternal.hpp"
#include "energy/modified_bases/ModHairpin.hpp"
#include "energy/modified_bases/ModInternal.hpp"
#include "energy/modified_bases/ModMultiloop.hpp"
#include "energy/modified_bases/ModStack.hpp"
#include "energy/rna_structure/RNAstructureFunctions.hpp"
#include "energy/vienna/ViennaFunctions.hpp"
#include "utils/common.hpp"

namespace knotergy {
int EnergyFunctions::stack(const std::size_t i, const std::size_t j, const std::size_t ci,
                           const std::size_t cj, const std::string& sequence,
                           const std::vector<std::string_view>& mod_sequence) {
    if (em_ == PKFreeEnergyModel::efn2 || em_ == PKFreeEnergyModel::enf2Simple) {
        return RNAstructureFunctions::stack_energy(i, j, ci, cj, sequence, mod_sequence, rsp_);
    }

    if (!mp_.empty()) {
        return ModStack::find_mod_stack_energy(i, j, ci, cj, sequence, mod_sequence, vp_, mp_);
    }

    if (em_ == PKFreeEnergyModel::ViennaRNA2) {
        return ViennaFunctions::stack_energy(i, j, ci, cj, sequence, vp_);
    }

    return NULL_ENERGY;
}

int EnergyFunctions::d5(const std::size_t i, const std::size_t j, const std::string& sequence,
                        const std::vector<std::string_view>& mod_sequence, const bool is_closing) {
    if (em_ == PKFreeEnergyModel::efn2 || em_ == PKFreeEnergyModel::enf2Simple) {
        return RNAstructureFunctions::dangle5_energy(i, j, sequence, mod_sequence, rsp_,
                                                     is_closing);
    }

    if (!mp_.empty()) {
        // The JSON files only have the dangles that involve a modified base.
        const int energy = ModBaseUtils::danlge5_mod_energy_lookup(i, j, mod_sequence, mp_,
                                                                   vp_.md.temperature, is_closing);
        return energy != NULL_ENERGY
                   ? energy
                   : ViennaFunctions::dangle5_energy(i, j, sequence, vp_, is_closing);
    }

    if (em_ == PKFreeEnergyModel::ViennaRNA2) {
        return ViennaFunctions::dangle5_energy(i, j, sequence, vp_, is_closing);
    }

    return NULL_ENERGY;
}

int EnergyFunctions::d3(const std::size_t i, const std::size_t j, const std::string& sequence,
                        const std::vector<std::string_view>& mod_sequence, const bool is_closing) {
    if (em_ == PKFreeEnergyModel::efn2 || em_ == PKFreeEnergyModel::enf2Simple) {
        return RNAstructureFunctions::dangle3_energy(i, j, sequence, mod_sequence, rsp_,
                                                     is_closing);
    }

    if (!mp_.empty()) {
        // The JSON files only have the dangles that involve a modified base.
        const int energy = ModBaseUtils::danlge3_mod_energy_lookup(i, j, mod_sequence, mp_,
                                                                   vp_.md.temperature, is_closing);
        return energy != NULL_ENERGY
                   ? energy
                   : ViennaFunctions::dangle3_energy(i, j, sequence, vp_, is_closing);
    }

    if (em_ == PKFreeEnergyModel::ViennaRNA2) {
        return ViennaFunctions::dangle3_energy(i, j, sequence, vp_, is_closing);
    }

    return NULL_ENERGY;
}

int EnergyFunctions::mismatch(const std::size_t i, const std::size_t j, const std::string& sequence,
                              const std::vector<std::string_view>& mod_sequence,
                              const bool is_external, const bool is_closing) {
    if (em_ == PKFreeEnergyModel::efn2 || em_ == PKFreeEnergyModel::enf2Simple) {
        return RNAstructureFunctions::mismatch_energy(i, j, sequence, mod_sequence, rsp_,
                                                      is_external, is_closing);
    }

    if (!mp_.empty()) {
        // The JSON files only have the mismatches that involve a modified base.
        const int energy = ModBaseUtils::get_mismatch_mod_energy(i, j, mod_sequence, mp_,
                                                                 vp_.md.temperature, is_closing);
        return energy != NULL_ENERGY
                   ? energy
                   : ViennaFunctions::mismatch_energy(i, j, sequence, vp_, is_external, is_closing);
    }

    if (em_ == PKFreeEnergyModel::ViennaRNA2) {
        return ViennaFunctions::mismatch_energy(i, j, sequence, vp_, is_external, is_closing);
    }

    return NULL_ENERGY;
}

int EnergyFunctions::terminal(const std::size_t i, const std::size_t j, const std::string& sequence,
                              const std::vector<std::string_view>& mod_sequence,
                              const bool is_closing) {
    if (em_ == PKFreeEnergyModel::efn2 || em_ == PKFreeEnergyModel::enf2Simple) {
        // The end of the helix as the loop sees it: from inside the pair when it closes the loop.
        return is_closing ? RNAstructureFunctions::helix_end(j, i, sequence, mod_sequence, rsp_)
                          : RNAstructureFunctions::helix_end(i, j, sequence, mod_sequence, rsp_);
    }

    if (!mp_.empty()) {
        // The JSON files only have the pairs that involve a modified base.
        const int energy = ModBaseUtils::get_terminalAU_mod_energy(i, j, mod_sequence, mp_,
                                                                   vp_.md.temperature, is_closing);
        return energy != NULL_ENERGY ? energy
                                     : ViennaFunctions::terminal_energy(i, j, sequence, vp_);
    }

    if (em_ == PKFreeEnergyModel::ViennaRNA2) {
        return ViennaFunctions::terminal_energy(i, j, sequence, vp_);
    }

    return NULL_ENERGY;
}

int EnergyFunctions::hairpin(const std::size_t i, const std::size_t j, const std::string& sequence,
                             const std::vector<std::string_view>& mod_sequence, bool& is_inf) {
    if (em_ == PKFreeEnergyModel::efn2 || em_ == PKFreeEnergyModel::enf2Simple) {
        return RNAstructureFunctions::hairpin_energy(i, j, sequence, mod_sequence, rsp_, is_inf);
    }

    if (!mp_.empty()) {
        return ModHairpin::find_mod_hairpin_energy(i, j, sequence, mod_sequence, vp_, mp_, is_inf);
    }

    if (em_ == PKFreeEnergyModel::ViennaRNA2) {
        return ViennaFunctions::hairpin_energy(i, j, sequence, is_inf, vp_);
    }

    return NULL_ENERGY;
}

int EnergyFunctions::internal(const std::size_t i, const std::size_t j, const std::size_t ci,
                              const std::size_t cj, const std::string& sequence,
                              const std::vector<std::string_view>& mod_sequence) {
    if (em_ == PKFreeEnergyModel::efn2 || em_ == PKFreeEnergyModel::enf2Simple) {
        return RNAstructureFunctions::internal_loop_energy(i, j, ci, cj, sequence, mod_sequence,
                                                           rsp_);
    }

    if (!mp_.empty()) {
        return ModInternal::find_mod_internal_energy(i, j, ci, cj, mod_sequence, sequence, vp_,
                                                     mp_);
    }

    if (em_ == PKFreeEnergyModel::ViennaRNA2) {
        return ViennaFunctions::internal_loop_energy(i, j, ci, cj, sequence, vp_);
    }

    return NULL_ENERGY;
}

int EnergyFunctions::multibranch(const LoopNode& node, const ProcessedRNAEntry& pRNA) {
    if (em_ == PKFreeEnergyModel::efn2 || em_ == PKFreeEnergyModel::enf2Simple) {
        bool simple_multibranch_rule = em_ == PKFreeEnergyModel::enf2Simple;
        return RNAstructureFunctions::multibranch_energy(node, pRNA, rsp_,
                                                         simple_multibranch_rule
                                                             ? rnastructure::MultibranchRule::Simple
                                                             : rnastructure::MultibranchRule::Efn2);
    }

    if (!mp_.empty()) {
        return ModMultiloop::find_mod_multiloop_energy(node, pRNA, vp_, mp_);
    }

    if (em_ == PKFreeEnergyModel::ViennaRNA2) {
        return ViennaFunctions::multibranch_energy(node, pRNA, vp_);
    }

    return NULL_ENERGY;
}

int EnergyFunctions::external(const std::vector<LoopNode*>& children,
                              const ProcessedRNAEntry& pRNA) {
    if (em_ == PKFreeEnergyModel::efn2 || em_ == PKFreeEnergyModel::enf2Simple) {
        return RNAstructureFunctions::external_energy(children, pRNA, rsp_);
    }

    if (!mp_.empty()) {
        return ModExternal::find_mod_external_energy(children, pRNA, vp_, mp_);
    }

    if (em_ == PKFreeEnergyModel::ViennaRNA2) {
        return ViennaFunctions::external_energy(children, pRNA, vp_);
    }

    return NULL_ENERGY;
}
}  // namespace knotergy