#include "ModInternal.hpp"

#include "energy/modified_bases/ModStack.hpp"
#include "utils/common.hpp"

// Modified-base corrections of internal loops, as in ViennaRNA (sc_PAIR_IL in
// constraints/sc_cb_mod.c):
//   - a stack or a bulge of one nucleotide gets the stacking correction;
//   - any other loop except a 1x1 internal loop gets the helix-end (AU/GU penalty) corrections of
//     both pairs;
//   - a loop with unpaired nucleotides on both sides gets the terminal-mismatch corrections of
//     both pairs.

namespace knotergy {

int ModInternal::find_mod_internal_energy(std::size_t i, std::size_t j, std::size_t ci,
                                          std::size_t cj,
                                          const std::vector<std::string_view>& mod_sequence,
                                          const std::string& sequence, vrna_md_param& vp,
                                          const all_mod_params& mp) {
    int unmod_energy = ViennaFunctions::internal_loop_energy(i, j, ci, cj, sequence, vp);
    // Both pairs and their terminal mismatches inside the loop
    std::vector<std::string_view> closing_unique_mod_bases =
        ModBaseUtils::unique_modified_bases_at_inner_edge(i, j, mod_sequence);
    std::vector<std::string_view> nested_unique_mod_bases =
        ModBaseUtils::unique_modified_bases_at_outer_edge(ci, cj, mod_sequence);
    if (closing_unique_mod_bases.empty() && nested_unique_mod_bases.empty()) {
        return unmod_energy;
    }
    double temperature = vp.md.temperature;
    int n1 = static_cast<int>(ci - i - 1);  // unpaired bases 5' side
    int n2 = static_cast<int>(j - cj - 1);  // unpaired bases 3' side
    int nl = std::max(n1, n2);
    int ns = std::min(n1, n2);
    if (nl == 0) {
        // This is a stack, so we can use the modified stack energy function
        return ModStack::find_mod_stack_energy(i, j, ci, cj, sequence, mod_sequence, vp, mp);
    }

    int e = unmod_energy;

    if (nl == 1 && ns == 0) {
        // This is a bulge of size 1, so we can use the modified stack energy function
        int unmod_stack_energy = ViennaFunctions::stack_energy(i, j, ci, cj, sequence, vp);
        int mod_stack_energy =
            ModStack::find_mod_stack_energy(i, j, ci, cj, sequence, mod_sequence, vp, mp);
        return e + mod_stack_energy - unmod_stack_energy;
    }

    // Helix ends of both pairs, except in 1x1 internal loops
    if (nl >= 2) {
        unsigned int type1 = ViennaUtils::get_pair_type(sequence[i], sequence[j], vp.md);
        unsigned int type2 = ViennaUtils::reverse_pair_type(sequence[ci], sequence[cj], vp.md);

        // Type 1
        int unmod_AU_penalty1 = type1 > 2 ? vp.p->TerminalAU : 0;
        std::string AU_key1 = ModBaseUtils::join_string_views({i, j}, mod_sequence);
        int mod_AU_penalty1 = ModBaseUtils::get_mod_energy(AU_key1, closing_unique_mod_bases, mp,
                                                           ModLookup::TerminalAU, temperature);

        int mod_AU1 = mod_AU_penalty1 != NULL_ENERGY ? mod_AU_penalty1 : unmod_AU_penalty1;
        e += mod_AU1 - unmod_AU_penalty1;

        // Type 2
        int unmod_AU_penalty2 = type2 > 2 ? vp.p->TerminalAU : 0;
        std::string AU_key2 = ModBaseUtils::join_string_views({ci, cj}, mod_sequence);
        int mod_AU_penalty2 = ModBaseUtils::get_mod_energy(AU_key2, nested_unique_mod_bases, mp,
                                                           ModLookup::TerminalAU, temperature);
        int mod_AU2 = mod_AU_penalty2 != NULL_ENERGY ? mod_AU_penalty2 : unmod_AU_penalty2;
        e += mod_AU2 - unmod_AU_penalty2;
    }

    // Terminal mismatches of both pairs, if there are unpaired bases on both sides
    if (ns >= 1) {
        e += ModBaseUtils::loop_mismatch_diff(i, j, sequence, mod_sequence, vp, mp, true);
        e += ModBaseUtils::loop_mismatch_diff(ci, cj, sequence, mod_sequence, vp, mp, false);
    }

    return e;
}

int ModInternal::find_mod_internal_energy(const PKBasePair& bp, const PKBasePair& next_bp,
                                          const ProcessedRNAEntry& pRNA, vrna_md_param& vp,
                                          const all_mod_params& mp) {
    return find_mod_internal_energy(bp.i, bp.j, next_bp.i, next_bp.j, pRNA.get_modified_sequence(),
                                    pRNA.get_sequence(), vp, mp);
}

int ModInternal::find_mod_internal_energy(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                          vrna_md_param& vp, const all_mod_params& mp) {
    return find_mod_internal_energy(node.begin, node.end, node.children[0]->begin,
                                    node.children[0]->end, pRNA.get_modified_sequence(),
                                    pRNA.get_sequence(), vp, mp);
}
}  // namespace knotergy
