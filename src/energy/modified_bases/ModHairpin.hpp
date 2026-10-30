#pragma once

#include "energy/modified_bases/ModBaseUtils.hpp"
#include "energy/vienna/ViennaFunctions.hpp"
#include "loop_tree/LoopNode.hpp"
#include "utils/common.hpp"

namespace knotergy {

class ModHairpin {
   public:
    // Corrects the AU penalty and terminal mismatch of hairpin
    static int find_mod_hairpin_energy(std::size_t i, std::size_t j, const std::string& sequence,
                                       const std::vector<std::string_view>& mod_sequence,
                                       vrna_md_param& vp, const all_mod_params& mp, bool& is_inf) {
        // Get the unmodified hairpin energy
        int unmod_energy = ViennaFunctions::hairpin_energy(i, j, sequence, is_inf, vp);

        // If the hairpin energy is infinite, return it directly
        if (is_inf) {
            return unmod_energy;
        }

        // No modified base in the closing pair or its terminal mismatch
        if (ModBaseUtils::unique_modified_bases_at_inner_edge(i, j, mod_sequence).empty()) {
            return unmod_energy;
        }

        int e = unmod_energy;

        // Helix end of the closing pair
        std::vector<std::string_view> unique_mod_bases =
            ModBaseUtils::unique_modified_bases_at_indices({i, j}, mod_sequence);
        if (!unique_mod_bases.empty()) {
            unsigned int type = ViennaUtils::get_pair_type(sequence[i], sequence[j], vp.md);
            int unmod_AU_penalty =
                type > 2 ? vp.p->TerminalAU : 0;  // Only apply AU penalty for non-GC pairs

            double temperature = vp.md.temperature;
            std::string key = ModBaseUtils::join_string_views({i, j}, mod_sequence);
            int modAU = ModBaseUtils::get_mod_energy(key, unique_mod_bases, mp,
                                                     ModLookup::TerminalAU, temperature);
            modAU = modAU != NULL_ENERGY ? modAU : unmod_AU_penalty;
            e += modAU - unmod_AU_penalty;
        }

        // Terminal mismatch inside the closing pair
        e += ModBaseUtils::loop_mismatch_diff(i, j, sequence, mod_sequence, vp, mp, true);

        return e;
    }

    static int find_mod_hairpin_energy(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                       const std::vector<std::string_view>& mod_sequence,
                                       vrna_md_param& vp, const all_mod_params& mp, bool& is_inf) {
        return find_mod_hairpin_energy(node.begin, node.end, pRNA.get_sequence(), mod_sequence, vp,
                                       mp, is_inf);
    }
};
}  // namespace knotergy