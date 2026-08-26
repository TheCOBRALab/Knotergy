#pragma once

#include "energy/modified_bases/ModBaseUtils.hpp"
#include "energy/vienna/ViennaFunctions.hpp"
#include "loop_tree/LoopNode.hpp"
#include "utils/common.hpp"

namespace knotergy {

class ModHairpin {
   public:
    // Adds an AU penalty to the modified hairpin energy
    static int find_mod_hairpin_energy(size_t i, size_t j, const std::string& sequence,
                                       const std::vector<std::string_view>& mod_sequence,
                                       vrna_md_param& vp, const all_mod_params& mp, bool& is_inf) {
        // Get the unmodified hairpin energy
        int unmod_energy = ViennaFunctions::hairpin_energy(i, j, sequence, is_inf, vp);

        // If the hairpin energy is infinite, return it directly
        if (is_inf) {
            return unmod_energy;
        }

        // Create a key for modified base lookup
        std::vector<std::string_view> unique_mod_bases =
            ModBaseUtils::unique_modified_bases_at_indices({i, j}, mod_sequence);

        if (unique_mod_bases.empty()) {
            // No modified bases found, return the unmodified energy
            return unmod_energy;
        }

        // Create a key for modified base lookup
        std::string key = ModBaseUtils::join_string_views({i, j}, mod_sequence);
        int modAU = ModBaseUtils::get_mod_energy(key, unique_mod_bases, mp, ModLookup::TerminalAU);
        modAU = modAU != NULL_ENERGY ? modAU : vp.p->TerminalAU;

        // Look up the modified hairpin energy
        return unmod_energy + modAU - vp.p->TerminalAU;
    }

    static int find_mod_hairpin_energy(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                       const std::vector<std::string_view>& mod_sequence,
                                       vrna_md_param& vp, const all_mod_params& mp, bool& is_inf) {
        return find_mod_hairpin_energy(node.begin, node.end, pRNA.get_sequence(), mod_sequence, vp,
                                       mp, is_inf);
    }
};
}  // namespace knotergy