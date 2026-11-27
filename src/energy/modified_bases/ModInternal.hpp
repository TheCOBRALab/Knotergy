#pragma once

#include "energy/modified_bases/ModBaseUtils.hpp"
#include "energy/vienna/ViennaFunctions.hpp"

namespace knotergy {
class ModInternal {
   public:
    static int find_mod_internal_energy(std::size_t i, std::size_t j, std::size_t ci,
                                        std::size_t cj,
                                        const std::vector<std::string_view>& mod_sequence,
                                        const std::string& sequence, const vrna_md_param& vp,
                                        const all_mod_params& mp);
    static int find_mod_internal_energy(const PKBasePair& bp, const PKBasePair& next_bp,
                                        const ProcessedRNAEntry& processed_rna,
                                        const vrna_md_param& vp, const all_mod_params& mp);
    static int find_mod_internal_energy(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                        const vrna_md_param& vp, const all_mod_params& mp);
};

}  // namespace knotergy