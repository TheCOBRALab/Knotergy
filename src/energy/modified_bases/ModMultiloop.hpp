#pragma once
#include "energy/dangles/Dangle1.hpp"
#include "energy/modified_bases/ModBaseUtils.hpp"
#include "energy/vienna/ViennaFunctions.hpp"

namespace knotergy {
class ModMultiloop {
   public:
    [[nodiscard]] static int find_mod_multiloop_energy(const LoopNode& node,
                                                       const ProcessedRNAEntry& pRNA,
                                                       vrna_md_param& vp, const all_mod_params& mp);

    [[nodiscard]] static int multiloop_dangle_0_2_energy(const LoopNode& node,
                                                         const ProcessedRNAEntry& pRNA,
                                                         vrna_md_param& vp,
                                                         const all_mod_params& mp);

    [[nodiscard]] static int multiloop_dangle_1_energy(const LoopNode& node,
                                                       const ProcessedRNAEntry& pRNA,
                                                       vrna_md_param& vp, const all_mod_params& mp);

    [[nodiscard]] static int multiloop_dangle_3_energy(const LoopNode& node,
                                                       const ProcessedRNAEntry& pRNA,
                                                       vrna_md_param& vp, const all_mod_params& mp);

    [[nodiscard]] static ModDiffs get_multiloop_diffs(const LoopNode& node,
                                                      const ProcessedRNAEntry& pRNA,
                                                      bool is_closing, vrna_md_param& vp,
                                                      const all_mod_params& mp);
};
}  // namespace knotergy