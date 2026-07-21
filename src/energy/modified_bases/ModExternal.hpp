#pragma once
#include "energy/dangles/Dangle1.hpp"
#include "energy/modified_bases/ModBaseUtils.hpp"

namespace knotergy {
class ModExternal {
   public:
    [[nodiscard]] static int find_mod_external_energy(
        const std::vector<std::unique_ptr<LoopNode>>& children, const ProcessedRNAEntry& pRNA,
        vrna_md_param& vp, const all_mod_params& mp);

    [[nodiscard]] static ModDiffs get_external_child_diffs(const LoopNode& child,
                                                           const ProcessedRNAEntry& pRNA,
                                                           vrna_md_param& vp,
                                                           const all_mod_params& mp);

    [[nodiscard]] static int external_dangle_0_2_energy(
        const std::vector<std::unique_ptr<LoopNode>>& children, const ProcessedRNAEntry& pRNA,
        vrna_md_param& vp, const all_mod_params& mp);

    [[nodiscard]] static int external_dangle_1_energy(
        const std::vector<std::unique_ptr<LoopNode>>& children, const ProcessedRNAEntry& pRNA,
        vrna_md_param& vp, const all_mod_params& mp);
};
}  // namespace knotergy