#pragma once
#include "energy/modified_bases/ModBaseUtils.hpp"

namespace knotergy {
class ModMultiloop {
   public:
    [[nodiscard]] static int find_mod_external_energy(
        const std::vector<std::unique_ptr<LoopNode>>& children, const ProcessedRNAEntry& pRNA,
        const std::vector<std::string_view>& mod_sequence, vrna_md_param& vp,
        const all_mod_params& mp);

    [[nodiscard]] static int find_mod_multiloop_energy(
        const LoopNode& node, const ProcessedRNAEntry& pRNA,
        const std::vector<std::string_view>& mod_sequence, vrna_md_param& vp,
        const all_mod_params& mp);
};
}  // namespace knotergy