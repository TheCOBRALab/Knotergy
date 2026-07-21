#include "ModMultiloop.hpp"

namespace knotergy {

namespace {
bool is_external = true;  // external loop energy correction applies to children of external loop
}  // namespace

ModDiffs ModMultiloop::get_multiloop_diffs(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                           bool is_closing, vrna_md_param& vp,
                                           const all_mod_params& mp) {
    const std::string& sequence = pRNA.get_sequence();
    const std::vector<std::string_view>& mod_sequence = pRNA.get_modified_sequence();
    std::vector<std::string_view> unique_mod_bases;

    if (is_closing) {
        unique_mod_bases =
            ModBaseUtils::unique_modified_bases_at_inner_edge(node.begin, node.end, mod_sequence);
    } else {
        unique_mod_bases =
            ModBaseUtils::unique_modified_bases_at_outer_edge(node.begin, node.end, mod_sequence);
    }

    if (unique_mod_bases.empty()) {
        return ModDiffs{};  // no modified bases in this child, return zero diffs
    }

    unsigned int type = ViennaUtils::get_pair_type(sequence[node.begin], sequence[node.end], vp.md);
    auto [n5d, n3d] = ViennaUtils::encode_outer_dangles(node.begin, node.end, pRNA, vp.md);

    ModDiffs diffs = ModBaseUtils::get_mod_diffs(node, n5d, n3d, type, unique_mod_bases,
                                                 mod_sequence, vp, mp, is_external, is_closing);
    return diffs;
}

}  // namespace knotergy
