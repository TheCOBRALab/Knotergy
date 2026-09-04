#include "ModExternal.hpp"

#include "energy/vienna/ViennaFunctions.hpp"

namespace knotergy {
namespace {
int is_closing = false;  // external loop energy correction does not apply to closing pair
int is_external = true;  // external loop energy correction applies to children of external loop
}  // namespace

int ModExternal::find_mod_external_energy(const std::vector<LoopNode*>& children,
                                          const ProcessedRNAEntry& pRNA, vrna_md_param& vp,
                                          const all_mod_params& mp) {
    if (vp.md.dangles == 1 || vp.md.dangles == 3) {
        return external_dangle_1_energy(children, pRNA, vp, mp);
    } else if (vp.md.dangles == 0 || vp.md.dangles == 2) {
        return external_dangle_0_2_energy(children, pRNA, vp, mp);
    } else {
        THROW_ERROR("Invalid dangle setting: " + std::to_string(vp.md.dangles));
    }
}

ModDiffs ModExternal::get_external_child_diffs(const LoopNode& child, const ProcessedRNAEntry& pRNA,
                                               vrna_md_param& vp, const all_mod_params& mp) {
    const std::string& sequence = pRNA.get_sequence();
    const std::vector<std::string_view>& mod_sequence = pRNA.get_modified_sequence();

    std::vector<std::string_view> unique_mod_bases =
        ModBaseUtils::unique_modified_bases_at_outer_edge(child.begin, child.end, mod_sequence);
    if (unique_mod_bases.empty()) {
        return ModDiffs{};  // no modified bases in this child, return zero diffs
    }
    unsigned int type =
        ViennaUtils::get_pair_type(sequence[child.begin], sequence[child.end], vp.md);
    auto [n5d, n3d] = ViennaUtils::encode_outer_dangles(child.begin, child.end, pRNA, vp.md);

    ModDiffs diffs = ModBaseUtils::get_mod_diffs(child, n5d, n3d, type, unique_mod_bases,
                                                 mod_sequence, vp, mp, is_external, is_closing);
    return diffs;
}

int ModExternal::external_dangle_0_2_energy(const std::vector<LoopNode*>& children,
                                            const ProcessedRNAEntry& pRNA, vrna_md_param& vp,
                                            const all_mod_params& mp) {
    if (vp.md.dangles != 0 && vp.md.dangles != 2) {
        THROW_ERROR("external_dangle_0_2_energy should only be called when dangles == 0 or 2");
    }
    const std::string& sequence = pRNA.get_sequence();
    int energy = 0;

    for (const LoopNode* child : children) {
        std::size_t i = child->begin;
        std::size_t j = child->end;
        unsigned int type = ViennaUtils::get_pair_type(sequence[i], sequence[j], vp.md);
        auto [n5d, n3d] = ViennaUtils::encode_outer_dangles(i, j, pRNA, vp.md);
        if (vp.md.dangles == 0) {
            n5d = -1;
            n3d = -1;
        }
        energy += ModExternal::mod_exterior_stem(*child, n5d, n3d, type,
                                                 pRNA.get_modified_sequence(), vp, mp);
    }
    return energy;
}

int ModExternal::external_dangle_1_energy(const std::vector<LoopNode*>& children,
                                          const ProcessedRNAEntry& pRNA, vrna_md_param& vp,
                                          const all_mod_params& mp) {
    if (vp.md.dangles != 1 && vp.md.dangles != 3) {
        THROW_ERROR("external_dangle_1_energy should only be called when dangles == 1 or 3");
    }

    std::vector<DangleSet> all_dangle_sets =
        Dangle1::populate_children_dangle_energies(children, pRNA, vp, is_external);
    int energy = 0;

    for (std::size_t idx = 0; idx < children.size(); ++idx) {
        const LoopNode& child = *children[idx];
        ModDiffs diffs = ModExternal::get_external_child_diffs(child, pRNA, vp, mp);
        DangleSet& current_set = all_dangle_sets[idx];
        ModBaseUtils::modify_dangle_set(current_set, diffs);
    }

    energy = Dangle1::get_external_dangle_1(children, all_dangle_sets);

    return energy;
}
int ModExternal::mod_exterior_stem(const LoopNode& node, int n5d, int n3d, unsigned int type,
                                   const std::vector<std::string_view>& mod_sequence,
                                   vrna_md_param& vp, const all_mod_params& mp) {
    int energy = 0;
    if (n5d >= 0 && n3d >= 0) {
        int mod_mm = ModBaseUtils::get_mismatch_mod_energy(node.begin, node.end, mod_sequence, mp);
        energy += mod_mm != NULL_ENERGY ? mod_mm : vp.p->mismatchExt[type][n5d][n3d];
    } else if (n5d >= 0) {
        int mod_d5 = ModBaseUtils::get_dangle5_mod_energy(node.begin, node.end, mod_sequence, mp);
        energy += mod_d5 != NULL_ENERGY ? mod_d5 : vp.p->dangle5[type][n5d];
    } else if (n3d >= 0) {
        int mod_d3 = ModBaseUtils::get_dangle3_mod_energy(node.begin, node.end, mod_sequence, mp);
        energy += mod_d3 != NULL_ENERGY ? mod_d3 : vp.p->dangle3[type][n3d];
    }
    if (type > 2) {
        int mod_terminal =
            ModBaseUtils::get_terminalAU_mod_energy(node.begin, node.end, mod_sequence, mp);
        energy += mod_terminal != NULL_ENERGY ? mod_terminal : vp.p->TerminalAU;
    }

    return energy;
}
}  // namespace knotergy
