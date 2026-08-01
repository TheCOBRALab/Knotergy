#include "ModExternal.hpp"

#include "energy/vienna/ViennaFunctions.hpp"

namespace knotergy {
namespace {
int is_closing = false;  // external loop energy correction does not apply to closing pair
int is_external = true;  // external loop energy correction applies to children of external loop
}  // namespace

int ModExternal::find_mod_external_energy(const std::vector<std::unique_ptr<LoopNode>>& children,
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

int ModExternal::external_dangle_0_2_energy(const std::vector<std::unique_ptr<LoopNode>>& children,
                                            const ProcessedRNAEntry& pRNA, vrna_md_param& vp,
                                            const all_mod_params& mp) {
    if (vp.md.dangles != 0 && vp.md.dangles != 2) {
        THROW_ERROR("external_dangle_0_2_energy should only be called when dangles == 0 or 2");
    }

    int energy = ViennaFunctions::external_energy(children, pRNA, vp);

    for (const std::unique_ptr<LoopNode>& child : children) {
        unsigned int type = ViennaUtils::get_pair_type(pRNA.get_sequence()[child->begin],
                                                       pRNA.get_sequence()[child->end], vp.md);
        auto [n5d, n3d] = ViennaUtils::encode_outer_dangles(child->begin, child->end, pRNA, vp.md);
        ModDiffs diffs = ModExternal::get_external_child_diffs(*child, pRNA, vp, mp);

        if (vp.md.dangles == 2) {
            if (n5d >= 0 && n3d >= 0)
                energy += diffs.mismatch;
            else if (n5d >= 0)
                energy += diffs.n5d;
            else if (n3d >= 0)
                energy += diffs.n3d;
        }
        if (type > 2) energy += diffs.terminalAU;
    }
    return energy;
}

int ModExternal::external_dangle_1_energy(const std::vector<std::unique_ptr<LoopNode>>& children,
                                          const ProcessedRNAEntry& pRNA, vrna_md_param& vp,
                                          const all_mod_params& mp) {
    if (vp.md.dangles != 1 && vp.md.dangles != 3) {
        THROW_ERROR("external_dangle_1_energy should only be called when dangles == 1 or 3");
    }

    std::vector<DangleSet> all_dangle_sets =
        Dangle1::populate_children_dangle_energies(children, pRNA, vp, is_external);
    int energy = 0;

    for (size_t idx = 0; idx < children.size(); ++idx) {
        const LoopNode& child = *children[idx];
        ModDiffs diffs = ModExternal::get_external_child_diffs(child, pRNA, vp, mp);
        DangleSet& current_set = all_dangle_sets[idx];
        ModBaseUtils::modify_dangle_set(current_set, diffs);
    }

    energy = Dangle1::get_external_dangle_1(children, all_dangle_sets);

    return energy;
}
}  // namespace knotergy
