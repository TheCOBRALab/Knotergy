#include "ModMultiloop.hpp"

#include "energy/vienna/ViennaFunctions.hpp"
// Add the appropriate CoaxialStacking include here if ModMultiloop.hpp
// does not already provide it.

namespace knotergy {
namespace {

constexpr bool kIsExternal = false;
constexpr bool kIsClosing = true;
constexpr bool kIsNotClosing = false;

/**
 * Apply modified-base corrections when ViennaRNA uses dangles 0 or 2.
 *
 * For dangles == 0, only terminal-AU contributes.
 * For dangles == 2, mismatch/dangle corrections and terminal-AU contribute.
 */
void add_dangle_0_2_diffs(int& energy, const ModDiffs& diffs, int n5d, int n3d, unsigned int type,
                          int dangles) {
    if (dangles == 2) {
        if (n5d >= 0 && n3d >= 0) {
            energy += diffs.mismatch;
        } else if (n5d >= 0) {
            energy += diffs.n5d;
        } else if (n3d >= 0) {
            energy += diffs.n3d;
        }
    }

    if (type > 2) {
        energy += diffs.terminalAU;
    }
}

}  // namespace

int ModMultiloop::find_mod_multiloop_energy(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                            vrna_md_param& vp, const all_mod_params& mp) {
    switch (vp.md.dangles) {
        case 0:
        case 2:  return multiloop_dangle_0_2_energy(node, pRNA, vp, mp);
        case 1:  return multiloop_dangle_1_energy(node, pRNA, vp, mp);
        case 3:  return multiloop_dangle_3_energy(node, pRNA, vp, mp);
        default: THROW_ERROR("Invalid dangle setting: " + std::to_string(vp.md.dangles));
    }
}

ModDiffs ModMultiloop::get_multiloop_diffs(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                           bool is_closing, vrna_md_param& vp,
                                           const all_mod_params& mp) {
    const std::string& sequence = pRNA.get_sequence();
    const std::vector<std::string_view>& mod_sequence = pRNA.get_modified_sequence();

    std::vector<std::string_view> unique_mod_bases;

    unsigned int type;
    int n5d;
    int n3d;

    if (is_closing) {
        unique_mod_bases =
            ModBaseUtils::unique_modified_bases_at_inner_edge(node.begin, node.end, mod_sequence);

        if (unique_mod_bases.empty()) {
            return ModDiffs{};
        }

        type = ViennaUtils::reverse_pair_type(sequence[node.begin], sequence[node.end], vp.md);

        std::tie(n5d, n3d) = ViennaUtils::encode_inner_dangles(node.begin, node.end, pRNA, vp.md);

        /*
         * Important:
         *
         * The closing pair is reversed relative to an ordinary outer stem.
         * Therefore n3d, n5d is passed into get_mod_diffs().
         */
        return ModBaseUtils::get_mod_diffs(node, n3d, n5d, type, unique_mod_bases, mod_sequence, vp,
                                           mp, kIsExternal, kIsClosing);
    }

    unique_mod_bases =
        ModBaseUtils::unique_modified_bases_at_outer_edge(node.begin, node.end, mod_sequence);

    if (unique_mod_bases.empty()) {
        return ModDiffs{};
    }

    type = ViennaUtils::get_pair_type(sequence[node.begin], sequence[node.end], vp.md);

    std::tie(n5d, n3d) = ViennaUtils::encode_outer_dangles(node.begin, node.end, pRNA, vp.md);

    return ModBaseUtils::get_mod_diffs(node, n5d, n3d, type, unique_mod_bases, mod_sequence, vp, mp,
                                       kIsExternal, kIsNotClosing);
}

int ModMultiloop::multiloop_dangle_0_2_energy(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                              vrna_md_param& vp, const all_mod_params& mp) {
    if (vp.md.dangles != 0 && vp.md.dangles != 2) {
        THROW_ERROR(
            "multiloop_dangle_0_2_energy should only be called when "
            "dangles == 0 or 2");
    }

    int energy = ViennaFunctions::multibranch_energy(node, pRNA, vp);
    const std::string& sequence = pRNA.get_sequence();

    // Correct the closing stem.
    {
        const unsigned int type =
            ViennaUtils::reverse_pair_type(sequence[node.begin], sequence[node.end], vp.md);

        const auto [n5d, n3d] =
            ViennaUtils::encode_inner_dangles(node.begin, node.end, pRNA, vp.md);

        const ModDiffs diffs = get_multiloop_diffs(node, pRNA, kIsClosing, vp, mp);

        add_dangle_0_2_diffs(energy, diffs, n5d, n3d, type, vp.md.dangles);
    }

    // Correct every child stem.
    for (const std::unique_ptr<LoopNode>& child_ptr : node.children) {
        const LoopNode& child = *child_ptr;

        const unsigned int type =
            ViennaUtils::get_pair_type(sequence[child.begin], sequence[child.end], vp.md);

        const auto [n5d, n3d] =
            ViennaUtils::encode_outer_dangles(child.begin, child.end, pRNA, vp.md);

        const ModDiffs diffs = get_multiloop_diffs(child, pRNA, kIsNotClosing, vp, mp);

        add_dangle_0_2_diffs(energy, diffs, n5d, n3d, type, vp.md.dangles);
    }

    return energy;
}

int ModMultiloop::multiloop_dangle_1_energy(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                            vrna_md_param& vp, const all_mod_params& mp) {
    if (vp.md.dangles != 1) {
        THROW_ERROR(
            "multiloop_dangle_1_energy should only be called when "
            "dangles == 1");
    }

    // Build the unmodified candidate dangle configurations.
    std::vector<DangleSet> children_dangle_sets =
        Dangle1::populate_children_dangle_energies(node.children, pRNA, vp, kIsExternal);
    DangleSet closing_set = Dangle1::get_ml_closing_dangle_energy(node, pRNA, vp);

    // Apply corrections to the closing stem's candidate configurations.
    {
        const ModDiffs closing_diffs = get_multiloop_diffs(node, pRNA, kIsClosing, vp, mp);
        ModBaseUtils::modify_dangle_set(closing_set, closing_diffs);
    }

    // Apply corrections to each child stem's candidate configurations.
    for (size_t idx = 0; idx < node.children.size(); ++idx) {
        const LoopNode& child = *node.children[idx];

        const ModDiffs child_diffs = get_multiloop_diffs(child, pRNA, kIsNotClosing, vp, mp);

        ModBaseUtils::modify_dangle_set(children_dangle_sets[idx], child_diffs);
    }

    /*
     * For dangles == 1, ViennaFunctions::multibranch_energy() is not used
     * because the dangle-selection routine calculates that portion.
     */
    int energy =
        vp.p->MLclosing + static_cast<int>(node.exclusive_unpaired_bases_count) * vp.p->MLbase;

    energy += Dangle1::get_multibranch_dangle_1(node, children_dangle_sets, closing_set);

    return energy;
}

int ModMultiloop::multiloop_dangle_3_energy(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                            vrna_md_param& vp, const all_mod_params& mp) {
    if (vp.md.dangles != 3) {
        THROW_ERROR(
            "multiloop_dangle_3_energy should only be called when "
            "dangles == 3");
    }

    std::vector<MultiloopStem> multiloop_stems =
        CoaxialStacking::populate_multiloop_stems(node, pRNA, vp, mp);

    /*
     * The closing stem is stored last.
     */
    {
        ModDiffs closing_diffs = get_multiloop_diffs(node, pRNA, kIsClosing, vp, mp);

        MultiloopStem& closing_stem = multiloop_stems.back();

        closing_stem.dangle5 += closing_diffs.n5d;
        closing_stem.dangle3 += closing_diffs.n3d;
    }

    /*
     * Child stems occupy the first node.children.size() entries.
     */
    for (size_t idx = 0; idx < node.children.size(); ++idx) {
        const LoopNode& child = *node.children[idx];

        ModDiffs child_diffs = get_multiloop_diffs(child, pRNA, kIsNotClosing, vp, mp);

        MultiloopStem& stem = multiloop_stems[idx];

        stem.dangle5 += child_diffs.n5d;
        stem.dangle3 += child_diffs.n3d;
    }

    int energy =
        vp.p->MLclosing + static_cast<int>(node.exclusive_unpaired_bases_count) * vp.p->MLbase;

    energy += CoaxialStacking::get_multibranch_dangle_3(node, multiloop_stems, pRNA, vp, mp);

    return energy;
}

}  // namespace knotergy