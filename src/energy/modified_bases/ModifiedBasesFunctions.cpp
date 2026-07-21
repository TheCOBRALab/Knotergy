#include "ModifiedBasesFunctions.hpp"

namespace knotergy {

struct ModNode {
    unsigned int type;
    int          n5d;
    int          n3d;
    ModDiffs     diffs;
};

int ModifiedBasesFunctions::find_mod_external_energy(
    const std::vector<std::unique_ptr<LoopNode>>& children, const ProcessedRNAEntry& pRNA,
    const std::vector<std::string_view>& mod_sequence, vrna_md_param& vp,
    const all_mod_params& mp) {
    bool is_external = true;
    bool is_closing  = false;  // external loop energy correction does not apply to closing pair
    const std::string& sequence = pRNA.get_sequence();

    std::vector<DangleSet> all_dangle_sets;
    DangleSet              empty_set;
    int                    energy = 0;

    if (vp.md.dangles == 1 || vp.md.dangles == 3) {
        all_dangle_sets =
            Dangle1::populate_children_dangle_energies(children, pRNA, vp, is_external);
    } else {
        energy = ViennaFunctions::external_energy(children, pRNA, vp);
    }

    // Initialize energy with unmodified external energy
    // If dangles == 1, will be replaced later with corrected dangle energies so don't double count

    for (size_t idx = 0; idx < children.size(); ++idx) {
        const LoopNode& child = *children[idx];

        std::vector<std::string_view> unique_mod_bases =
            ModBaseUtils::unique_modified_bases_at_outer_edge(child.begin, child.end, mod_sequence);
        if (unique_mod_bases.empty()) continue;  // no modified bases in this child, skip to next

        // Get encoding of pair type and dangles for this child loop
        unsigned int type =
            ViennaUtils::get_pair_type(sequence[child.begin], sequence[child.end], vp.md);
        auto [n5d, n3d] = ViennaUtils::encode_outer_dangles(child.begin, child.end, pRNA, vp.md);

        // Get energy differences between modified and unmodified energies for this child loop
        ModDiffs diffs = get_mod_diffs(child, n5d, n3d, type, unique_mod_bases, mod_sequence, vp,
                                       mp, is_external, is_closing);

        if (vp.md.dangles == 1 || vp.md.dangles == 3) {
            DangleSet& current_set = all_dangle_sets[idx];
            modify_dangle_set(current_set, diffs);
        } else {
            energy += update_energy(diffs, n5d, n3d, type, vp);
        }
    }

    if (vp.md.dangles == 1 || vp.md.dangles == 3) {
        energy = Dangle1::get_external_dangle_1(children, all_dangle_sets);
    }

    return energy;
}

int ModifiedBasesFunctions::update_energy(const ModDiffs& diffs, int n5d, int n3d,
                                          unsigned int type, vrna_md_param& vp) {
    int total_energy_diff = 0;
    if (vp.md.dangles == 2) {
        if (n5d >= 0 && n3d >= 0)
            total_energy_diff += diffs.mismatch;
        else if (n5d >= 0)
            total_energy_diff += diffs.n5d;
        else if (n3d >= 0)
            total_energy_diff += diffs.n3d;
        if (type > 2) total_energy_diff += diffs.terminalAU;
    } else if (vp.md.dangles == 0) {
        if (type > 2) total_energy_diff += diffs.terminalAU;
    }
    return total_energy_diff;
}

// Gets the modified energy for a multiloop
int ModifiedBasesFunctions::find_mod_multiloop_energy(
    const LoopNode& node, const ProcessedRNAEntry& pRNA,
    const std::vector<std::string_view>& mod_sequence, vrna_md_param& vp,
    const all_mod_params& mp) {
    bool is_external    = false;
    bool is_closing     = true;   // multiloop energy correction only applies to closing pair
    bool is_not_closing = false;  // used for child loops

    // Get dangle energies for all children and closing pair if dangles == 1
    std::vector<DangleSet> children_dangle_sets;
    DangleSet              closing_set;
    DangleSet              empty_set;  // Used when dangles != 1 to avoid checking condition in loop
    // Dangle3
    std::vector<MultiloopStem> multiloop_stems;

    if (vp.md.dangles == 1) {
        children_dangle_sets =
            Dangle1::populate_children_dangle_energies(node.children, pRNA, vp, is_external);
        closing_set = Dangle1::get_ml_closing_dangle_energy(node, pRNA, vp);
    } else if (vp.md.dangles == 3) {
        multiloop_stems = CoaxialStacking::populate_multiloop_stems(node, pRNA, vp, mp);
    }

    int energy = 0;
    // Start with multiloop penalties
    if (vp.md.dangles != 1 && vp.md.dangles != 3) {
        energy = ViennaFunctions::multibranch_energy(node, pRNA, vp);
    } else {
        // Multiloop initialization penalties
        int ml_init = vp.p->MLclosing + node.exclusive_unpaired_bases_count * vp.p->MLbase;
        energy      = ml_init;
    }

    size_t                        i        = node.begin;
    size_t                        j        = node.end;
    const std::string&            sequence = pRNA.get_sequence();
    std::vector<std::string_view> unique_mod_bases =
        ModBaseUtils::unique_modified_bases_at_inner_edge(i, j, mod_sequence);

    if (!unique_mod_bases.empty()) {
        unsigned int pair_type = ViennaUtils::reverse_pair_type(sequence[i], sequence[j], vp.md);
        auto [n5d, n3d]        = ViennaUtils::encode_inner_dangles(i, j, pRNA, vp.md);
        ModDiffs diffs = get_mod_diffs(node, n3d, n5d, pair_type, unique_mod_bases, mod_sequence,
                                       vp, mp, is_external, is_closing);

        if (vp.md.dangles == 1) {
            modify_dangle_set(closing_set, diffs);
        } else if (vp.md.dangles == 3) {
            // closing stem is last in the vector
            MultiloopStem& closing_stem = multiloop_stems.back();
            closing_stem.dangle5 += diffs.n5d;
            closing_stem.dangle3 += diffs.n3d;
        } else {
            energy += update_energy(diffs, n5d, n3d, pair_type, vp);
        }
    }

    // Add energy corrections for modified bases in child loops
    for (size_t idx = 0; idx < node.children.size(); ++idx) {
        size_t                        ci = node.children[idx]->begin;
        size_t                        cj = node.children[idx]->end;
        std::vector<std::string_view> child_unique_mod_bases =
            ModBaseUtils::unique_modified_bases_at_outer_edge(ci, cj, mod_sequence);
        if (child_unique_mod_bases.empty())
            continue;  // no modified bases in this child, skip to next

        unsigned int c_pair_type    = ViennaUtils::get_pair_type(sequence[ci], sequence[cj], vp.md);
        auto [child_n5d, child_n3d] = ViennaUtils::encode_outer_dangles(ci, cj, pRNA, vp.md);
        ModDiffs child_diffs = get_mod_diffs(*node.children[idx], child_n5d, child_n3d, c_pair_type,
                                             child_unique_mod_bases, mod_sequence, vp, mp,
                                             is_external, is_not_closing);
        if (vp.md.dangles == 1) {
            DangleSet& current_set = (vp.md.dangles == 1) ? children_dangle_sets[idx] : closing_set;
            modify_dangle_set(current_set, child_diffs);
        } else if (vp.md.dangles == 3) {
            // Update the dangle energies for this child stem in the multiloop
            MultiloopStem& stem = multiloop_stems[idx];
            stem.dangle5 += child_diffs.n5d;
            stem.dangle3 += child_diffs.n3d;
        } else {
            energy += update_energy(child_diffs, child_n5d, child_n3d, c_pair_type, vp);
        }
    }

    // Compute the energy of dangle 1 of multiloop
    if (vp.md.dangles == 1) {
        energy += Dangle1::get_multibranch_dangle_1(node, children_dangle_sets, closing_set);
    }
    if (vp.md.dangles == 3) {
        energy += CoaxialStacking::get_multibranch_dangle_3(node, multiloop_stems, pRNA, vp, mp);
    }

    return energy;
}

void ModifiedBasesFunctions::modify_dangle_set(DangleSet& dangle_set, ModDiffs diffs) {
    dangle_set.both_dangle += diffs.mismatch;
    dangle_set.left_dangle += diffs.n5d;
    dangle_set.right_dangle += diffs.n3d;
    dangle_set += diffs.terminalAU;  // adds terminalAU diff to all configurations
}

ModDiffs ModifiedBasesFunctions::get_mod_diffs(
    const LoopNode& node, int n5d, int n3d, unsigned int type,
    const std::vector<std::string_view>& unique_mod_bases,
    const std::vector<std::string_view>& mod_sequence, vrna_md_param& vp, const all_mod_params& mp,
    bool is_external, bool is_closing) {
    if (is_external && is_closing) {
        THROW_ERROR("An external loop cannot be a closing pair, check loop tree construction");
    }

    vrna_param_t* P        = vp.p;
    int           mismatch = 0;
    int           dangle5  = 0;
    int           dangle3  = 0;
    int           terminal = 0;
    std::string   mismatch_key, dangle5_key, dangle3_key, terminal_key;
    // Unmodified energies for exterior stem (mismatch, dangle5, dangle3, terminalAU)
    if (n5d >= 0 && n3d >= 0) {
        mismatch = is_external ? P->mismatchExt[type][n5d][n3d] : P->mismatchM[type][n5d][n3d];
        mismatch_key =
            is_closing ? ModBaseUtils::join_string_views(
                             {node.end, node.end - 1, node.begin, node.begin + 1}, mod_sequence)
                       : ModBaseUtils::join_string_views(
                             {node.begin, node.begin - 1, node.end, node.end + 1}, mod_sequence);
    }

    if (n5d >= 0) {
        dangle5     = P->dangle5[type][n5d];
        dangle5_key = is_closing ? ModBaseUtils::join_string_views(
                                       {node.begin, node.end, node.end - 1}, mod_sequence)
                                 : ModBaseUtils::join_string_views(
                                       {node.begin, node.end, node.begin - 1}, mod_sequence);
    }

    if (n3d >= 0) {
        dangle3     = P->dangle3[type][n3d];
        dangle3_key = is_closing ? ModBaseUtils::join_string_views(
                                       {node.begin, node.end, node.begin + 1}, mod_sequence)
                                 : ModBaseUtils::join_string_views(
                                       {node.begin, node.end, node.end + 1}, mod_sequence);
    }

    if (type > 2) {
        terminal     = P->TerminalAU;
        terminal_key = is_closing
                           ? ModBaseUtils::join_string_views({node.end, node.begin}, mod_sequence)
                           : ModBaseUtils::join_string_views({node.begin, node.end}, mod_sequence);
    }

    // Get modified energies
    int modMismatch = ModBaseUtils::get_mod_energy(mismatch_key, unique_mod_bases, mp, mismatch,
                                                   ModLookup::Mismatch);
    int modDangle5  = ModBaseUtils::get_mod_energy(dangle5_key, unique_mod_bases, mp, dangle5,
                                                   ModLookup::Dangle5);
    int modDangle3  = ModBaseUtils::get_mod_energy(dangle3_key, unique_mod_bases, mp, dangle3,
                                                   ModLookup::Dangle3);
    int modTerminal = ModBaseUtils::get_mod_energy(terminal_key, unique_mod_bases, mp, terminal,
                                                   ModLookup::TerminalAU);

    // Calculate differences
    if (vp.md.dangles == 0) {
        return ModDiffs(modTerminal - terminal, 0, 0, 0);
    }

    int diffMismatch = modMismatch - mismatch;
    int diff5        = modDangle5 - dangle5;
    int diff3        = modDangle3 - dangle3;
    int diffTerminal = modTerminal - terminal;

    return ModDiffs(diffTerminal, diffMismatch, diff5, diff3);
}

}  // namespace knotergy