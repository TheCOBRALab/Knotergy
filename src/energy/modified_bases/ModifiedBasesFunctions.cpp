#include "ModifiedBasesFunctions.hpp"

namespace knotergy {

// Returns a vector of unique modified bases found at the specified indices
std::vector<std::string_view> ModifiedBasesFunctions::unique_modified_bases_at_indices(
    std::vector<size_t> indices, const std::vector<std::string_view>& mod_sequence) {
    std::vector<std::string_view> modified;
    modified.reserve(indices.size());
    for (size_t idx : indices) {
        if (!RNAProcessor::is_unmod_base(mod_sequence[idx]) &&
            std::find(modified.begin(), modified.end(), mod_sequence[idx]) == modified.end()) {
            modified.push_back(mod_sequence[idx]);
        }
    }
    return modified;
}

std::vector<std::string_view> ModifiedBasesFunctions::unique_modified_bases_at_inner_edge(
    size_t i, size_t j, const std::vector<std::string_view>& mod_sequence) {
    return unique_modified_bases_at_indices({i, j, i + 1, j - 1}, mod_sequence);
}

std::vector<std::string_view> ModifiedBasesFunctions::unique_modified_bases_at_outer_edge(
    size_t i, size_t j, const std::vector<std::string_view>& mod_sequence) {
    std::vector<size_t> indices;
    indices.reserve(4);
    indices.push_back(i);
    indices.push_back(j);
    if (i > 0) indices.push_back(i - 1);
    if (j + 1 < mod_sequence.size()) indices.push_back(j + 1);
    return unique_modified_bases_at_indices(indices, mod_sequence);
}

// Gets the modified energy of a stack
int ModifiedBasesFunctions::find_mod_stack_energy(size_t i, size_t j, size_t ci, size_t cj,
                                                  const std::string& sequence,
                                                  const std::vector<std::string_view>& mod_sequence,
                                                  vrna_md_param& vp, const all_mod_params& mp) {
    // Get unmodified energy first to use as fallback if no modified nucleotides are found
    int unmod_energy = ViennaFunctions::stack_energy(i, j, ci, cj, sequence, vp);

    // Find all modified bases at the inner edge of the stack (i, j, i+1, j-1)
    std::vector<std::string_view> unique_mod_bases =
        unique_modified_bases_at_inner_edge(i, j, mod_sequence);
    if (unique_mod_bases.empty()) return unmod_energy;

    // Used to look up stacking energies in modified base parameters
    std::string l_key = join_string_views({i, ci, j, cj}, mod_sequence);

    // Get energy correction for modified bases (returns original energy if no modifications found)
    int e = get_mod_energy(l_key, unique_mod_bases, mp, unmod_energy, ModLookup::Stacking);

    // // uncomment if you want values to match RNAfold, but note this is likely a bug in RNAfold
    // // If no key found for the stack, try looking up the reverse stack (ci, i, cj, j)
    // if (e == unmod_energy) {
    //     std::string r_key = join_string_views({cj, j, ci, i}, mod_sequence);
    //     e = get_mod_energy(r_key, unique_mod_bases, mp, unmod_energy,
    //     ModLookup::Stacking);
    // }

    return e;
}

struct ModNode {
    unsigned int type;
    int n5d;
    int n3d;
    ModDiffs diffs;
};

int ModifiedBasesFunctions::find_mod_external_energy(
    const std::vector<std::unique_ptr<LoopNode>>& children, const ProcessedRNAEntry& pRNA,
    const std::vector<std::string_view>& mod_sequence, vrna_md_param& vp,
    const all_mod_params& mp) {
    bool is_external = true;
    bool is_closing = false;  // external loop energy correction does not apply to closing pair
    const std::string& sequence = pRNA.get_sequence();

    std::vector<DangleSet> all_dangle_sets;
    DangleSet empty_set;
    int energy = 0;

    if (vp.md.dangles == 1) {
        all_dangle_sets =
            ViennaDangles::populate_children_dangle_energies(children, pRNA, vp, is_external);
    } else {
        energy = ViennaFunctions::external_energy(children, pRNA, vp);
    }

    // Initialize energy with unmodified external energy
    // If dangles == 1, will be replaced later with corrected dangle energies so don't double count

    for (size_t idx = 0; idx < children.size(); ++idx) {
        const LoopNode& child = *children[idx];

        std::vector<std::string_view> unique_mod_bases =
            unique_modified_bases_at_outer_edge(child.begin, child.end, mod_sequence);
        if (unique_mod_bases.empty()) continue;  // no modified bases in this child, skip to next

        // Get encoding of pair type and dangles for this child loop
        unsigned int type =
            ViennaUtils::get_pair_type(sequence[child.begin], sequence[child.end], vp.md);
        auto [n5d, n3d] = ViennaUtils::encode_outer_dangles(child.begin, child.end, pRNA, vp.md);

        // Get energy differences between modified and unmodified energies for this child loop
        ModDiffs diffs = get_mod_diffs(child, n5d, n3d, type, unique_mod_bases, mod_sequence, vp,
                                       mp, is_external, is_closing);

        if (vp.md.dangles == 1) {
            DangleSet& current_set = all_dangle_sets[idx];
            modify_dangle_set(current_set, diffs);
        } else {
            energy += update_energy(diffs, n5d, n3d, type, vp);
        }
    }

    if (vp.md.dangles == 1) {
        energy = ViennaDangles::get_external_dangle_1(children, all_dangle_sets);
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
    bool is_external = false;
    bool is_closing = true;       // multiloop energy correction only applies to closing pair
    bool is_not_closing = false;  // used for child loops

    // Get dangle energies for all children and closing pair if dangles == 1
    std::vector<DangleSet> children_dangle_sets;
    DangleSet closing_set;
    DangleSet empty_set;  // Used when dangles != 1 to avoid checking condition in loop
    int energy = 0;
    size_t i = node.begin;
    size_t j = node.end;
    const std::string& sequence = pRNA.get_sequence();
    std::vector<std::string_view> unique_mod_bases =
        unique_modified_bases_at_inner_edge(i, j, mod_sequence);

    if (!unique_mod_bases.empty()) {
        unsigned int pair_type = ViennaUtils::reverse_pair_type(sequence[i], sequence[j], vp.md);
        auto [n5d, n3d] = ViennaUtils::encode_inner_dangles(i, j, pRNA, vp.md);
        ModDiffs diffs = get_mod_diffs(node, n3d, n5d, pair_type, unique_mod_bases, mod_sequence,
                                       vp, mp, is_external, is_closing);

        if (vp.md.dangles == 1) {
            children_dangle_sets = ViennaDangles::populate_children_dangle_energies(
                node.children, pRNA, vp, is_external);
            closing_set = ViennaDangles::get_ml_closing_dangle_energy(node, pRNA, vp);
            modify_dangle_set(closing_set, diffs);
        } else if (vp.md.dangles != 1) {
            energy = ViennaFunctions::multibranch_energy(node, pRNA, vp);
            energy += update_energy(diffs, n5d, n3d, pair_type, vp);
        }
    }

    // Add energy corrections for modified bases in child loops
    for (size_t idx = 0; idx < node.children.size(); ++idx) {
        size_t ci = node.children[idx]->begin;
        size_t cj = node.children[idx]->end;
        std::vector<std::string_view> child_unique_mod_bases =
            unique_modified_bases_at_outer_edge(ci, cj, mod_sequence);
        if (child_unique_mod_bases.empty())
            continue;  // no modified bases in this child, skip to next

        unsigned int c_pair_type = ViennaUtils::get_pair_type(sequence[ci], sequence[cj], vp.md);
        auto [child_n5d, child_n3d] = ViennaUtils::encode_outer_dangles(ci, cj, pRNA, vp.md);
        ModDiffs child_diffs = get_mod_diffs(*node.children[idx], child_n5d, child_n3d, c_pair_type,
                                             child_unique_mod_bases, mod_sequence, vp, mp,
                                             is_external, is_not_closing);
        if (vp.md.dangles == 1) {
            DangleSet& current_set = (vp.md.dangles == 1) ? children_dangle_sets[idx] : closing_set;
            modify_dangle_set(current_set, child_diffs);
        } else {
            energy += update_energy(child_diffs, child_n5d, child_n3d, c_pair_type, vp);
        }
    }

    // Compute the energy of dangle 1 of multiloop
    if (vp.md.dangles == 1) {
        energy = ViennaDangles::get_multibranch_dangle_1(node, children_dangle_sets, closing_set);
    }

    return energy;
}

// Uses lookup key to get the energy from the modified base parameter file
int ModifiedBasesFunctions::get_mod_energy(const std::string& key,
                                           const std::vector<std::string_view>& unique_mod_bases,
                                           const all_mod_params& mp, int unmod_energy,
                                           ModLookup lookup_type) {
    // If no modified bases are present, return the unmodified energy
    if (unique_mod_bases.empty()) {
        return static_cast<int>(unmod_energy);
    }

    // Get the pointer to the correct energy map based on lookup type
    for (const std::string_view& mod_base : unique_mod_bases) {
        const modified_base_param* param = mp.get_modified_base_param(std::string(mod_base));

        if (!param) {
            THROW_ERROR("Modified base '" + std::string(mod_base) +
                        "' found in sequence but no parameters provided for it.");
        }

        const std::map<std::string, float>* energy_lookup = nullptr;
        switch (lookup_type) {
            case ModLookup::Stacking:
                energy_lookup = &param->stacking_energies();
                break;
            case ModLookup::Terminal:
                energy_lookup = &param->terminal_energies();
                break;
            case ModLookup::Mismatch:
                energy_lookup = &param->mismatch_energies();
                break;
            case ModLookup::Dangle5:
                energy_lookup = &param->dangle5_energies();
                break;
            case ModLookup::Dangle3:
                energy_lookup = &param->dangle3_energies();
                break;
            default:
                THROW_ERROR("Invalid ModLookup type: " +
                            std::to_string(static_cast<int>(lookup_type)));
                break;
        }

        // If the map exists, look up the energy for this key and return if found
        if (energy_lookup) {
            auto it = energy_lookup->find(key);
            if (it != energy_lookup->end()) {
                int mod_energy = static_cast<int>(it->second * 100);
                // std::cout << "Modified base energy found for key: " << key << " -> " <<
                // mod_energy
                //           << " Diff: " << mod_energy - unmod_energy
                //           << std::endl;
                return mod_energy;
            }
        }
    }

    return static_cast<int>(unmod_energy);
}

// Joins string views at specified indices into a single string
// Move to general utilities
std::string ModifiedBasesFunctions::join_string_views(
    std::vector<size_t> indices, const std::vector<std::string_view>& mod_sequence) {
    std::string key;
    size_t total = 0;

    // Calculate total size needed
    for (size_t idx : indices) total += mod_sequence[idx].size();
    key.reserve(total);

    // Concatenate string views
    for (size_t idx : indices) {
        key.append(mod_sequence[idx]);
    }
    return key;
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

    vrna_param_t* P = vp.p;
    int mismatch = 0;
    int dangle5 = 0;
    int dangle3 = 0;
    int terminal = 0;
    std::string mismatch_key, dangle5_key, dangle3_key, terminal_key;
    // Unmodified energies for exterior stem (mismatch, dangle5, dangle3, terminalAU)
    if (n5d >= 0 && n3d >= 0) {
        mismatch = is_external ? P->mismatchExt[type][n5d][n3d] : P->mismatchM[type][n5d][n3d];
        mismatch_key = is_closing
                           ? join_string_views({node.end, node.end - 1, node.begin, node.begin + 1},
                                               mod_sequence)
                           : join_string_views({node.begin, node.begin - 1, node.end, node.end + 1},
                                               mod_sequence);
    }

    if (n5d >= 0) {
        dangle5 = P->dangle5[type][n5d];
        dangle5_key = is_closing
                          ? join_string_views({node.begin, node.end, node.end - 1}, mod_sequence)
                          : join_string_views({node.begin, node.end, node.begin - 1}, mod_sequence);
    }

    if (n3d >= 0) {
        dangle3 = P->dangle3[type][n3d];
        dangle3_key = is_closing
                          ? join_string_views({node.begin, node.end, node.begin + 1}, mod_sequence)
                          : join_string_views({node.begin, node.end, node.end + 1}, mod_sequence);
    }

    if (type > 2) {
        terminal = P->TerminalAU;
        terminal_key = is_closing ? join_string_views({node.end, node.begin}, mod_sequence)
                                  : join_string_views({node.begin, node.end}, mod_sequence);
    }

    // Get modified energies
    int modMismatch =
        get_mod_energy(mismatch_key, unique_mod_bases, mp, mismatch, ModLookup::Mismatch);
    int modDangle5 = get_mod_energy(dangle5_key, unique_mod_bases, mp, dangle5, ModLookup::Dangle5);
    int modDangle3 = get_mod_energy(dangle3_key, unique_mod_bases, mp, dangle3, ModLookup::Dangle3);
    int modTerminal =
        get_mod_energy(terminal_key, unique_mod_bases, mp, terminal, ModLookup::Terminal);

    // Calculate differences
    if (vp.md.dangles == 0) {
        return ModDiffs(modTerminal - terminal, 0, 0, 0);
    }

    int diffMismatch = modMismatch - mismatch;
    int diff5 = modDangle5 - dangle5;
    int diff3 = modDangle3 - dangle3;
    int diffTerminal = modTerminal - terminal;

    return ModDiffs(diffTerminal, diffMismatch, diff5, diff3);
}

}  // namespace knotergy