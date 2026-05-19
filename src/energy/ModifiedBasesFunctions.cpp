#include "ModifiedBasesFunctions.hpp"

namespace knotergy {

// Gets the modified energy of a stack
int ModifiedBasesFunctions::find_mod_stack_energy(size_t i, size_t j, size_t ci, size_t cj,
                                                  const std::string& sequence,
                                                  const std::vector<std::string_view>& mod_sequence,
                                                  vrna_md_param& vp, const all_mod_params& mp) {
    // Get Vienna stacking energy for unmodified bases
    int unmod_energy = ViennaFunctions::stack_energy(i, j, ci, cj, sequence, vp);

    // Get the grapheme views of modified bases in the relevant positions (if any)
    std::vector<std::string_view> unique_mod_bases =
        unique_modified_bases_at_indices({i, j, ci, cj}, mod_sequence);
    if (unique_mod_bases.empty()) return unmod_energy;

    // Used to look up stacking energies in modified base parameters
    std::string l_key = join_string_views({i, ci, j, cj}, mod_sequence);
    std::string r_key = join_string_views({cj, j, ci, i}, mod_sequence);

    // Get energy correction for modified bases (returns original energy if no modifications found)
    // Try the Vienna-defined ordering first, then the swapped orientation
    int e = get_mod_energy(l_key, unique_mod_bases, mp, unmod_energy, ModLookup::Stacking);

    // // uncomment if you want values to match RNAfold, but note this is likely a bug in RNAfold
    // if (e == unmod_energy) {
    //     e = get_mod_energy(r_key, unique_mod_bases, mp, unmod_energy,
    //     ModLookup::Stacking);
    // }

    return e;
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
    if (vp.md.dangles == 1) {
        children_dangle_sets =
            ViennaDangles::populate_children_dangle_energies(node.children, pRNA, vp, is_external);
        closing_set = ViennaDangles::get_ml_closing_dangle_energy(node, pRNA, vp);
    }

    // Initialize energy with unmodified multiloop energy if not dangle 1
    int energy = 0;
    if (vp.md.dangles != 1) {
        energy = ViennaFunctions::multibranch_energy(node, pRNA, vp);
    }

    // Add energy corrections for modified bases in closing pair and adjacent nucleotides
    energy += update_energy(node, pRNA, mod_sequence, vp, mp, closing_set, is_external, is_closing);

    // Add energy corrections for modified bases in child loops
    for (size_t i = 0; i < node.children.size(); ++i) {
        DangleSet& current_set = (vp.md.dangles == 1) ? children_dangle_sets[i] : closing_set;
        energy += update_energy(*node.children[i], pRNA, mod_sequence, vp, mp, current_set,
                                is_external, is_not_closing);
    }

    // Compute the energy of dangle 1 of multiloop
    if (vp.md.dangles == 1) {
        energy = ViennaDangles::get_multibranch_dangle_1(node, children_dangle_sets, closing_set);
    }

    return energy;
}

int ModifiedBasesFunctions::find_mod_external_energy(
    const std::vector<std::shared_ptr<LoopNode>>& children, const ProcessedRNAEntry& pRNA,
    const std::vector<std::string_view>& mod_sequence, vrna_md_param& vp,
    const all_mod_params& mp) {
    bool is_external = true;
    std::vector<DangleSet> all_dangle_sets;
    DangleSet empty_set;  // Used when dangles != 1 to avoid checking condition in loop

    if (vp.md.dangles == 1) {
        all_dangle_sets =
            ViennaDangles::populate_children_dangle_energies(children, pRNA, vp, is_external);
    }

    // Initialize energy with unmodified external energy
    // If dangles == 1, will be replaced later with corrected dangle energies so don't double count
    int energy = 0;
    if (vp.md.dangles != 1) {
        energy = ViennaFunctions::external_energy(children, pRNA, vp);
    }

    for (size_t i = 0; i < children.size(); ++i) {
        DangleSet& current_set = (vp.md.dangles == 1) ? all_dangle_sets[i] : empty_set;
        energy += update_energy(*children[i], pRNA, mod_sequence, vp, mp, current_set, is_external);
    }

    if (vp.md.dangles == 1) {
        energy = ViennaDangles::get_external_dangle_1(children, all_dangle_sets);
    }
    return energy;
}

// Updates the energy of multiloop & external loop components
int ModifiedBasesFunctions::update_energy(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                          const std::vector<std::string_view>& mod_sequence,
                                          vrna_md_param& vp, const all_mod_params& mp,
                                          DangleSet& current_set, bool is_external,
                                          bool is_closing) {
    if (is_external && is_closing) {
        THROW_ERROR("An excternal loop cannot be a closing pair, check loop tree construction");
    }
    const std::string& sequence = pRNA.get_sequence();
    int total_energy_diff = 0;

    // Get the indices of the node and its adjacent nucleotides
    // This is to check if modified bases are present in these relevant positions.
    std::vector<size_t> indices;
    indices.push_back(node.begin);
    indices.push_back(node.end);
    if (node.begin > 0) indices.push_back(node.begin - 1);
    if (node.end + 1 < sequence.size()) indices.push_back(node.end + 1);

    // Get unique modified bases in this node
    std::vector<std::string_view> unique_mod_bases =
        unique_modified_bases_at_indices(indices, mod_sequence);
    if (unique_mod_bases.empty()) return 0;  // no modified bases in this child

    // Get pair type and encoded dangle nucleotides
    unsigned int type = ViennaUtils::get_pair_type(sequence[node.begin], sequence[node.end], vp.md);

    // Closing pair of multiloop uses inner dangles, otherwise use outer dangles
    int n5d, n3d;
    if (is_closing) {
        std::tie(n3d, n5d) = ViennaUtils::encode_inner_dangles(node.begin, node.end, pRNA, vp.md);
    } else {
        std::tie(n5d, n3d) = ViennaUtils::encode_outer_dangles(node.begin, node.end, pRNA, vp.md);
    }

    // Get modified dangle and mismatch energy differences from unmodified energies
    ModDiffs diffs = get_mod_dangle_energy_diffs(node, n5d, n3d, type, unique_mod_bases,
                                                 mod_sequence, vp, mp, is_external);

    // Apply energy corrections based on dangle settings
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
    } else if (vp.md.dangles == 1) {
        modify_dangle_set(current_set, diffs);
        return 0;
    } else {
        THROW_ERROR("Invalid dangle setting: " + std::to_string(vp.md.dangles));
    }
    return total_energy_diff;
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

// Joins string views at specified indices into a single string
std::string ModifiedBasesFunctions::join_string_views(
    std::vector<size_t> indices, const std::vector<std::string_view>& mod_sequence) {
    std::string key;
    size_t total = 0;

    // Calculate total size needed
    for (size_t idx : indices) total += mod_sequence[idx].size();
    key.reserve(total);

    // Concatenate string views
    for (size_t idx : indices) key.append(mod_sequence[idx]);
    return key;
}

void ModifiedBasesFunctions::modify_dangle_set(DangleSet& dangle_set, ModDiffs diffs) {
    dangle_set.both_dangle += diffs.mismatch;
    dangle_set.left_dangle += diffs.n5d;
    dangle_set.right_dangle += diffs.n3d;
    dangle_set.no_dangle += diffs.terminalAU;
    dangle_set.left_dangle += diffs.terminalAU;
    dangle_set.right_dangle += diffs.terminalAU;
    dangle_set.both_dangle += diffs.terminalAU;
}

// Gets the energy difference between modified and unmodified for a given key and lookup type
int ModifiedBasesFunctions::get_mod_energy_difference(
    const std::string& key, const std::vector<std::string_view>& unique_mod_bases,
    const all_mod_params& mp, int unmod_energy, ModLookup lookup_type) {
    int mod_energy = get_mod_energy(key, unique_mod_bases, mp, unmod_energy, lookup_type);
    return mod_energy - unmod_energy;
}

//
ModDiffs ModifiedBasesFunctions::get_mod_dangle_energy_diffs(
    const LoopNode& node, const int n5d, const int n3d, const unsigned int type,
    const std::vector<std::string_view>& unique_mod_bases,
    const std::vector<std::string_view>& mod_sequence, vrna_md_param& vp, const all_mod_params& mp,
    bool is_external) {
    // Stores the difference between modified and unmodified energies
    int diffTerminal = 0;  // terminal AU penalty
    int diffMM = 0;        // both neighbors (terminal mismatch)
    int diff5 = 0;         // 5' neighbor only
    int diff3 = 0;         // 3' neighbor only

    // Correct energies for dangling ends and mismatches
    if (n5d >= 0 && n3d >= 0) {
        std::string mismatch_key =
            join_string_views({node.begin, node.begin - 1, node.end, node.end + 1}, mod_sequence);
        int unmod_energy =
            is_external ? vp.p->mismatchExt[type][n5d][n3d] : vp.p->mismatchM[type][n5d][n3d];

        diffMM = get_mod_energy_difference(mismatch_key, unique_mod_bases, mp, unmod_energy,
                                           ModLookup::Mismatch);
    }

    // Dangling 5' end only
    if (n5d >= 0) {
        std::string dangle5_key =
            join_string_views({node.begin, node.end, node.begin - 1}, mod_sequence);
        int unmod_energy = vp.p->dangle5[type][n5d];
        diff5 = get_mod_energy_difference(dangle5_key, unique_mod_bases, mp, unmod_energy,
                                          ModLookup::Dangle5);
    }

    // Dangling 3' end only
    if (n3d >= 0) {
        std::string dangle3_key =
            join_string_views({node.begin, node.end, node.end + 1}, mod_sequence);
        int unmod_energy = vp.p->dangle3[type][n3d];
        diff3 = get_mod_energy_difference(dangle3_key, unique_mod_bases, mp, unmod_energy,
                                          ModLookup::Dangle3);
    }

    // Terminal AU penalty
    if (type > 2) {
        std::string terminal_key = join_string_views({node.begin, node.end}, mod_sequence);
        int unmod_energy = vp.p->TerminalAU;
        diffTerminal = get_mod_energy_difference(terminal_key, unique_mod_bases, mp, unmod_energy,
                                                 ModLookup::Terminal);
    }

    return ModDiffs(diffTerminal, diffMM, diff5, diff3);
}

}  // namespace knotergy