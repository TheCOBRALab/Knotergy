#include "ModifiedBasesFunctions.hpp"

namespace knotergy {

int ModifiedBasesFunctions::find_mod_stack_energy(
    const size_t& i, const size_t& j, const size_t& ci, const size_t& cj, std::string sequence,
    const std::vector<std::string_view>& mod_sequence,
    const std::vector<modified_base_params>& mod_params) {

    // Get Vienna stacking energy for unmodified bases
    int unmod_energy = ViennaFunctions::stack_energy(i, j, ci, cj, sequence);
    
    std::vector<std::string_view> unique_mod_bases = unique_modified_bases_at_indices({i, j, ci, cj}, mod_sequence);
    if (unique_mod_bases.empty()) return unmod_energy;

    // Used to look up stacking energies in modified base parameters
    std::string key_right = join_string_views({cj, j,  ci, i }, mod_sequence);
    std::string key_left  = join_string_views({i,  ci, j,  cj}, mod_sequence);
    
    // // This makes no logical sense, but tests show that ViennaRNA uses either orientation if there's a modified base on the 5' side
    // // and only left orientation if modified base is only on 3' side
    // if (!unique_modified_bases_at_indices({i, ci}, mod_sequence).empty()) {
    //     double e = get_mod_energy(key_right, unique_mod_bases, mod_params, unmod_energy, ModLookup::Stacking);
    //     if (e == unmod_energy) {
    //         e = get_mod_energy(key_left, unique_mod_bases, mod_params, unmod_energy, ModLookup::Stacking);
    //     }
    //     return e;
    // }
    // return get_mod_energy(key_right, unique_mod_bases, mod_params, unmod_energy, ModLookup::Stacking);
    

    
    // Get energy correction for modified bases (returns original energy if no modifications found)
    // Try the Vienna-defined ordering first, then the swapped orientation
    int e = get_mod_energy(key_right, unique_mod_bases, mod_params, unmod_energy, ModLookup::Stacking);
    if (e == unmod_energy) {
        e = get_mod_energy(key_left, unique_mod_bases, mod_params, unmod_energy, ModLookup::Stacking);
    }

    return e;
}

int ModifiedBasesFunctions::find_mod_external_energy(
    const std::vector<std::shared_ptr<LoopNode>>& children, const std::string& sequence,
    const std::vector<std::string_view>& mod_sequence,
    const std::vector<modified_base_params>& mod_params) {
    
    bool is_external = true;
    std::vector<DangleSet> dangle_energies;

    if (ViennaParams::md.dangles == 1) {
        dangle_energies = ViennaDangles::populate_children_dangle_energies(children, sequence, is_external);
    }   
    
    // Initialize energy with unmodified external energy
    // If dangles == 1, will be replaced later with corrected dangle energies so don't double count
     int energy;
     if (ViennaParams::md.dangles != 1) {
         energy = ViennaFunctions::external_energy(children, sequence);
     }

     for (size_t i = 0; i < children.size(); ++i){
        std::shared_ptr<LoopNode> c = children[i];

        // Gather indices to check for modified bases
        std::vector<size_t> indices = {c->begin, c->end};
        if (c->begin > 0) indices.push_back(c->begin - 1);
        if (c->end + 1 < sequence.size()) indices.push_back(c->end + 1);

        // Get unique modified bases in this child
        std::vector<std::string_view> unique_mod_bases = unique_modified_bases_at_indices(indices, mod_sequence);
        if (unique_mod_bases.empty()) continue; // no modified bases in this child

        // Get pair type and encoded dangle nucleotides
        unsigned int type = ViennaUtils::get_pair_type(sequence[c->begin], sequence[c->end]);
        auto [n5d, n3d] = ViennaUtils::encode_outer_dangles(c->begin, c->end, sequence);
        
        // --------------------------- Dangles ---------------------------
        if (ViennaParams::md.dangles == 0) {
            n5d = -1; n3d = -1;
        } else if (ViennaParams::md.dangles == 1) {
            DangleSet& current_set = dangle_energies[i];
            modify_dangle_set(current_set, c, type, n5d, n3d, unique_mod_bases, mod_sequence, mod_params, is_external);
            continue; // modified energies handled in modify_dangle_set
        } 

        // --------------------------- Dangle 0 or 2 ---------------------------
        // Correct energies for dangling ends and mismatches
        if (n5d >=0 && n3d >=0) {
            // Dangle key: XYZW = X pairs with Z, Y is 5' dangle, W is 3' dangle
            std::string mismatch_key = join_string_views({c->begin, c->begin - 1, c->end, c->end + 1}, mod_sequence);
            energy += get_mod_energy_difference(mismatch_key, unique_mod_bases, mod_params, ViennaParams::p->mismatchExt[type][n5d][n3d], ModLookup::Mismatch);
        } else if (n5d >=0) {
            // 5' Dangle key: XYZ = XY are the pair, Z is 5' dangle
            std::string dangle5_key = join_string_views({c->begin, c->end, c->begin - 1}, mod_sequence);
            energy += get_mod_energy_difference(dangle5_key, unique_mod_bases, mod_params, ViennaParams::p->dangle5[type][n5d], ModLookup::Dangle5);
        } else if (n3d >=0) {
            // 3' Dangle key: XYZ = XY are the pair, Z is 3' dangle
            std::string dangle3_key = join_string_views({c->begin, c->end, c->end + 1}, mod_sequence);
            energy += get_mod_energy_difference(dangle3_key, unique_mod_bases, mod_params, ViennaParams::p->dangle3[type][n3d], ModLookup::Dangle3);
        }

        if (type > 2) {
            // Terminal key: XY = X pairs with Y
            std::string terminal_key = join_string_views({c->begin, c->end}, mod_sequence);
            energy += get_mod_energy_difference(terminal_key, unique_mod_bases, mod_params, ViennaParams::p->TerminalAU, ModLookup::Terminal);
        }
     }

     if (ViennaParams::md.dangles == 1) {
         energy = ViennaDangles::get_external_dangle_1(children, dangle_energies);
     }
     return energy;
}

int ModifiedBasesFunctions::get_mod_energy(
    const std::string& key, const std::vector<std::string_view>& unique_mod_bases,
    const std::vector<modified_base_params>& mod_params, int unmod_energy, ModLookup lookup_type) {
    if (unique_mod_bases.empty()) {
        return static_cast<int>(unmod_energy);
    }
    
    unsigned int modified_found = 0;
    for (const modified_base_params& param : mod_params) {
        if (std::find(unique_mod_bases.begin(), unique_mod_bases.end(), param.modified_base) == unique_mod_bases.end()) {
            continue;
        }
        ++modified_found;

        const std::map<std::string, float>* energy_lookup = nullptr;
        switch (lookup_type) {
            case ModLookup::Stacking: energy_lookup = &param.stacking_energies; break;
            case ModLookup::Terminal: energy_lookup = &param.terminal_energies; break;
            case ModLookup::Mismatch: energy_lookup = &param.mismatch_energies; break;
            case ModLookup::Dangle5:  energy_lookup = &param.dangle5_energies; break;
            case ModLookup::Dangle3:  energy_lookup = &param.dangle3_energies; break;
            default:
                THROW_ERROR("Invalid ModLookup type: " + std::to_string(static_cast<int>(lookup_type)));
                break;
        }

        // look for matching energy entry
        if (energy_lookup) {
            if (auto it = energy_lookup->find(key); it != energy_lookup->end()) {
                std::cout << "Modified base energy found for key: " << key << " -> " << it->second << " Diff: " << static_cast<int>(it->second * 100) - unmod_energy << std::endl;
                return static_cast<int>(it->second * 100);
            }
        }

        // early exit if all modified bases have been found and no matches exist
        if (modified_found >= unique_mod_bases.size()) {
            break;
        }
    }
    std::cout << "No modified base energy found for key: " << key << ", using unmodified energy: " << unmod_energy << std::endl;
    std::cout << "ModLookup type: " << static_cast<int>(lookup_type) << std::endl;
    return static_cast<int>(unmod_energy);
}

int ModifiedBasesFunctions::get_mod_energy_difference(
    const std::string& key, const std::vector<std::string_view>& unique_mod_bases,
    const std::vector<modified_base_params>& mod_params, int unmod_energy, ModLookup lookup_type) {
    int mod_energy = get_mod_energy(key, unique_mod_bases, mod_params, unmod_energy, lookup_type);
    return mod_energy - unmod_energy;
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

void ModifiedBasesFunctions::modify_dangle_set( DangleSet& dangle_set,
                                    const std::shared_ptr<LoopNode>& c, unsigned int type,
                                    int n5d, int n3d, std::vector<std::string_view> unique_mod_bases,
                                    const std::vector<std::string_view>& mod_sequence,
                                    const std::vector<modified_base_params>& mod_params,
                                    bool is_external){

    if (n5d >= 0 && n3d >= 0) {
        std::string mismatch_key = join_string_views({c->begin, c->begin - 1, c->end, c->end + 1}, mod_sequence);
        int unmod_energy = is_external ? ViennaParams::p->mismatchExt[type][n5d][n3d] : ViennaParams::p->mismatchM[type][n5d][n3d];
        dangle_set.both_dangle = get_mod_energy_difference(mismatch_key, unique_mod_bases, mod_params, unmod_energy, ModLookup::Mismatch);
    }
    else if (n5d >= 0) {
        std::string dangle5_key = join_string_views({c->begin, c->end, c->begin - 1}, mod_sequence);
        int unmod_energy = is_external ? ViennaParams::p->dangle5[type][n5d] : ViennaParams::p->dangle5[type][n5d];
        dangle_set.left_dangle = get_mod_energy_difference(dangle5_key, unique_mod_bases, mod_params, unmod_energy, ModLookup::Dangle5);
    }
    else if (n3d >= 0) {
        std::string dangle3_key = join_string_views({c->begin, c->end, c->end + 1}, mod_sequence);
        int unmod_energy = is_external ? ViennaParams::p->dangle3[type][n3d] : ViennaParams::p->dangle3[type][n3d];
        dangle_set.right_dangle = get_mod_energy_difference(dangle3_key, unique_mod_bases, mod_params, unmod_energy, ModLookup::Dangle3);
    }
    if (type > 2) {
        std::string terminal_key = join_string_views({c->begin, c->end}, mod_sequence);
        int unmod_energy = is_external ? ViennaParams::p->TerminalAU : ViennaParams::p->TerminalAU;
        int terminal_au_diff = get_mod_energy_difference(terminal_key, unique_mod_bases, mod_params, unmod_energy, ModLookup::Terminal);
        
        dangle_set.no_dangle += terminal_au_diff;
        dangle_set.left_dangle += terminal_au_diff;
        dangle_set.right_dangle += terminal_au_diff;
        dangle_set.both_dangle += terminal_au_diff;
    }

}

}