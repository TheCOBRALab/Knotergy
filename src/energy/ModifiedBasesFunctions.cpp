#include "ModifiedBasesFunctions.hpp"

namespace knotergy
{

double ModifiedBasesFunctions::find_mod_stack_energy(
    const size_t& i, const size_t& j, const size_t& ci, const size_t& cj, std::string sequence, const std::vector<std::string_view>& mod_sequence, const std::vector<modified_base_params>& mod_params
) {
    double energy = ViennaFunctions::stack_energy(i, j, ci, cj, sequence);

    std::vector<std::string_view> modified = ModifiedBasesFunctions::modified_bases({i, j, ci, cj}, mod_sequence);
    if (modified.empty()) return energy;
    
    // Stacking key: XYZW = (X–Y) stacked on (Z–W); e.g. "PUAA" => P–U on A–A
    std::string key = ModifiedBasesFunctions::join_string_views({i, ci, j, cj}, mod_sequence);

    for (const modified_base_params& param : mod_params) {
        if (std::find(modified.begin(), modified.end(), param.modified_base) == modified.end()) continue;

        auto it = param.stacking_energies.find(key);
        if (it != param.stacking_energies.end()) {
            // TODO: Consider checking for duplicates in other modified_base_params entries
            energy = static_cast<double>(it->second);
            return energy * 100;
        }
    }
    return energy;
}

double ModifiedBasesFunctions::find_mod_external_energy(const std::vector<std::shared_ptr<LoopNode>>& children, const std::string& sequence, std::vector<std::string_view> mod_sequence, const std::vector<modified_base_params>& mod_params) {
    
    // --------------- Dangle type 1 ---------------
    if (ViennaParams::md.dangles == 1) {
        bool is_external = true;
        std::vector<DangleSet> dangle_energies = ViennaDangles::populate_children_dangle_energies(children, sequence, is_external);
        
        // Update existing dangle energies with modified base corrections
        for (size_t i = 0; i < children.size(); ++i) {
            std::shared_ptr<LoopNode> c = children[i];
            DangleSet& current_set = dangle_energies[i];
            std::vector<std::string_view> modified = ModifiedBasesFunctions::modified_bases({c->begin, c->end}, mod_sequence);

            unsigned int type = ViennaUtils::get_pair_type(sequence[c->begin], sequence[c->end]);
            int n5d = c->begin > 0 ? vrna_nucleotide_encode(sequence[c->begin - 1], &ViennaParams::md) : -1;
            int n3d = c->end < sequence.size() - 1 ? vrna_nucleotide_encode(sequence[c->end + 1], &ViennaParams::md) : -1;
            
            dangle_energies[i] = ModifiedBasesFunctions::correct_dangle_set(current_set, c, type, n5d, n3d, modified, mod_sequence, mod_params, is_external);
        }
        
        return ViennaDangles::get_external_dangle_1(children, dangle_energies);
    }
    
    // --------------- No dangles or dangle type 2 ---------------
    double energy = ViennaFunctions::external_energy(children, sequence);
    for (std::shared_ptr<LoopNode> c : children) {
        std::vector<std::string_view> modified = ModifiedBasesFunctions::modified_bases({c->begin, c->end}, mod_sequence);
        if (modified.empty()) continue;
        
        unsigned int type = ViennaUtils::get_pair_type(sequence[c->begin], sequence[c->end]);
        int n5d = c->begin > 0 ? vrna_nucleotide_encode(sequence[c->begin - 1], &ViennaParams::md) : -1;
        int n3d = c->end < sequence.size() - 1 ? vrna_nucleotide_encode(sequence[c->end + 1], &ViennaParams::md) : -1;

        if (ViennaParams::md.dangles == 0) {
            n5d = -1;
            n3d = -1;
        }

        if (n5d >= 0 && n3d >= 0) {
            std::string mismatch_key = ModifiedBasesFunctions::join_string_views({c->begin, c->begin-1, c->end, c->end+1}, mod_sequence);
            energy += ModifiedBasesFunctions::get_corrected_energy(mismatch_key, modified, mod_params, ViennaParams::p->mismatchExt[type][n5d][n3d]);
        } else if (n5d >= 0) {
            std::string dangle5_key  = ModifiedBasesFunctions::join_string_views({c->begin, c->end, c->begin - 1}, mod_sequence);
            energy += ModifiedBasesFunctions::get_corrected_energy(dangle5_key, modified, mod_params, ViennaParams::p->dangle5[type][n5d]);
        } else if (n3d >= 0) {
            std::string dangle3_key  = ModifiedBasesFunctions::join_string_views({c->begin, c->end, c->end + 1}, mod_sequence);
            energy += ModifiedBasesFunctions::get_corrected_energy(dangle3_key, modified, mod_params, ViennaParams::p->dangle3[type][n3d]);
        }

        if (type > 2){
            std::string terminal_key = ModifiedBasesFunctions::join_string_views({c->end, c->begin}, mod_sequence);
            energy += ModifiedBasesFunctions::get_corrected_energy(terminal_key, modified, mod_params, ViennaParams::p->TerminalAU);
        }
    }

    return energy;
}


double ModifiedBasesFunctions::find_mod_multibranch_energy(const LoopNode& node, const std::string& sequence, std::vector<std::string_view> mod_sequence, const std::vector<modified_base_params>& mod_params) {
    
    std::vector<std::string_view> modified = ModifiedBasesFunctions::modified_bases({node.begin, node.begin + 1, node.end, node.end - 1}, mod_sequence);
    double energy = 0.0;
    // ------------------ Dangle 1 Energy ------------------
    if (ViennaParams::md.dangles == 1) {
        bool is_external = false;

        energy = ViennaParams::p->MLclosing; // closing penalty
        energy += node.exclusive_unpaired_bases_count * ViennaParams::p->MLbase; // unpaired bases penalty

        DangleSet closing = ViennaDangles::get_ml_dangle_energy(node, sequence);
        closing = ModifiedBasesFunctions::correct_dangle_set(closing, std::make_shared<LoopNode>(node), 
                                                            ViennaUtils::reverse_pair_type(sequence[node.begin], sequence[node.end]),
                                                            vrna_nucleotide_encode(sequence[node.begin + 1], &ViennaParams::md),
                                                            vrna_nucleotide_encode(sequence[node.end - 1], &ViennaParams::md),
                                                            modified,
                                                            mod_sequence,
                                                            mod_params,
                                                            is_external);


        std::vector<DangleSet> dangle_energies = ViennaDangles::populate_children_dangle_energies(node.children, sequence);
        // Update existing dangle energies with modified base corrections
        for (size_t i = 0; i < node.children.size(); ++i) {
            std::shared_ptr<LoopNode> c = node.children[i];
            DangleSet& current_set = dangle_energies[i];
            std::vector<std::string_view> c_modified = ModifiedBasesFunctions::modified_bases({c->begin, c->end}, mod_sequence);

            unsigned int type = ViennaUtils::get_pair_type(sequence[c->begin], sequence[c->end]);
            int n5d = c->begin > 0 ? vrna_nucleotide_encode(sequence[c->begin - 1], &ViennaParams::md) : -1;
            int n3d = c->end < sequence.size() - 1 ? vrna_nucleotide_encode(sequence[c->end + 1], &ViennaParams::md) : -1;

            dangle_energies[i] = ModifiedBasesFunctions::correct_dangle_set(current_set, c, type, n5d, n3d, c_modified, mod_sequence, mod_params, is_external);
        }

        return ViennaDangles::get_multibranch_dangle_1(node, dangle_energies, closing) + energy;
    }

    energy += ViennaFunctions::multibranch_energy(node, sequence); // base multibranch energy

    size_t i = node.begin;
    size_t j = node.end;

    unsigned int pair_type = ViennaUtils::reverse_pair_type(sequence[i], sequence[j]);
    auto [n5d, n3d] = ViennaUtils::encode_nucleotides(sequence[i + 1], sequence[j - 1]);

    if (ViennaParams::md.dangles == 0) {
        n5d = -1;
        n3d = -1;
    }

    energy += get_corrected_branch_energy(node, pair_type, n3d, n5d, ModifiedBasesFunctions::modified_bases({i + 1, j - 1}, mod_sequence), mod_sequence, mod_params);

    for (const std::shared_ptr<LoopNode>& child : node.children) {
        size_t ci = child->begin;
        size_t cj = child->end;
        std::vector<std::string_view> c_modified = ModifiedBasesFunctions::modified_bases({ci, ci - 1, cj, cj + 1}, mod_sequence);
        unsigned int child_pair_type = ViennaUtils::get_pair_type(sequence[ci], sequence[cj]);
        auto [child_n5d, child_n3d]  = ViennaUtils::encode_nucleotides(sequence[ci - 1], sequence[cj + 1]);

        if (ViennaParams::md.dangles == 0) {
            child_n5d = -1;
            child_n3d = -1;
        }

        energy += get_corrected_branch_energy(*child, child_pair_type, child_n5d, child_n3d, c_modified, mod_sequence, mod_params);
    }

    return energy;

}


// Returns a vector of unique modified bases found at the specified indices
std::vector<std::string_view> ModifiedBasesFunctions::modified_bases(std::initializer_list<size_t> indices, const std::vector<std::string_view>& mod_sequence) {
    std::vector<std::string_view> modified;
    modified.reserve(indices.size());
    for (size_t idx : indices) {
        if (!RNAProcessor::is_unmod_base(mod_sequence[idx]) && std::find(modified.begin(), modified.end(), mod_sequence[idx]) == modified.end()) {
            modified.push_back(mod_sequence[idx]);
        }
    }
    return modified;
}

// Joins string views at specified indices into a single string
std::string ModifiedBasesFunctions::join_string_views(std::initializer_list<size_t> indices, const std::vector<std::string_view>& mod_sequence) {
    std::string key;
    size_t total = 0;

    // Calculate total size needed
    for (size_t idx : indices) total += mod_sequence[idx].size();
    key.reserve(total);

    // Concatenate string views
    for (size_t idx : indices) key.append(mod_sequence[idx]);
    return key;
}

double ModifiedBasesFunctions::get_corrected_energy(const std::string& key, const std::vector<std::string_view>& modified, const std::vector<modified_base_params>& mod_params, double vienna_energy) {
    unsigned int modified_found = 0;
    for (const modified_base_params& param : mod_params) {
        if (std::find(modified.begin(), modified.end(), param.modified_base) == modified.end()) continue;
        ++modified_found;

        // look for matching energy entry
        auto it = param.stacking_energies.find(key);
        if (it != param.stacking_energies.end()) {
            return static_cast<double>(it->second) - vienna_energy;
        }

        // early exit if all modified bases have been found and no matches exist
        if (modified_found >= modified.size()) {
            break;
        }
    }

    return 0.0;
}

double ModifiedBasesFunctions::get_corrected_branch_energy(const LoopNode& node, unsigned int type, int n5d, int n3d, std::vector<std::string_view> modified, const std::vector<std::string_view>& mod_sequence, const std::vector<modified_base_params>& mod_params) {
    double energy = 0.0;
    if (n5d >= 0 && n3d >= 0) {
        std::string mismatch_key = ModifiedBasesFunctions::join_string_views({node.begin, node.begin-1, node.end, node.end+1}, mod_sequence);
        energy += ModifiedBasesFunctions::get_corrected_energy(mismatch_key, modified, mod_params, ViennaParams::p->mismatchM[type][n5d][n3d]);
    } else if (n5d >= 0) {
        std::string dangle5_key  = ModifiedBasesFunctions::join_string_views({node.begin, node.end, node.begin - 1}, mod_sequence);
        energy += ModifiedBasesFunctions::get_corrected_energy(dangle5_key, modified, mod_params, ViennaParams::p->dangle5[type][n5d]);
    } else if (n3d >= 0) {
        std::string dangle3_key  = ModifiedBasesFunctions::join_string_views({node.begin, node.end, node.end + 1}, mod_sequence);
        energy += ModifiedBasesFunctions::get_corrected_energy(dangle3_key, modified, mod_params, ViennaParams::p->dangle3[type][n3d]);
    }

    if (type > 2){
        std::string terminal_key = ModifiedBasesFunctions::join_string_views({node.end, node.begin}, mod_sequence);
        energy += ModifiedBasesFunctions::get_corrected_energy(terminal_key, modified, mod_params, ViennaParams::p->TerminalAU);
    }
    return energy;
}

DangleSet ModifiedBasesFunctions::correct_dangle_set(const DangleSet& original, const std::shared_ptr<LoopNode>& c, unsigned int type, int n5d, int n3d, std::vector<std::string_view> modified,  const std::vector<std::string_view>& mod_sequence, const std::vector<modified_base_params>& mod_params, bool is_external) {
    DangleSet corrected = original;
    
    if (n5d >= 0 && n3d >= 0) {
        std::string mismatch_key = ModifiedBasesFunctions::join_string_views({c->begin, c->begin-1, c->end, c->end+1}, mod_sequence);
        int vienna_energy = is_external ? ViennaParams::p->mismatchExt[type][n5d][n3d] : ViennaParams::p->mismatchM[type][n5d][n3d];
        corrected.both_dangle += static_cast<int>(ModifiedBasesFunctions::get_corrected_energy(mismatch_key, modified, mod_params, vienna_energy) * 100);
    } 
    if (n5d >= 0) {
        std::string dangle5_key  = ModifiedBasesFunctions::join_string_views({c->begin, c->end, c->begin - 1}, mod_sequence);
        corrected.left_dangle += static_cast<int>(ModifiedBasesFunctions::get_corrected_energy(dangle5_key, modified, mod_params, ViennaParams::p->dangle5[type][n5d]) * 100);
    } 
    if (n3d >= 0) {
        std::string dangle3_key  = ModifiedBasesFunctions::join_string_views({c->begin, c->end, c->end + 1}, mod_sequence);
        corrected.right_dangle += static_cast<int>(ModifiedBasesFunctions::get_corrected_energy(dangle3_key, modified, mod_params, ViennaParams::p->dangle3[type][n3d]) * 100);
    }

    if (type > 2){
        std::string terminal_key = ModifiedBasesFunctions::join_string_views({c->end, c->begin}, mod_sequence);
        int terminal_au = static_cast<int>(ModifiedBasesFunctions::get_corrected_energy(terminal_key, modified, mod_params, ViennaParams::p->TerminalAU) * 100);
        corrected.no_dangle += terminal_au;
        corrected.left_dangle += terminal_au;
        corrected.right_dangle += terminal_au;
        corrected.both_dangle += terminal_au;
    }

    return corrected;

}

} // namespace knotergy
    
