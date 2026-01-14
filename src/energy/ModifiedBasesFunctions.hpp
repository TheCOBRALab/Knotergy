#pragma once

#include "../preprocessing/RNAProcessor.hpp"
#include "./ViennaFunctions.hpp"
#include <utility>

namespace knotergy {

    class ModifiedBasesFunctions {
    public:
        static double find_mod_stack_energy(
            const size_t& i, const size_t& j, const size_t& ci, const size_t& cj, const std::vector<std::string_view>& mod_sequence, const std::vector<modified_base_params>& mod_params
        ) {
            double energy = INF;

            // Find all modified bases involved in the stack
            std::vector<std::string_view> modified = modified_bases({i, j, ci, cj}, mod_sequence);
            if (modified.empty()) return energy; // Return INF
            
            // Construct key for lookup
            std::string key = join_string_views({i, ci, j, cj}, mod_sequence);

            // Search through modified base parameters to find stacking energy
            for (const modified_base_params& param : mod_params) {
                // Skip if this modified base is not involved in the current stack
                if (std::find(modified.begin(), modified.end(), param.modified_base) == modified.end()) continue;

                auto it = param.stacking_energies.find(key);
                if (it != param.stacking_energies.end()) {
                    // TODO: Consider checking for duplicates in other modified_base_params entries
                    energy = static_cast<double>(it->second);
                    return energy * 100;
                }
            }
            
            return energy; // return INF if no modified stacking energy found
        }

        static double find_mod_external_energy(const std::vector<std::shared_ptr<LoopNode>>& children,
                                               const std::string& sequence,
                                               const std::vector<std::string_view>& mod_sequence,
                                               const std::vector<modified_base_params>& mod_params) {
            double energy = 0.0;
            // Collect all modified bases in the external loop

            if (ViennaParams::md.dangles == 1) {
                    // TODO: Implement modified base handling for dangle 1 in external loops
            }

            for (std::shared_ptr<LoopNode> c : children) {
                if (c->loop_type != LoopType::Pseudoknot) {

                    // sanity check for indices
                    if (c->begin >= sequence.size() - 1 || c->end >= sequence.size()) {
                        THROW_ERROR("Invalid indices for external loop energy calculation.");
                    }

                    // Get pair type and nucleotides
                    unsigned int pair_type = ViennaUtils::get_pair_type(sequence[c->begin], sequence[c->end], ViennaParams::md);
                    int n5d = vrna_nucleotide_encode(sequence[c->begin - 1], &ViennaParams::md);
                    int n3d = vrna_nucleotide_encode(sequence[c->end + 1], &ViennaParams::md);

                    if (ViennaParams::md.dangles == 0) {
                        n5d = -1;
                        n3d = -1;
                    }

                    std::vector<std::string_view> modified = ModifiedBasesFunctions::modified_bases({c->begin, c->end}, mod_sequence);

                    return mod_exterior_stem(
                        pair_type,
                        n5d,
                        n3d,
                        modified,
                        mod_sequence,
                        mod_params,
                        ViennaParams::P
                    );
                }
            }

            return energy; // return INF if no modified external energy found
            
        }

        static double mod_exterior_stem(unsigned int type,
                                        int          n5d,
                                        int          n3d,
                                        const std::vector<std::string_view>& modified, 
                                        const std::vector<std::string_view>& mod_sequence,
                                        const std::vector<modified_base_params>& mod_params,
                                        vrna_param_t *p = ViennaParams::P) {
            double energy = 0.0;
            // Base Case: No modified bases involved
            std::vector<std::string_view> modified = ModifiedBasesFunctions::modified_bases({c->begin, c->end}, mod_sequence);
            if (modified.empty()) {
                return vrna_E_exterior_stem(type, n5d, n3d, p);
            }

            std::string key = join_string_views({c->begin, c->end}, mod_sequence);
            for (const modified_base_params& param : mod_params) {
                // Skip if this modified base is not involved in the current pair
                if (std::find(modified.begin(), modified.end(), param.modified_base) == modified.end()) continue;

                // Check mismatch energies
                if (n5d >=0 && n3d >=0) {
                    auto mismatch_it = param.mismatch_energies.find(key);
                    if (mismatch_it != param.mismatch_energies.end()) {
                        energy += static_cast<int>(mismatch_it->second * 100);
                    } else {
                        energy += ViennaParams::P->mismatchExt[type][n5d][n3d];
                    }
                } 
                // Check 5' dangle energies
                else if (n5d >=0) {
                    
                    auto dangle5_it = param.dangle5_energies.find(key);
                    if (dangle5_it != param.dangle5_energies.end()) {
                        energy += static_cast<int>(dangle5_it->second * 100);
                    } else {
                        energy += ViennaParams::P->dangle5[type][n5d];
                    }
                } 
                // Check 3' dangle energies
                else if (n3d >=0) {
                    auto dangle3_it = param.dangle3_energies.find(key);
                    if (dangle3_it != param.dangle3_energies.end()) {
                        energy += static_cast<int>(dangle3_it->second * 100);
                    } else {
                        energy += ViennaParams::P->dangle3[type][n3d];
                    }
                }

                return energy;
            }
        }

    private:
        static std::vector<std::string_view> modified_bases(std::initializer_list<size_t> indices, const std::vector<std::string_view>& mod_sequence) {
            std::vector<std::string_view> modified;
            modified.reserve(indices.size());
            for (size_t idx : indices) {
                // Check if base is modified and not already in the list
                if (!RNAProcessor::is_unmod_base(mod_sequence[idx]) && std::find(modified.begin(), modified.end(), mod_sequence[idx]) == modified.end()) {
                    modified.push_back(mod_sequence[idx]);
                }
            }
            return modified;
        }

        static std::string join_string_views(std::initializer_list<size_t> indices, const std::vector<std::string_view>& mod_sequence) {
            std::string key;
            size_t total = 0;

            // Calculate total size needed
            for (size_t idx : indices) total += mod_sequence[idx].size();
            key.reserve(total);

            // Concatenate string views
            for (size_t idx : indices) key.append(mod_sequence[idx]);
            return key;
        }
    };





} // namespace knotergy