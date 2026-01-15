#pragma once

#include "../preprocessing/RNAProcessor.hpp"
#include "./ViennaFunctions.hpp"
#include <utility>

namespace knotergy {

    class ModifiedBasesFunctions {
    public:
        static double find_mod_stack_energy(
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

    static double find_mod_external_energy(const std::vector<std::shared_ptr<LoopNode>>& children, const std::string& sequence, std::vector<std::string_view> mod_sequence, const std::vector<modified_base_params>& mod_params) {
            if (ViennaParams::md.dangles == 1) {
                return INF; // TODO
            }
            
            double energy = ViennaFunctions::external_energy(children, sequence);
            for (std::shared_ptr<LoopNode> c : children) {
                std::vector<std::string_view> modified = ModifiedBasesFunctions::modified_bases({c->begin, c->end}, mod_sequence);
                if (modified.empty()) continue;

                std::string key = ModifiedBasesFunctions::join_string_views({c->begin, c->end}, mod_sequence);

                unsigned int type = ViennaUtils::get_pair_type(sequence[c->begin], sequence[c->end]);
                int n5d = c->begin > 0 ? vrna_nucleotide_encode(sequence[c->begin - 1], &ViennaParams::md) : -1;
                int n3d = c->end < sequence.size() - 1 ? vrna_nucleotide_encode(sequence[c->end + 1], &ViennaParams::md) : -1;

                if (ViennaParams::md.dangles == 0) {
                    n5d = -1;
                    n3d = -1;
                }

                if (n5d >= 0 && n3d >= 0) {
                    for (const modified_base_params& param : mod_params) {
                        if (std::find(modified.begin(), modified.end(), param.modified_base) == modified.end()) continue;

                        auto it = param.mismatch_energies.find(key);
                        if (it != param.mismatch_energies.end()) {
                            energy += static_cast<double>(it->second) - ViennaParams::P->mismatchExt[type][n5d][n3d];
                        }
                    }
                } else if (n5d >= 0) {
                    for (const modified_base_params& param : mod_params) {
                        if (std::find(modified.begin(), modified.end(), param.modified_base) == modified.end()) continue;

                        auto it = param.dangle5_energies.find(key);
                        if (it != param.dangle5_energies.end()) {
                            energy += static_cast<double>(it->second) - ViennaParams::P->dangle5[type][n5d];
                        }
                    }
                } else if (n3d >= 0) {
                    for (const modified_base_params& param : mod_params) {
                        if (std::find(modified.begin(), modified.end(), param.modified_base) == modified.end()) continue;

                        auto it = param.dangle3_energies.find(key);
                        if (it != param.dangle3_energies.end()) {
                            energy += static_cast<double>(it->second) - ViennaParams::P->dangle3[type][n3d];
                        }
                    }
                }
            }

            return energy;
        }


    private:
        // Returns a vector of unique modified bases found at the specified indices
        static std::vector<std::string_view> modified_bases(std::initializer_list<size_t> indices, const std::vector<std::string_view>& mod_sequence) {
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