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