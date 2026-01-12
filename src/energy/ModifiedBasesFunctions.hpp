#pragma once

#include "../preprocessing/RNAProcessor.hpp"
#include <utility>

namespace knotergy {

    class ModifiedBasesFunctions {
    public:
            static std::pair<bool, double> find_mod_stack_energy(
                const size_t& i, const size_t& j, const size_t& ci, const size_t& cj, const std::vector<std::string_view>& mod_sequence, const std::vector<modified_base_params>& mod_params
            ) {
                bool energy_found = false;
                double energy = 0.0;
                std::string_view modified;

                // Determine if any of the four bases is unmodified and use it to find the relevant modified_base_params
                if      (!RNAProcessor::is_unmod_base(mod_sequence[i]))  modified = mod_sequence[i];
                else if (!RNAProcessor::is_unmod_base(mod_sequence[ci]))  modified = mod_sequence[ci];
                else if (!RNAProcessor::is_unmod_base(mod_sequence[j])) modified = mod_sequence[j];
                else if (!RNAProcessor::is_unmod_base(mod_sequence[cj])) modified = mod_sequence[cj];
                else return {energy_found, energy};
                

                // Construct key for lookup
                std::string key;
                key.reserve(mod_sequence[i].size() + mod_sequence[ci].size() +
                            mod_sequence[cj].size() + mod_sequence[j].size());
                key.append(mod_sequence[i]);
                key.append(mod_sequence[ci]);
                key.append(mod_sequence[j]);
                key.append(mod_sequence[cj]);

                // Search through modified base parameters to find stacking energy
                for (const modified_base_params& param : mod_params) {
                    if (param.modified_base != modified) continue;

                    auto it = param.stacking_energies.find(key);
                    if (it != param.stacking_energies.end()) {
                        energy_found = true;
                        energy = static_cast<double>(it->second);
                        return {energy_found, energy * 100};
                    }
                }
                
                return {energy_found, energy};
            }
    };



} // namespace knotergy