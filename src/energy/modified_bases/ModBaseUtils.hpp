#pragma once
#include "loop_tree/bands/BasePair.hpp"
#include "preprocessing/RNAProcessor.hpp"

#include <vector>

namespace knotergy {
/**
 * @brief Enumeration of modified base energy lookup types.
 *
 * Specifies which type of energy parameter to look up for modified bases.
 */
enum class ModLookup { Stacking, Terminal, Mismatch, Dangle5, Dangle3 };

// Stores the differences in energy contributions due to modified bases
struct ModDiffs {
    ModDiffs(int terminal_diff, int mismatch_diff, int n5d_diff, int n3d_diff)
        : terminalAU{terminal_diff}, mismatch{mismatch_diff}, n5d{n5d_diff}, n3d{n3d_diff} {}
    const int terminalAU;
    const int mismatch;
    const int n5d;
    const int n3d;
};

class ModBaseUtils {
   public:
    /**
     * @brief Join string views at specified indices into a single string.
     *
     * @param indices Vector of sequence indices.
     * @param mod_sequence The modified RNA sequence (grapheme views).
     * @return Concatenated string of bases at the specified indices.
     *
     * NOTE: Should move to general utils if used outside of modified base energy calculations
     */
    [[nodiscard]] static std::string join_string_views(
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

    /**
     * @brief Find unique modified bases at specified sequence positions.
     *
     * @param indices Vector of sequence indices to check.
     * @param mod_sequence The modified RNA sequence (grapheme views).
     * @return Vector of unique modified base string views found at those positions.
     */
    [[nodiscard]] static std::vector<std::string_view> unique_modified_bases_at_indices(
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

    [[nodiscard]] static std::vector<std::string_view> unique_modified_bases_at_inner_edge(
        size_t i, size_t j, const std::vector<std::string_view>& mod_sequence) {
        return unique_modified_bases_at_indices({i, j, i + 1, j - 1}, mod_sequence);
    }

    [[nodiscard]] static std::vector<std::string_view> unique_modified_bases_at_outer_edge(
        size_t i, size_t j, const std::vector<std::string_view>& mod_sequence) {
        std::vector<size_t> indices;
        indices.reserve(4);
        indices.push_back(i);
        indices.push_back(j);
        if (i > 0) indices.push_back(i - 1);
        if (j + 1 < mod_sequence.size()) indices.push_back(j + 1);
        return unique_modified_bases_at_indices(indices, mod_sequence);
    }

    /**
     * @brief Get modified energy or fall back to unmodified energy.
     *
     * Looks up the modified energy parameter for a given key and modified bases.
     * Returns the unmodified energy if no modified parameter is found.
     *
     * @param key The parameter key to look up.
     * @param modified Vector of modified base identifiers.
     * @param mp Vector of modified base parameters.
     * @param unmod_energy The unmodified energy to use as fallback.
     * @param lookup_type Type of energy lookup (Stacking, Terminal, etc.).
     * @return Energy value in centicalories.
     */
    [[nodiscard]] static int get_mod_energy(const std::string& key,
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
};
}  // namespace knotergy
