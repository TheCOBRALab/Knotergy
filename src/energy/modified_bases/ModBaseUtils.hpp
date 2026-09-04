#pragma once
#include "energy/dangles/Dangle1.hpp"
#include "energy/vienna/ViennaUtils.hpp"
#include "preprocessing/RNAProcessor.hpp"
#include "utils/common.hpp"

namespace viennarna = thermorna::viennarna;

#include <vector>

namespace knotergy {
/**
 * @brief Enumeration of modified base energy lookup types.
 *
 * Specifies which type of energy parameter to look up for modified bases.
 */
enum class ModLookup { Stacking, TerminalAU, Mismatch, Dangle5, Dangle3 };

// Stores the differences in energy contributions due to modified bases
struct ModDiffs {
    ModDiffs() : terminalAU{0}, mismatch{0}, n5d{0}, n3d{0} {}
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
        std::vector<std::size_t> indices, const std::vector<std::string_view>& mod_sequence) {
        std::string key;
        std::size_t total = 0;

        // Calculate total size needed
        for (std::size_t idx : indices) {
            if (idx >= mod_sequence.size()) {
                THROW_ERROR("Index " + std::to_string(idx) +
                            " is out of bounds for modified sequence of size " +
                            std::to_string(mod_sequence.size()));
            }
            total += mod_sequence[idx].size();
        }
        key.reserve(total);

        // Concatenate string views
        for (std::size_t idx : indices) {
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
        std::vector<std::size_t> indices, const std::vector<std::string_view>& mod_sequence) {
        std::vector<std::string_view> modified;
        if (mod_sequence.empty())
            return modified;  // no modified sequence provided, return empty vector
        modified.reserve(indices.size());
        for (std::size_t idx : indices) {
            if (!RNAProcessor::is_unmodified_base(mod_sequence[idx]) &&
                std::find(modified.begin(), modified.end(), mod_sequence[idx]) == modified.end()) {
                modified.push_back(mod_sequence[idx]);
            }
        }
        return modified;
    }

    [[nodiscard]] static std::vector<std::string_view> unique_mod_bases_in_string(
        const std::string& str) {
        std::vector<std::string_view> graphemes = ProcessedRNAEntry::parse_modified_sequence(str);
        std::vector<std::string_view> unique_mod_bases;

        for (const std::string_view& g : graphemes) {
            if (!RNAProcessor::is_unmodified_base(g) &&
                std::find(unique_mod_bases.begin(), unique_mod_bases.end(), g) ==
                    unique_mod_bases.end()) {
                unique_mod_bases.push_back(g);
            }
        }
        return unique_mod_bases;
    }

    [[nodiscard]] static std::vector<std::string_view> unique_mod_bases_at_inner_edge(
        std::size_t i, std::size_t j, const std::vector<std::string_view>& mod_sequence) {
        return unique_modified_bases_at_indices({i, j, i + 1, j - 1}, mod_sequence);
    }

    [[nodiscard]] static std::vector<std::string_view> unique_modified_bases_at_outer_edge(
        std::size_t i, std::size_t j, const std::vector<std::string_view>& mod_sequence) {
        std::vector<std::size_t> indices;
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
                                            const all_mod_params& mp, ModLookup lookup_type) {
        // If no modified bases are present, return the unmodified energy
        if (unique_mod_bases.empty()) {
            return NULL_ENERGY;
        }

        // Get the pointer to the correct energy map based on lookup type
        for (const std::string_view& mod_base : unique_mod_bases) {
            const modified_base_param* param = mp.get_modified_base_param(std::string(mod_base));

            if (!param) {
                THROW_ERROR("Modified base '" + std::string(mod_base) +
                            "' found in sequence but no parameters provided for it.");
            }

            const param_map* energy_lookup = nullptr;
            switch (lookup_type) {
                case ModLookup::Stacking:   energy_lookup = &param->stacking_energies(); break;
                case ModLookup::TerminalAU: energy_lookup = &param->terminal_energies(); break;
                case ModLookup::Mismatch:   energy_lookup = &param->mismatch_energies(); break;
                case ModLookup::Dangle5:    energy_lookup = &param->dangle5_energies(); break;
                case ModLookup::Dangle3:    energy_lookup = &param->dangle3_energies(); break;
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

        return NULL_ENERGY;
    }

    [[nodiscard]] static int get_mismatch_mod_energy(
        std::size_t i, std::size_t j, const std::vector<std::string_view>& mod_sequence,
        const all_mod_params& mp, bool is_closing = false) {
        if (mod_sequence.empty()) {
            return NULL_ENERGY;
        }
        std::string mismatch_key = get_mismatch_key(i, j, mod_sequence, is_closing);
        std::vector<std::string_view> unique_mod_bases = unique_mod_bases_in_string(mismatch_key);
        int mod_energy = get_mod_energy(mismatch_key, unique_mod_bases, mp, ModLookup::Mismatch);
        return mod_energy;
    }

    [[nodiscard]] static int get_dangle5_mod_energy(
        std::size_t i, std::size_t j, const std::vector<std::string_view>& mod_sequence,
        const all_mod_params& mp, bool is_closing = false) {
        if (mod_sequence.empty()) {
            return NULL_ENERGY;
        }
        std::string dangle5_key = get_dangle5_key(i, j, mod_sequence, is_closing);
        std::vector<std::string_view> unique_mod_bases = unique_mod_bases_in_string(dangle5_key);
        int mod_energy = get_mod_energy(dangle5_key, unique_mod_bases, mp, ModLookup::Dangle5);
        return mod_energy;
    }

    [[nodiscard]] static int get_dangle3_mod_energy(
        std::size_t i, std::size_t j, const std::vector<std::string_view>& mod_sequence,
        const all_mod_params& mp, bool is_closing = false) {
        if (mod_sequence.empty()) {
            return NULL_ENERGY;
        }
        std::string dangle3_key = get_dangle3_key(i, j, mod_sequence, is_closing);
        std::vector<std::string_view> unique_mod_bases = unique_mod_bases_in_string(dangle3_key);
        int mod_energy = get_mod_energy(dangle3_key, unique_mod_bases, mp, ModLookup::Dangle3);
        return mod_energy;
    }

    [[nodiscard]] static int get_terminalAU_mod_energy(
        std::size_t i, std::size_t j, const std::vector<std::string_view>& mod_sequence,
        const all_mod_params& mp, bool is_closing = false) {
        if (mod_sequence.empty()) {
            return NULL_ENERGY;
        }
        std::string terminal_key = get_terminal_key(i, j, mod_sequence, is_closing);
        std::vector<std::string_view> unique_mod_bases = unique_mod_bases_in_string(terminal_key);
        int mod_energy = get_mod_energy(terminal_key, unique_mod_bases, mp, ModLookup::TerminalAU);
        return mod_energy;
    }

    [[nodiscard]] static std::string get_mismatch_key(
        std::size_t i, std::size_t j, const std::vector<std::string_view>& mod_sequence,
        bool is_closing) {
        if (is_closing) {
            return ModBaseUtils::join_string_views({j, j - 1, i, i + 1}, mod_sequence);
        } else {
            return ModBaseUtils::join_string_views({i, i - 1, j, j + 1}, mod_sequence);
        }
    }

    [[nodiscard]] static std::string get_dangle3_key(
        std::size_t i, std::size_t j, const std::vector<std::string_view>& mod_sequence,
        bool is_closing) {
        if (is_closing) {
            return ModBaseUtils::join_string_views({j, i, i + 1}, mod_sequence);
        } else {
            return ModBaseUtils::join_string_views({i, j, j + 1}, mod_sequence);
        }
    }

    [[nodiscard]] static std::string get_dangle5_key(
        std::size_t i, std::size_t j, const std::vector<std::string_view>& mod_sequence,
        bool is_closing) {
        if (is_closing) {
            return ModBaseUtils::join_string_views({j, i, j - 1}, mod_sequence);
        } else {
            return ModBaseUtils::join_string_views({i, j, i - 1}, mod_sequence);
        }
    }

    [[nodiscard]] static std::string get_terminal_key(
        std::size_t i, std::size_t j, const std::vector<std::string_view>& mod_sequence,
        bool is_closing) {
        if (is_closing) {
            return ModBaseUtils::join_string_views({j, i}, mod_sequence);
        } else {
            return ModBaseUtils::join_string_views({i, j}, mod_sequence);
        }
    }

    [[nodiscard]] static ModDiffs get_mod_diffs(
        const LoopNode& node, int n5d, int n3d, unsigned int type,
        const std::vector<std::string_view>& unique_mod_bases,
        const std::vector<std::string_view>& mod_sequence, vrna_md_param& vp,
        const all_mod_params& mp, bool is_external, bool is_closing) {
        if (is_external && is_closing) {
            THROW_ERROR("An external loop cannot be a closing pair, check loop tree construction");
        }

        viennarna::vrna_param_t* P = vp.p;
        int mismatch = 0;
        int dangle5 = 0;
        int dangle3 = 0;
        int terminal = 0;
        std::string mismatch_key, dangle5_key, dangle3_key, terminal_key;
        // Unmodified energies for exterior stem (mismatch, dangle5, dangle3, terminalAU)
        if (n5d >= 0 && n3d >= 0) {
            mismatch = is_external ? P->mismatchExt[type][n5d][n3d] : P->mismatchM[type][n5d][n3d];
            mismatch_key = get_mismatch_key(node.begin, node.end, mod_sequence, is_closing);
        }

        if (n5d >= 0) {
            dangle5 = P->dangle5[type][n5d];
            dangle5_key = get_dangle5_key(node.begin, node.end, mod_sequence, is_closing);
        }

        if (n3d >= 0) {
            dangle3 = P->dangle3[type][n3d];
            dangle3_key = get_dangle3_key(node.begin, node.end, mod_sequence, is_closing);
        }

        if (type > 2) {
            terminal = P->TerminalAU;
            terminal_key = get_terminal_key(node.begin, node.end, mod_sequence, is_closing);
        }

        // Get modified energies
        int modMismatch =
            ModBaseUtils::get_mod_energy(mismatch_key, unique_mod_bases, mp, ModLookup::Mismatch);
        int modDangle5 =
            ModBaseUtils::get_mod_energy(dangle5_key, unique_mod_bases, mp, ModLookup::Dangle5);
        int modDangle3 =
            ModBaseUtils::get_mod_energy(dangle3_key, unique_mod_bases, mp, ModLookup::Dangle3);
        int modTerminalAU =
            ModBaseUtils::get_mod_energy(terminal_key, unique_mod_bases, mp, ModLookup::TerminalAU);

        modMismatch = modMismatch != NULL_ENERGY ? modMismatch : mismatch;
        modDangle5 = modDangle5 != NULL_ENERGY ? modDangle5 : dangle5;
        modDangle3 = modDangle3 != NULL_ENERGY ? modDangle3 : dangle3;
        modTerminalAU = modTerminalAU != NULL_ENERGY ? modTerminalAU : terminal;

        // Calculate differences
        if (vp.md.dangles == 0) {
            return ModDiffs(modTerminalAU - terminal, 0, 0, 0);
        }

        int diffMismatch = modMismatch - mismatch;
        int diff5 = modDangle5 - dangle5;
        int diff3 = modDangle3 - dangle3;
        int diffTerminal = modTerminalAU - terminal;

        return ModDiffs(diffTerminal, diffMismatch, diff5, diff3);
    }
    static void modify_dangle_set(DangleSet& dangle_set, ModDiffs diffs) {
        dangle_set.both_dangle += diffs.mismatch;
        dangle_set.left_dangle += diffs.n5d;
        dangle_set.right_dangle += diffs.n3d;
        dangle_set += diffs.terminalAU;  // adds terminalAU diff to all configurations
    }
};
}  // namespace knotergy
