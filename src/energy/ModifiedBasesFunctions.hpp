#pragma once

#include <utility>

#include "../preprocessing/RNAProcessor.hpp"
#include "./ViennaDangles.hpp"
#include "./ViennaFunctions.hpp"

namespace knotergy {

/**
 * @brief Enumeration of modified base energy lookup types.
 *
 * Specifies which type of energy parameter to look up for modified bases.
 */
enum class ModLookup { Stacking, Terminal, Mismatch, Dangle5, Dangle3 };

// Stores the differences in energy contributions due to modified bases
struct ModDiffs {
   ModDiffs(int terminal_diff, int mismatch_diff, int n5d_diff, int n3d_diff):
         terminalAU{terminal_diff},
         mismatch{mismatch_diff},
         n5d{n5d_diff},
         n3d{n3d_diff} {}
    const int terminalAU;
    const int mismatch;
    const int n5d;
    const int n3d;
};

class ModifiedBasesFunctions {
   public:
    /**
     * @brief Calculate stacking energy for base pairs with modified bases.
     *
     * @param i 5' position of outer base pair.
     * @param j 3' position of outer base pair.
     * @param ci 5' position of inner base pair.
     * @param cj 3' position of inner base pair.
     * @param sequence The unmodified RNA nucleotide sequence.
     * @param mod_sequence The modified RNA sequence (grapheme views).
     * @param mod_params Vector of modified base parameters.
     * @return Stacking energy in centicalories, accounting for modified bases.
     */
    static int find_mod_stack_energy(const size_t& i, const size_t& j, const size_t& ci,
                                        const size_t& cj, const std::string& sequence,
                                        const std::vector<std::string_view>& mod_sequence,
                                        vrna_md_param& vp,
                                        const std::vector<modified_base_param>& mod_params);
    
    /**
     * @brief Calculate external loop energy with modified bases.
     *
     * @param children Vector of child loop nodes in the external loop.
     * @param sequence The unmodified RNA nucleotide sequence.
     * @param mod_sequence The modified RNA sequence (grapheme views).
     * @param mod_params Vector of modified base parameters.
     * @return External loop energy in centicalories, accounting for modified bases.
     */
    static int find_mod_external_energy(const std::vector<std::shared_ptr<LoopNode>>& children, const std::string& sequence,
                                           const std::vector<std::string_view>& mod_sequence,
                                           vrna_md_param& vp,
                                           const std::vector<modified_base_param>& mod_params);

   private:
    /**
     * @brief Find unique modified bases at specified sequence positions.
     *
     * @param indices Vector of sequence indices to check.
     * @param mod_sequence The modified RNA sequence (grapheme views).
     * @return Vector of unique modified base string views found at those positions.
     */
    static std::vector<std::string_view> unique_modified_bases_at_indices(std::vector<size_t> indices, const std::vector<std::string_view>& mod_sequence);

    /**
     * @brief Join string views at specified indices into a single string.
     *
     * @param indices Vector of sequence indices.
     * @param mod_sequence The modified RNA sequence (grapheme views).
     * @return Concatenated string of bases at the specified indices.
     */
    static std::string join_string_views(std::vector<size_t> indices,
                                         const std::vector<std::string_view>& mod_sequence);

    /**
     * @brief Get modified energy or fall back to unmodified energy.
     *
     * Looks up the modified energy parameter for a given key and modified bases.
     * Returns the unmodified energy if no modified parameter is found.
     *
     * @param key The parameter key to look up.
     * @param modified Vector of modified base identifiers.
     * @param mod_params Vector of modified base parameters.
     * @param unmod_energy The unmodified energy to use as fallback.
     * @param lookup_type Type of energy lookup (Stacking, Terminal, etc.).
     * @return Energy value in centicalories.
     */
    static int get_mod_energy(const std::string& key,
                                 const std::vector<std::string_view>& modified,
                                 const std::vector<modified_base_param>& mod_params,
                                 int unmod_energy, ModLookup lookup_type);
    
    static void modify_dangle_set( DangleSet& original_set, ModDiffs diffs);
   
   // Returns the difference between modified and unmodified energy based on the lookup type and key
    static int get_mod_energy_difference(const std::string& key,
                                          const std::vector<std::string_view>& modified,
                                          const std::vector<modified_base_param>& mod_params,
                                          int unmod_energy, ModLookup lookup_type);
   
   static ModDiffs get_mod_dangle_energy_diffs(const std::shared_ptr<LoopNode>& c, const int n5d,
                                               const int n3d, const unsigned int type, [[maybe_unused]] const unsigned int r_type,
                                               const std::vector<std::string_view>& unique_mod_bases,
                                               const std::vector<std::string_view>& mod_sequence,
                                                vrna_md_param& vp,
                                               const std::vector<modified_base_param>& mod_params, bool is_external) ;
};

}  // namespace knotergy