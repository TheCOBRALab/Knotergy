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

/**
 * @brief Energy calculation functions for RNA structures with modified bases.
 *
 * This class provides methods to compute energy contributions when RNA sequences
 * contain modified nucleotides. It uses custom energy parameters for modified bases
 * while falling back to standard ViennaRNA calculations when appropriate.
 */
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
                                        const size_t& cj, std::string sequence,
                                        const std::vector<std::string_view>& mod_sequence,
                                        const std::vector<modified_base_params>& mod_params);
    
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
                                           const std::vector<std::string_view>& mod_sequence, const std::vector<modified_base_params>& mod_params);

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
                                 const std::vector<modified_base_params>& mod_params,
                                 int unmod_energy, ModLookup lookup_type);
    
    /**
     * @brief Get the energy difference between modified and unmodified parameters.
     *
     * @param key The parameter key to look up.
     * @param modified Vector of modified base identifiers.
     * @param mod_params Vector of modified base parameters.
     * @param unmod_energy The unmodified energy value.
     * @param lookup_type Type of energy lookup (Stacking, Terminal, etc.).
     * @return Energy difference in centicalories.
     */
    static int get_mod_energy_difference(const std::string& key,
                                       const std::vector<std::string_view>& modified,
                                       const std::vector<modified_base_params>& mod_params,
                                       int unmod_energy, ModLookup lookup_type);
   
    /**
     * @brief Modify a DangleSet with modified base energy parameters.
     *
     * Adjusts dangle energies in the set based on modified base parameters.
     *
     * @param original_set The DangleSet to modify (modified in place).
     * @param c Shared pointer to the child loop node.
     * @param type Pair type identifier.
     * @param n5d 5' dangle nucleotide encoding.
     * @param n3d 3' dangle nucleotide encoding.
     * @param unique_mod_bases Vector of unique modified bases.
     * @param mod_sequence The modified RNA sequence (grapheme views).
     * @param mod_params Vector of modified base parameters.
     * @param is_external Whether this is for an external loop.
     */
    static void modify_dangle_set( DangleSet& original_set,
                                        const std::shared_ptr<LoopNode>& c, unsigned int type,
                                        int n5d, int n3d, std::vector<std::string_view> unique_mod_bases,
                                        const std::vector<std::string_view>& mod_sequence,
                                        const std::vector<modified_base_params>& mod_params,
                                        bool is_external);
};

}  // namespace knotergy