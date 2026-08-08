#pragma once

// This file is separate from ViennaFunctions.hpp to avoid circular dependencies
// And helps seperate concerns of energy calculation and ViennaRNA encoding details

#include "io/parameters/ViennaParams.hpp"
#include "loop_tree/LoopNode.hpp"
#include "preprocessing/ProcessedRNAEntry.hpp"

#include <ViennaRNA/sequences/alphabet.hpp>

#include <tuple>

namespace viennarna = thermorna::viennarna;

namespace knotergy {
/**
 * @brief Utility functions for ViennaRNA nucleotide encoding and pair type operations.
 *
 * This class provides helper functions to convert nucleotides to ViennaRNA's internal
 * encoding format and determine base pair types. All functions work with ViennaRNA's
 * model parameters.
 */
class ViennaUtils {
   public:
    ViennaUtils() = default;
    ~ViennaUtils() = default;

    /**
     * @brief Encode two nucleotides for ViennaRNA functions.
     *
     * @param i First nucleotide character.
     * @param j Second nucleotide character.
     * @param md ViennaRNA model details for encoding.
     * @return Tuple of (encoded_i, encoded_j) in ViennaRNA format.
     */
    [[nodiscard]] static std::tuple<int, int> encode_nucleotides(const char& i, const char& j) {
        int encoded_i = fast_nucleotide_encode(i);
        int encoded_j = fast_nucleotide_encode(j);
        return std::make_tuple(encoded_i, encoded_j);
    }

    /**
     * @brief Encode outer dangling nucleotides adjacent to a base pair.
     *
     * @param i 5' position of the base pair.
     * @param j 3' position of the base pair.
     * @param sequence The RNA nucleotide sequence.
     * @param pair_table Base-pair indices for the RNA sequence (NULL_INDEX for unpaired).
     * @return Tuple of (encoded 5' dangle, encoded 3' dangle). Returns -1 if out of bounds.
     */
    [[nodiscard]] static std::tuple<int, int> encode_outer_dangles(
        const size_t i, const size_t j, const std::string& sequence,
        const std::vector<size_t>& pair_table, viennarna::vrna_md_t& md) {
        bool has_5d_dangle_out =
            i > 0 && (pair_table[i - 1] == NULL_INDEX || (md.dangles != 1 && md.dangles != 3));

        bool has_3d_dangle_out = j + 1 < sequence.size() && (pair_table[j + 1] == NULL_INDEX ||
                                                             (md.dangles != 1 && md.dangles != 3));

        int n5d_dangle = has_5d_dangle_out ? fast_nucleotide_encode(sequence[i - 1]) : -1;
        int n3d_dangle = has_3d_dangle_out ? fast_nucleotide_encode(sequence[j + 1]) : -1;
        return std::make_tuple(n5d_dangle, n3d_dangle);
    }

    /**
     * @brief Overload of encode_outer_dangles that takes a ProcessedRNAEntry for convenience.
     *
     * @param i 5' position of the base pair.
     * @param j 3' position of the base pair.
     * @param entry The ProcessedRNAEntry containing the sequence and pair_table.
     * @param md ViennaRNA model details for encoding.
     * @return Tuple of (encoded 5' dangle, encoded 3' dangle). Returns -1 if out of bounds.
     */
    [[nodiscard]] static std::tuple<int, int> encode_outer_dangles(const size_t i, const size_t j,
                                                                   const ProcessedRNAEntry& entry,
                                                                   viennarna::vrna_md_t& md) {
        return encode_outer_dangles(i, j, entry.get_sequence(), entry.get_pair_table(), md);
    }

    /**
     * @brief Encode inner nucleotides within a base pair.
     *
     * @param i 5' position of the base pair.
     * @param j 3' position of the base pair.
     * @param sequence The RNA nucleotide sequence.
     * @param pair_table Base-pair indices for the RNA sequence (NULL_INDEX for unpaired).
     * @param md ViennaRNA model details for encoding.
     * @return Tuple of (encoded nucleotide at i+1, encoded nucleotide at j-1).
     */
    [[nodiscard]] static std::tuple<int, int> encode_inner_dangles(
        const size_t i, const size_t j, const std::string& sequence,
        const std::vector<size_t>& pair_table, viennarna::vrna_md_t& md) {
        bool has_5d_dangle_in =
            (pair_table[i + 1] == NULL_INDEX || (md.dangles != 1 && md.dangles != 3));
        bool has_3d_dangle_in =
            (pair_table[j - 1] == NULL_INDEX || (md.dangles != 1 && md.dangles != 3));

        int encoded_i = has_5d_dangle_in ? fast_nucleotide_encode(sequence[i + 1]) : -1;
        int encoded_j = has_3d_dangle_in ? fast_nucleotide_encode(sequence[j - 1]) : -1;
        return std::make_tuple(encoded_i, encoded_j);
    }

    /**
     * @brief Overload of encode_inner_dangles that takes a ProcessedRNAEntry for convenience.
     *
     * @param i 5' position of the base pair.
     * @param j 3' position of the base pair.
     * @param entry The ProcessedRNAEntry containing the sequence and pair_table.
     * @param md ViennaRNA model details for encoding.
     * @return Tuple of (encoded nucleotide at i+1, encoded nucleotide at j-1).
     */
    [[nodiscard]] static std::tuple<int, int> encode_inner_dangles(const size_t i, const size_t j,
                                                                   const ProcessedRNAEntry& entry,
                                                                   viennarna::vrna_md_t& md) {
        return encode_inner_dangles(i, j, entry.get_sequence(), entry.get_pair_table(), md);
    }

    /**
     * @brief Get the ViennaRNA pair type for two nucleotides.
     *
     * @param i First nucleotide character.
     * @param j Second nucleotide character.
     * @return ViennaRNA pair type identifier (0 for non-canonical pairs).
     */
    [[nodiscard]] static unsigned int get_pair_type(const char& i, const char& j,
                                                    viennarna::vrna_md_t& md) {
        auto [encoded_i, encoded_j] = encode_nucleotides(i, j);
        return vrna_get_ptype_md(encoded_i, encoded_j, &md);
    }

    /**
     * @brief Get the reverse (complementary) pair type.
     *
     * @param type ViennaRNA pair type identifier.
     * @return Reversed pair type identifier.
     */
    [[nodiscard]] static unsigned int reverse_pair_type(unsigned int type,
                                                        viennarna::vrna_md_t& md) {
        return static_cast<unsigned int>(md.rtype[type]);
    }

    /**
     * @brief Get the reverse pair type for two nucleotides.
     *
     * @param i First nucleotide character.
     * @param j Second nucleotide character.
     * @return Reversed pair type identifier.
     */
    [[nodiscard]] static unsigned int reverse_pair_type(const char& i, const char& j,
                                                        viennarna::vrna_md_t& md) {
        return reverse_pair_type(get_pair_type(i, j, md), md);
    }

    /**
     * @brief A faster version of ViennaRNA's vrna_nucleotide_encode.
     *
     * @param c The nucleotide character to encode.s.
     * @return The encoded nucleotide.
     */
    [[nodiscard]] static constexpr int fast_nucleotide_encode(char c) {
        switch (c) {
            case 'A': return 1;
            case 'C': return 2;
            case 'G': return 3;
            case 'T':
            case 'U': return 4;
            case 'N': return 0;
            default:  return 0;
        }
    }
};
}  // namespace knotergy