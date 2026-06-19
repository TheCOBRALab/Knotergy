#pragma once

// This file is separate from ViennaFunctions.hpp to avoid circular dependencies (ViennaDangles.hpp
// needs ViennaUtils.hpp)

#include "io/parameters/ViennaParams.hpp"
#include "loop_tree/LoopNode.hpp"
#include "preprocessing/ProcessedRNAEntry.hpp"

#include <tuple>

extern "C" {
#include <ViennaRNA/sequences/alphabet.h>
}

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
     * @brief Populate ViennaRNA encoding fields for a LoopNode based on its position and the RNA
     * sequence.
     *
     * This function encodes the nucleotides at the loop boundaries and their adjacent positions
     * according to ViennaRNA's encoding scheme. It also determines the pair type and reverse pair
     * type for the loop's closing pair. The encoding of dangling nucleotides is determined based on
     * ViennaRNA's dangle settings and whether adjacent positions are paired or unpaired.
     *
     * @param node The LoopNode to populate with ViennaRNA encodings.
     * @param pRNA The ProcessedRNAEntry containing the RNA sequence and pairing information.
     * @param vp The ViennaRNA model parameters for encoding.
     */
    static void populate_node_encodings(LoopNode& node, const ProcessedRNAEntry& pRNA,
                                        vrna_md_param& vp) {
        if (node.loop_type == LoopType::External) return;  // Doesn't have pairings
        const std::string& sequence = pRNA.get_sequence();
        const std::vector<size_t>& pair_table = pRNA.get_pair_table();
        vrna_md_s& md = vp.md;

        size_t i = node.begin;
        size_t j = node.end;
        char ni = sequence[i];  // nucleotide i
        char nj = sequence[j];  // nucleotide j

        // Encode the pair's nucleotides
        std::tie(node.i_encoded, node.j_encoded) = encode_nucleotides(ni, nj, md);

        // Encode the pair type and reverse pair type
        node.pair_type = vrna_get_ptype_md(node.i_encoded, node.j_encoded, &md);
        node.r_pair_type = reverse_pair_type(node.pair_type, md);

        // NOTE: Not using this file's helper functions since encoding nucleotides multiple times
        // is surpiseingly expensive.

        // Encode dangling nucleotides based on ViennaRNA's dangle settings
        if (md.dangles != 0) {
            std::tie(node.n5d_outer, node.n3d_outer) =
                encode_outer_dangles(i, j, sequence, pair_table, md);

            // Stacked pairs don't have inner dangles even in dangle mode 2
            if (node.loop_type != LoopType::Stack) {
                std::tie(node.n5d_inner, node.n3d_inner) =
                    encode_inner_dangles(i, j, sequence, pair_table, md);
            }
        } else if (node.loop_type == LoopType::Hairpin) {
            // Hairpins always have inner dangles even in dangle mode 0
            std::tie(node.n5d_inner, node.n3d_inner) =
                encode_inner_dangles(i, j, sequence, pair_table, md);
        }
    }

    /**
     * @brief Encode two nucleotides for ViennaRNA functions.
     *
     * @param i First nucleotide character.
     * @param j Second nucleotide character.
     * @return Tuple of (encoded_i, encoded_j) in ViennaRNA format.
     */
    [[nodiscard]] static std::tuple<int, int> encode_nucleotides(const char& i, const char& j,
                                                                 vrna_md_t& md) {
        int encoded_i = vrna_nucleotide_encode(i, &md);
        int encoded_j = vrna_nucleotide_encode(j, &md);
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
        const std::vector<size_t>& pair_table, vrna_md_t& md) {
        if (md.dangles == 0) return std::make_tuple(-1, -1);

        bool has_5d_dangle_out =
            i > 0 && (pair_table[i - 1] == NULL_INDEX || (md.dangles != 1 && md.dangles != 3));

        bool has_3d_dangle_out = j + 1 < sequence.size() && (pair_table[j + 1] == NULL_INDEX ||
                                                             (md.dangles != 1 && md.dangles != 3));

        int n5d_dangle = has_5d_dangle_out ? vrna_nucleotide_encode(sequence[i - 1], &md) : -1;
        int n3d_dangle = has_3d_dangle_out ? vrna_nucleotide_encode(sequence[j + 1], &md) : -1;
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
                                                                   vrna_md_t& md) {
        return ViennaUtils::encode_outer_dangles(i, j, entry.get_sequence(), entry.get_pair_table(),
                                                 md);
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
        const std::vector<size_t>& pair_table, vrna_md_t& md) {
        bool has_5d_dangle_in =
            (pair_table[i + 1] == NULL_INDEX || (md.dangles != 1 && md.dangles != 3));
        bool has_3d_dangle_in =
            (pair_table[j - 1] == NULL_INDEX || (md.dangles != 1 && md.dangles != 3));

        int encoded_i = has_5d_dangle_in ? vrna_nucleotide_encode(sequence[i + 1], &md) : -1;
        int encoded_j = has_3d_dangle_in ? vrna_nucleotide_encode(sequence[j - 1], &md) : -1;
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
                                                                   vrna_md_t& md) {
        return ViennaUtils::encode_inner_dangles(i, j, entry.get_sequence(), entry.get_pair_table(),
                                                 md);
    }

    /**
     * @brief Get the ViennaRNA pair type for two nucleotides.
     *
     * @param i First nucleotide character.
     * @param j Second nucleotide character.
     * @return ViennaRNA pair type identifier (0 for non-canonical pairs).
     */
    [[nodiscard]] static unsigned int get_pair_type(const char& i, const char& j, vrna_md_t& md) {
        auto [encoded_i, encoded_j] = ViennaUtils::encode_nucleotides(i, j, md);
        return vrna_get_ptype_md(encoded_i, encoded_j, &md);
    }

    /**
     * @brief Get the reverse (complementary) pair type.
     *
     * @param type ViennaRNA pair type identifier.
     * @return Reversed pair type identifier.
     */
    [[nodiscard]] static unsigned int reverse_pair_type(unsigned int type, vrna_md_t& md) {
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
                                                        vrna_md_t& md) {
        return ViennaUtils::reverse_pair_type(ViennaUtils::get_pair_type(i, j, md), md);
    }
};
}  // namespace knotergy