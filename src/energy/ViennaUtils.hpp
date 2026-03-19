#pragma once

// This file is separate from ViennaFunctions.hpp to avoid circular dependencies (ViennaDangles.hpp
// needs ViennaUtils.hpp)

#include <tuple>

#include "../io/ViennaParams.hpp"

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
     * @brief Encode two nucleotides for ViennaRNA functions.
     *
     * @param i First nucleotide character.
     * @param j Second nucleotide character.
     * @return Tuple of (encoded_i, encoded_j) in ViennaRNA format.
     */
    static std::tuple<int, int> encode_nucleotides(const char& i, const char& j, vrna_md_t& md) {
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
     * @return Tuple of (encoded 5' dangle, encoded 3' dangle). Returns -1 if out of bounds.
     */
    static std::tuple<int, int> encode_outer_dangles(const size_t i, const size_t j, const std::string& sequence, vrna_md_t& md) {
        int encoded_i = i > 0 ? vrna_nucleotide_encode(sequence[i - 1], &md) : -1;
        int encoded_j = j + 1 < sequence.size() ? vrna_nucleotide_encode(sequence[j + 1], &md) : -1;
        return std::make_tuple(encoded_i, encoded_j);
    }

    /**
     * @brief Encode inner nucleotides within a base pair.
     *
     * @param i 5' position of the base pair.
     * @param j 3' position of the base pair.
     * @param sequence The RNA nucleotide sequence.
     * @return Tuple of (encoded nucleotide at i+1, encoded nucleotide at j-1).
     */
    static std::tuple<int, int> encode_inner_dangles(const size_t i, const size_t j, const std::string& sequence, vrna_md_t& md) {
        int encoded_i = vrna_nucleotide_encode(sequence[i + 1], &md);
        int encoded_j = vrna_nucleotide_encode(sequence[j - 1], &md);
        return std::make_tuple(encoded_i, encoded_j);
    }

    /**
     * @brief Get the ViennaRNA pair type for two nucleotides.
     *
     * @param i First nucleotide character.
     * @param j Second nucleotide character.
     * @return ViennaRNA pair type identifier (0 for non-canonical pairs).
     */
    static unsigned int get_pair_type(const char& i, const char& j, vrna_md_t& md) {
        auto [encoded_i, encoded_j] = ViennaUtils::encode_nucleotides(i, j, md);
        return vrna_get_ptype_md(encoded_i, encoded_j, &md);
    }

    /**
     * @brief Get the reverse (complementary) pair type.
     *
     * @param type ViennaRNA pair type identifier.
     * @return Reversed pair type identifier.
     */
    static unsigned int reverse_pair_type(unsigned int type, vrna_md_t& md) {
        return static_cast<unsigned int>(md.rtype[type]);
    }

    /**
     * @brief Get the reverse pair type for two nucleotides.
     *
     * @param i First nucleotide character.
     * @param j Second nucleotide character.
     * @return Reversed pair type identifier.
     */
    static unsigned int reverse_pair_type(const char& i, const char& j, vrna_md_t& md) {
        return ViennaUtils::reverse_pair_type(ViennaUtils::get_pair_type(i, j, md), md);
    }
};
}  // namespace knotergy