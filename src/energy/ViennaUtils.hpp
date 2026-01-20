#pragma once

// This file is separate from ViennaFunctions.hpp to avoid circular dependencies (ViennaDangles.hpp needs ViennaUtils.hpp)

#include <tuple>
#include "../io/ViennaParams.hpp"

extern "C" {
#include <ViennaRNA/sequences/alphabet.h>
}

namespace knotergy {
    class ViennaUtils {
    public:
        ViennaUtils() = default;
        ~ViennaUtils() = default;

        static std::tuple<int, int> encode_nucleotides(const char& i, const char& j) {
            int encoded_i = vrna_nucleotide_encode(i, &ViennaParams::md);
            int encoded_j = vrna_nucleotide_encode(j, &ViennaParams::md);
            return std::make_tuple(encoded_i, encoded_j);
        }

        static unsigned int get_pair_type(const char& i, const char& j) {
            auto [encoded_i, encoded_j] = ViennaUtils::encode_nucleotides(i, j);
            return vrna_get_ptype_md(encoded_i, encoded_j, &ViennaParams::md);
        }

        static unsigned int reverse_pair_type(unsigned int type) {
            return static_cast<unsigned int>(ViennaParams::md.rtype[type]);
        }

        static unsigned int reverse_pair_type(const char& i, const char& j) {
            return ViennaUtils::reverse_pair_type(ViennaUtils::get_pair_type(i, j));
        }
    };
} 