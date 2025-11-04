#pragma once

extern "C" {
#include <ViennaRNA/sequences/alphabet.h>
}

namespace knotergy {
    enum DangleIdx { None = 0, Left = 1, Right = 2, Both = 3 };
    class ViennaUtils {
    public:
        ViennaUtils() = default;
        ~ViennaUtils() = default;
        static unsigned int get_pair_type(const char& i, const char& j, vrna_md_t& md) {
            int encoded_i = vrna_nucleotide_encode(i, &md);
            int encoded_j = vrna_nucleotide_encode(j, &md);
            return vrna_get_ptype_md(encoded_i, encoded_j, &md);
        }

        static unsigned int reverse_pair_type(unsigned int type, vrna_md_t& md) {
            return static_cast<unsigned int>(md.rtype[type]);
        }
    };
} 