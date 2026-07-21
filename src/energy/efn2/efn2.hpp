#pragma once

#include "io/parameters/ViennaParams.hpp"

#include <cmath>
#include <cstddef>
#include <string>

namespace knotergy {

// This is currently not used. But it fills the gap between ViennaRNA
// and RNAstructure for single-bulge corrections.
int efn2_single_bulge_correction(size_t i, size_t j, size_t ci, size_t cj,
                                 const std::string& sequence, const vrna_md_param& vp) {
    const size_t n1 = ci - i - 1;
    const size_t n2 = j - cj - 1;

    // Only applies to 1x0 or 0x1 bulges
    if (!((n1 == 1 && n2 == 0) || (n1 == 0 && n2 == 1))) {
        return 0;
    }

    int  count = 1;
    char bulged_base;

    if (n1 == 1) {
        bulged_base = sequence[i + 1];

        // Scan from outer pair toward 5'
        for (int k = static_cast<int>(i); k >= 0 && sequence[k] == bulged_base; --k) {
            ++count;
        }

        // Scan from inner pair toward 3'
        for (size_t k = ci; k < sequence.size() && sequence[k] == bulged_base; ++k) {
            ++count;
        }
    } else {
        bulged_base = sequence[j - 1];

        // Scan from inner pair toward 5'
        for (int k = static_cast<int>(cj); k >= 0 && sequence[k] == bulged_base; --k) {
            ++count;
        }

        // Scan from outer pair toward 3'
        for (size_t k = j; k < sequence.size() && sequence[k] == bulged_base; ++k) {
            ++count;
        }
    }

    // Vienna energies use 0.01 kcal/mol units.
    // RNAstructure rounds this correction to 0.1 kcal/mol.
    constexpr double R_KCAL             = 0.00198720425864083;
    const double     temperature_kelvin = vp.p->temperature + 273.15;

    int correction = -10 * static_cast<int>(lround(10.0 * R_KCAL * temperature_kelvin *
                                                   log(static_cast<double>(count))));

    // RNAstructure single-C movable-bulge bonus: -0.90 kcal/mol
    if (toupper(static_cast<unsigned char>(bulged_base)) == 'C' && count > 1) {
        correction -= 90;
    }

    return correction;
}
}  // namespace knotergy