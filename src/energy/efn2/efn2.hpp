#pragma once

/**
 * RNAstructure's efn2 corrections for ViennaRNA.
 * This is still incomplete, but it is a start. The efn2 corrections are only applied to
 * single-bulge loops. Multiloop corrections are not yet implemented, but is planned for the future.
 *
 * Note: Knotergy should also support efn2's --simple mode.
 */

#include "io/parameters/ViennaParams.hpp"

#include <cmath>
#include <cstddef>
#include <string>

namespace knotergy {

// Fills the gap between ViennaRNA and RNAstructure for single-bulge corrections.
int efn2_single_bulge_correction(std::size_t i, std::size_t j, std::size_t ci, std::size_t cj,
                                 const std::string& sequence, const vrna_md_param& vp) {
    const std::size_t n1 = ci - i - 1;
    const std::size_t n2 = j - cj - 1;

    // Only applies to 1x0 or 0x1 bulges
    if (!((n1 == 1 && n2 == 0) || (n1 == 0 && n2 == 1))) {
        return 0;
    }

    int count = 1;
    char bulged_base;

    if (n1 == 1) {
        bulged_base = sequence[i + 1];

        // Scan from outer pair toward 5'
        for (std::size_t k = i; sequence[k] == bulged_base; --k) {
            ++count;
            if (k == 0) break;  // Prevent underflow
        }

        // Scan from inner pair toward 3'
        for (std::size_t k = ci; k < sequence.size() && sequence[k] == bulged_base; ++k) {
            ++count;
        }
    } else {
        bulged_base = sequence[j - 1];

        // Scan from inner pair toward 5'
        for (std::size_t k = cj; sequence[k] == bulged_base; --k) {
            ++count;
            if (k == 0) break;  // Prevent underflow
        }

        // Scan from outer pair toward 3'
        for (std::size_t k = j; k < sequence.size() && sequence[k] == bulged_base; ++k) {
            ++count;
        }
    }

    // Vienna energies use 0.01 kcal/mol units.
    // RNAstructure rounds this correction to 0.1 kcal/mol.
    constexpr double R_KCAL = 0.00198720425864083;
    const double temperature_kelvin = vp.p->temperature + 273.15;

    int correction = -10 * static_cast<int>(lround(10.0 * R_KCAL * temperature_kelvin *
                                                   log(static_cast<double>(count))));

    // RNAstructure single-C movable-bulge bonus: -0.90 kcal/mol
    if (toupper(static_cast<unsigned char>(bulged_base)) == 'C' && count > 1) {
        correction -= 90;
    }

    return correction;
}
}  // namespace knotergy