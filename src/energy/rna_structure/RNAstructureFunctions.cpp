#include "energy/rna_structure/RNAstructureFunctions.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <optional>

namespace knotergy {

namespace {

// efn2's penalty(a, b), in RNAstructure's units (tenths of kcal/mol).
int raw_helix_end(std::size_t a, std::size_t b, const std::string& sequence,
                  const std::vector<std::string_view>& mod_sequence,
                  const rnastructure_param& rsp) {
    const rnastructure::Datatable& d = rsp.data;
    const auto code = [&](std::size_t k) {
        return RNAstructureUtils::rnastructure_code(k, sequence, mod_sequence, rsp);
    };
    if (b >= 1 && a + 1 < sequence.size()) {
        return d.helix_ends(code(a), code(b), code(a + 1), code(b - 1));
    }
    return d.helix_ends(code(b), code(a), code(b + 1), code(a - 1));  // at an end of the sequence
}

}  // namespace

// --------------------------------------------
//               Stack Energy
// --------------------------------------------
int RNAstructureFunctions::stack_energy(std::size_t i, std::size_t j, std::size_t ci,
                                        std::size_t cj, const std::string& sequence,
                                        const std::vector<std::string_view>& mod_sequence,
                                        const rnastructure_param& rsp) {
    const rnastructure::Datatable& d = rsp.data;
    rnastructure::Code l_code =
        RNAstructureUtils::rnastructure_code(i, sequence, mod_sequence, rsp);
    rnastructure::Code r_code =
        RNAstructureUtils::rnastructure_code(j, sequence, mod_sequence, rsp);
    rnastructure::Code cl_code =
        RNAstructureUtils::rnastructure_code(ci, sequence, mod_sequence, rsp);
    rnastructure::Code cr_code =
        RNAstructureUtils::rnastructure_code(cj, sequence, mod_sequence, rsp);

    return RNAstructureUtils::from_rnastructure(d.stack(l_code, r_code, cl_code, cr_code));
}

// --------------------------------------------
//               Helix End
// --------------------------------------------
int RNAstructureFunctions::helix_end(std::size_t a, std::size_t b, const std::string& sequence,
                                     const std::vector<std::string_view>& mod_sequence,
                                     const rnastructure_param& rsp) {
    return RNAstructureUtils::from_rnastructure(raw_helix_end(a, b, sequence, mod_sequence, rsp));
}

// --------------------------------------------
//             Hairpin Energy
// --------------------------------------------
int RNAstructureFunctions::hairpin_energy(std::size_t i, std::size_t j, const std::string& sequence,
                                          const std::vector<std::string_view>& mod_sequence,
                                          const rnastructure_param& rsp, bool& is_inf) {
    if (j <= i || j >= sequence.size()) {
        THROW_ERROR("Invalid indices for hairpin energy calculation: i = " + std::to_string(i) +
                    ", j = " + std::to_string(j));
    }

    const rnastructure::Datatable& d = rsp.data;
    const rnastructure::MiscLoop& m = d.misc();
    const auto code = [&](std::size_t k) {
        return RNAstructureUtils::rnastructure_code(k, sequence, mod_sequence, rsp);
    };

    // loop size
    const std::size_t size = j - i - 1;

    if (size < 3) {
        std::cerr
            << "Warning: Hairpin loop size is too small (usually < 3). Infinite Energy. Pairing: "
            << sequence[i] << "-" << sequence[j] << " i: " << i << ", j: " << j << std::endl;
        is_inf = true;
        return INF;
    }
    if (!d.can_pair(code(i), code(j))) {
        is_inf = true;
        return INF;
    }

    const int end = raw_helix_end(j, i, sequence, mod_sequence, rsp);

    // Special hairpins are looked up by the whole loop, closing pair included, and replace
    // every other term.
    if (size == 3 || size == 4 || size == 6) {
        std::array<rnastructure::Code, 8> loop{};
        for (std::size_t k = 0; k < size + 2; ++k) loop[k] = code(i + k);
        const std::optional<int> special = size == 3   ? d.triloop(loop.data())
                                           : size == 4 ? d.tloop(loop.data())
                                                       : d.hexaloop(loop.data());
        if (special) return RNAstructureUtils::from_rnastructure(*special + end);
    }

    // Loop initiation; beyond 30 nucleotides it is extrapolated logarithmically.
    int energy = m.eparam[4];  // always 0 in RNAstructure
    if (size <= 30) {
        energy += d.hairpin()[size];
    } else {
        energy +=
            d.hairpin()[30] +
            static_cast<int>(std::round(m.prelog * std::log(static_cast<double>(size) / 30.0)));
    }

    // First mismatch; loops of 3 have none.
    if (size > 3) energy += d.tstkh(code(i), code(j), code(i + 1), code(j - 1));

    // GU closure preceded by two Gs (5'-GGG...U-3')
    if (i >= 2 && d.has_symbol(code(i), "G") && d.has_symbol(code(j), "U") &&
        d.has_symbol(code(i - 1), "G") && d.has_symbol(code(i - 2), "G")) {
        energy += m.gubonus;
    }

    energy += end;

    // All-C loops
    const rnastructure::Code C = d.code("C");
    for (std::size_t k = i + 1; k < j; ++k) {
        if (code(k) != C) return RNAstructureUtils::from_rnastructure(energy);
    }
    energy += size == 3 ? m.c3 : m.cint + static_cast<int>(size) * m.cslope;
    return RNAstructureUtils::from_rnastructure(energy);
}

}  // namespace knotergy
