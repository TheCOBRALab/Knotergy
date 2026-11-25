#include "energy/rna_structure/RNAstructureFunctions.hpp"

#include <algorithm>
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

// The last index of RNAstructure's dangle table: which side of the helix the nucleotide is on.
constexpr std::size_t kDangle3 = 1;
constexpr std::size_t kDangle5 = 2;

// efn2's erg4. The table has the pair as the loop sees it: first the nucleotide whose 3'
// neighbor is in the loop (the 3' dangle), then the one whose 5' neighbor is (the 5' dangle).
// For a branch of the loop that is j then i; for the pair that closes it, i then j.
int dangle_energy(std::size_t i, std::size_t j, std::size_t side, bool is_closing,
                  const std::string& sequence, const std::vector<std::string_view>& mod_sequence,
                  const rnastructure_param& rsp) {
    const bool has_nucleotide =
        is_closing ? j > i + 1 : (side == kDangle5 ? i > 0 : j + 1 < sequence.size());
    if (j <= i || j >= sequence.size() || !has_nucleotide) {
        THROW_ERROR(std::string("Invalid indices for ") + (side == kDangle5 ? "5'" : "3'") +
                    " dangle energy calculation: i = " + std::to_string(i) +
                    ", j = " + std::to_string(j) +
                    (is_closing ? ", dangle inside the pair" : ", dangle outside the pair"));
    }

    const rnastructure::Datatable& d = rsp.data;
    const auto code = [&](std::size_t k) {
        return RNAstructureUtils::rnastructure_code(k, sequence, mod_sequence, rsp);
    };
    const std::size_t first = is_closing ? i : j;
    const std::size_t second = is_closing ? j : i;
    const std::size_t dangling = side == kDangle3 ? first + 1 : second - 1;
    return RNAstructureUtils::from_rnastructure(
        d.dangle(code(first), code(second), code(dangling), side));
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

    // To Knotergy's units. A forbidden table entry (N beside a closing helix of one pair, for
    // example) makes the hairpin infinite as well.
    const auto finish = [&](int raw) {
        const int converted = RNAstructureUtils::from_rnastructure(raw);
        if (converted == INF) is_inf = true;
        return converted;
    };

    // Special hairpins are looked up by the whole loop, closing pair included, and replace
    // every other term.
    if (size == 3 || size == 4 || size == 6) {
        std::array<rnastructure::Code, 8> loop{};
        for (std::size_t k = 0; k < size + 2; ++k) loop[k] = code(i + k);
        const std::optional<int> special = size == 3   ? d.triloop(loop.data())
                                           : size == 4 ? d.tloop(loop.data())
                                                       : d.hexaloop(loop.data());
        if (special) return finish(*special + end);
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
        if (code(k) != C) return finish(energy);
    }
    energy += size == 3 ? m.c3 : m.cint + static_cast<int>(size) * m.cslope;
    return finish(energy);
}

// --------------------------------------------
//            Internal Loop Energy
// --------------------------------------------
int RNAstructureFunctions::internal_loop_energy(std::size_t i, std::size_t j, std::size_t ci,
                                                std::size_t cj, const std::string& sequence,
                                                const std::vector<std::string_view>& mod_sequence,
                                                const rnastructure_param& rsp) {
    // c = child or nested bp
    if (j <= i || cj <= ci || ci <= i || j <= cj || j >= sequence.size()) {
        THROW_ERROR("Invalid indices for internal loop energy calculation: i = " +
                    std::to_string(i) + ", j = " + std::to_string(j) +
                    ", ci = " + std::to_string(ci) + ", cj = " + std::to_string(cj));
    }

    const std::size_t n1 = ci - i - 1;  // unpaired bases 5' side
    const std::size_t n2 = j - cj - 1;  // unpaired bases 3' side
    const std::size_t size = n1 + n2;
    if (size == 0) return stack_energy(i, j, ci, cj, sequence, mod_sequence, rsp);

    const rnastructure::Datatable& d = rsp.data;
    const rnastructure::MiscLoop& m = d.misc();
    const auto code = [&](std::size_t k) {
        return RNAstructureUtils::rnastructure_code(k, sequence, mod_sequence, rsp);
    };
    if (!d.can_pair(code(i), code(j)) || !d.can_pair(code(ci), code(cj))) return INF;

    // Loop initiation by size; beyond 30 nucleotides it is extrapolated logarithmically.
    const auto initiation = [&](const std::vector<rnastructure::Energy>& by_size) {
        if (size <= 30) return static_cast<int>(by_size[size]);
        return by_size[30] +
               static_cast<int>(std::round(m.prelog * std::log(static_cast<double>(size) / 30.0)));
    };
    // The helix ends that face the loop: the outer pair from inside, the inner from outside.
    const auto helix_ends = [&] {
        return raw_helix_end(j, i, sequence, mod_sequence, rsp) +
               raw_helix_end(ci, cj, sequence, mod_sequence, rsp);
    };

    if (n1 == 0 || n2 == 0) {
        // ---------------- Bulge ----------------
        if (size > 1) {
            return RNAstructureUtils::from_rnastructure(initiation(d.bulge()) + m.eparam[2] +
                                                        helix_ends());
        }

        // A single bulged nucleotide does not interrupt the helix: the two pairs stack.
        int energy = d.stack(code(i), code(j), code(ci), code(cj)) + d.bulge()[1] + m.eparam[2];

        // The bulge can sit at any of the identical nucleotides next to it, paired or not;
        // count those positions, 5' of the bulge and 3' of it.
        const std::size_t before = n1 == 1 ? i : cj;
        const std::size_t after = n1 == 1 ? ci : j;
        const rnastructure::Code bulged = code(before + 1);
        int count = 1;
        for (std::size_t k = before; code(k) == bulged; --k) {
            ++count;
            if (k == 0) break;
        }
        for (std::size_t k = after; k < sequence.size() && code(k) == bulged; ++k) {
            ++count;
        }
        if (d.has_symbol(bulged, "C") && count > 1) energy += m.singlecbulge;
        energy -= static_cast<int>(std::round(d.RT() * 10 * std::log(static_cast<double>(count))));
        return RNAstructureUtils::from_rnastructure(energy);
    }

    // ---------------- Internal loop ----------------
    int energy = 0;
    if (n1 == 2 && n2 == 2) {
        energy = d.iloop22(code(i), code(ci), code(j), code(cj), code(i + 1), code(i + 2),
                           code(j - 1), code(j - 2));
    } else if (n1 == 1 && n2 == 2) {
        energy =
            d.iloop21(code(i), code(j), code(i + 1), code(j - 1), code(cj + 1), code(ci), code(cj));
    } else if (n1 == 2 && n2 == 1) {
        energy = d.iloop21(code(cj), code(ci), code(cj + 1), code(ci - 1), code(i + 1), code(j),
                           code(i));
    } else if (size == 2) {
        energy = d.iloop11(code(i), code(i + 1), code(ci), code(j), code(j - 1), code(cj));
    } else {
        // Mismatches next to both pairs; 1xn and 2x3 loops have tables of their own.
        rnastructure::Table4D mismatch = rnastructure::Table4D::Tstki;
        if (n1 == 1 || n2 == 1) {
            mismatch = rnastructure::Table4D::Tstki1n;
        } else if ((n1 == 2 && n2 == 3) || (n1 == 3 && n2 == 2)) {
            mismatch = rnastructure::Table4D::Tstki23;
        }
        energy = d.table4(mismatch, code(i), code(j), code(i + 1), code(j - 1)) +
                 d.table4(mismatch, code(cj), code(ci), code(cj + 1), code(ci - 1)) +
                 initiation(d.inter()) + m.eparam[3];

        // Asymmetry (Ninio) penalty
        const int lopsided = static_cast<int>(n1 > n2 ? n1 - n2 : n2 - n1);
        const std::size_t shorter = std::min<std::size_t>(2, std::min(n1, n2));
        energy += std::min<int>(m.maxpen, lopsided * m.poppen[shorter]);
    }
    return RNAstructureUtils::from_rnastructure(energy + helix_ends());
}

// --------------------------------------------
//              Dangle Energy
// --------------------------------------------
int RNAstructureFunctions::dangle5_energy(std::size_t i, std::size_t j, const std::string& sequence,
                                          const std::vector<std::string_view>& mod_sequence,
                                          const rnastructure_param& rsp, bool is_closing) {
    return dangle_energy(i, j, kDangle5, is_closing, sequence, mod_sequence, rsp);
}

int RNAstructureFunctions::dangle3_energy(std::size_t i, std::size_t j, const std::string& sequence,
                                          const std::vector<std::string_view>& mod_sequence,
                                          const rnastructure_param& rsp, bool is_closing) {
    return dangle_energy(i, j, kDangle3, is_closing, sequence, mod_sequence, rsp);
}

}  // namespace knotergy
