#include "energy/rna_structure/RNAstructureFunctions.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
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

// One step through an exterior or multibranch loop, 5' to 3': an unpaired nucleotide, a helix
// as the loop sees it (see dangle_energy), or a pseudoknot, which nothing stacks on.
struct LoopElement {
    enum class Type { Unpaired, Helix, Pseudoknot };
    Type type = Type::Unpaired;
    rnastructure::Code first = 0;   // the nucleotide, or where the loop leaves the helix
    rnastructure::Code second = 0;  // where the loop reaches the helix
};

// The minimisation of efn2's ergexterior and ergmulti, in RNAstructure's units: the lowest
// energy over the ways the helices and the nucleotides next to them can stack. best[n] is the
// lowest energy of the first n elements; each option uses up the elements it stacks.
int lowest_stacking_energy(const std::vector<LoopElement>& elements, rnastructure::Table4D mismatch,
                           const rnastructure::Datatable& d) {
    std::vector<int> best(elements.size() + 1, 0);
    for (std::size_t n = 1; n <= elements.size(); ++n) {
        const std::size_t k = n - 1;
        const LoopElement& e = elements[k];
        const auto helix = [&](std::size_t back) {
            return n > back && elements[k - back].type == LoopElement::Type::Helix;
        };
        const auto unpaired = [&](std::size_t back) {
            return n > back && elements[k - back].type == LoopElement::Type::Unpaired;
        };
        int lowest = best[n - 1];  // nothing stacks here
        const auto option = [&](std::size_t used, int energy) {
            lowest = std::min(lowest, best[n - used] + energy);
        };

        if (e.type == LoopElement::Type::Unpaired) {
            if (helix(1)) {
                const LoopElement& h = elements[k - 1];
                // 3' dangle on the helix before this nucleotide
                option(2, d.dangle(h.first, h.second, e.first, kDangle3));
                // Terminal mismatch with the nucleotide before that helix
                if (unpaired(2)) {
                    option(3,
                           d.table4(mismatch, h.first, h.second, e.first, elements[k - 2].first));
                }
                // The two helices before this nucleotide stack across the nucleotide between
                // them, which forms a mismatch with this one.
                if (unpaired(2) && helix(3)) {
                    const LoopElement& between = elements[k - 2];
                    const LoopElement& h2 = elements[k - 3];
                    option(4, d.coaxstack(h2.first, h2.second, between.first, e.first) +
                                  d.tstackcoax(h.first, h.second, e.first, between.first));
                }
            }
        } else if (e.type == LoopElement::Type::Helix) {
            // 5' dangle of the nucleotide before this helix
            if (unpaired(1)) {
                option(2, d.dangle(e.first, e.second, elements[k - 1].first, kDangle5));
            }
            // Flush coaxial stack on the helix right before this one
            if (helix(1)) {
                const LoopElement& h = elements[k - 1];
                option(2, d.coax(h.first, h.second, e.second, e.first));
            }
            // This helix and the one before it stack across the nucleotide between them, which
            // forms a mismatch with the nucleotide before that helix.
            if (unpaired(1) && helix(2) && unpaired(3)) {
                const LoopElement& between = elements[k - 1];
                const LoopElement& h = elements[k - 2];
                const LoopElement& before = elements[k - 3];
                option(4, d.coaxstack(e.first, e.second, before.first, between.first) +
                              d.tstackcoax(h.first, h.second, between.first, before.first));
            }
        }
        best[n] = lowest;
    }
    return best.back();
}

// Goes through a loop from position `next` up to `end` (not included) and adds what it finds
// to `elements`: unpaired nucleotides and the loop's branches. Returns the helix ends of the
// branches, in RNAstructure's units.
int add_loop_elements(std::vector<LoopElement>& elements, const std::vector<LoopNode*>& children,
                      std::size_t next, std::size_t end, const char* loop_name,
                      const ProcessedRNAEntry& pRNA, const rnastructure_param& rsp) {
    const std::string& sequence = pRNA.get_sequence();
    const std::vector<std::string_view>& mod_sequence = pRNA.get_modified_sequence();
    const auto code = [&](std::size_t k) {
        return RNAstructureUtils::rnastructure_code(k, sequence, mod_sequence, rsp);
    };
    const auto add_unpaired_up_to = [&](std::size_t position) {
        for (; next < position; ++next) {
            elements.push_back({LoopElement::Type::Unpaired, code(next), 0});
        }
    };

    int helix_ends = 0;
    for (const LoopNode* child : children) {
        if (child->begin < next || child->end <= child->begin || child->end >= end) {
            THROW_ERROR(std::string("Invalid child of the ") + loop_name + ": begin = " +
                        std::to_string(child->begin) + ", end = " + std::to_string(child->end));
        }
        add_unpaired_up_to(child->begin);
        next = child->end + 1;

        if (child->loop_type == LoopType::Pseudoknot) {
            elements.push_back({LoopElement::Type::Pseudoknot, 0, 0});
            continue;
        }
        elements.push_back({LoopElement::Type::Helix, code(child->end), code(child->begin)});
        helix_ends += raw_helix_end(child->begin, child->end, sequence, mod_sequence, rsp);
    }
    add_unpaired_up_to(end);
    return helix_ends;
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

// --------------------------------------------
//          Terminal Mismatch Energy
// --------------------------------------------
int RNAstructureFunctions::mismatch_energy(std::size_t i, std::size_t j,
                                           const std::string& sequence,
                                           const std::vector<std::string_view>& mod_sequence,
                                           const rnastructure_param& rsp, bool is_external,
                                           bool is_closing) {
    if (is_external && is_closing) {
        THROW_ERROR("An external loop has no closing pair to put a mismatch inside of.");
    }
    const bool has_nucleotides = is_closing ? j > i + 2 : (i > 0 && j + 1 < sequence.size());
    if (j <= i || j >= sequence.size() || !has_nucleotides) {
        THROW_ERROR("Invalid indices for mismatch energy calculation: i = " + std::to_string(i) +
                    ", j = " + std::to_string(j) +
                    (is_closing ? ", mismatch inside the pair" : ", mismatch outside the pair"));
    }

    const rnastructure::Datatable& d = rsp.data;
    const auto code = [&](std::size_t k) {
        return RNAstructureUtils::rnastructure_code(k, sequence, mod_sequence, rsp);
    };
    // The pair as the loop sees it, as for the dangles: then its 3' and its 5' dangle.
    const std::size_t first = is_closing ? i : j;
    const std::size_t second = is_closing ? j : i;
    const rnastructure::Table4D table =
        is_external ? rnastructure::Table4D::Tstack : rnastructure::Table4D::Tstkm;
    return RNAstructureUtils::from_rnastructure(
        d.table4(table, code(first), code(second), code(first + 1), code(second - 1)));
}

// --------------------------------------------
//            Multibranch Energy
// --------------------------------------------
int RNAstructureFunctions::multibranch_energy(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                              const rnastructure_param& rsp) {
    const std::string& sequence = pRNA.get_sequence();
    const std::vector<std::string_view>& mod_sequence = pRNA.get_modified_sequence();
    if (node.loop_type != LoopType::Multibranch || node.end <= node.begin ||
        node.end >= sequence.size()) {
        THROW_ERROR("multibranch_energy should only be called for multibranch loops: begin = " +
                    std::to_string(node.begin) + ", end = " + std::to_string(node.end));
    }

    const rnastructure::Datatable& d = rsp.data;
    const rnastructure::MiscLoop& m = d.misc();
    const auto code = [&](std::size_t k) {
        return RNAstructureUtils::rnastructure_code(k, sequence, mod_sequence, rsp);
    };

    // The loop is a circle. Start it at the closing pair, which the loop sees from inside.
    std::vector<LoopElement> elements;
    elements.push_back({LoopElement::Type::Helix, code(node.begin), code(node.end)});
    int energy = raw_helix_end(node.end, node.begin, sequence, mod_sequence, rsp) +
                 add_loop_elements(elements, node.children, node.begin + 1, node.end,
                                   "multibranch loop", pRNA, rsp);

    // Unpaired nucleotides before each helix, going round the loop; the closing pair is last.
    std::vector<int> gaps;
    int gap = 0;
    for (std::size_t k = 1; k < elements.size(); ++k) {
        if (elements[k].type == LoopElement::Type::Unpaired) {
            ++gap;
        } else {
            gaps.push_back(gap);
            gap = 0;
        }
    }
    gaps.push_back(gap);
    const int helices = static_cast<int>(gaps.size());
    const int unpaired = static_cast<int>(elements.size()) - helices;

    // Offset, penalty per helix and penalty per unpaired nucleotide. Beyond 8 nucleotides efn2
    // switches to a logarithm with a fixed coefficient of 1.1 kcal/mol.
    energy += m.efn2a + helices * m.efn2c;
    if (unpaired > 8) {
        energy += 8 * m.efn2b + static_cast<int>(11.0 * std::log(unpaired / 8.0) + 0.5);
    } else {
        energy += unpaired * m.efn2b;
    }

    // Strain: an odd number of helices with at most one unpaired nucleotide among them
    if (helices % 2 != 0 && unpaired <= 1) energy += m.strain;

    // Asymmetry: how much the gaps on the two sides of a helix differ, averaged over the
    // helices and capped at 2. In single precision like efn2, so that it rounds the same way.
    int differences = 0;
    for (std::size_t k = 0; k < gaps.size(); ++k) {
        differences += std::abs(gaps[k] - gaps[(k + 1) % gaps.size()]);
    }
    const float asymmetry =
        std::min(2.0F, static_cast<float>(differences) / static_cast<float>(helices));
    energy += static_cast<int>(static_cast<float>(m.mlasym) * asymmetry + 0.5);

    // Stacking, as in an exterior loop. The loop has to be opened somewhere to go through it,
    // and what stacks together covers up to 4 elements, so efn2 opens it in 4 places next to
    // each other and keeps the best: one of them does not cut through a stack.
    int stacking = 0;
    for (int opening = 0; opening < 4; ++opening) {
        stacking =
            std::min(stacking, lowest_stacking_energy(elements, rnastructure::Table4D::Tstkm, d));
        std::rotate(elements.begin(), elements.begin() + 1, elements.end());
    }

    return RNAstructureUtils::from_rnastructure(energy + stacking);
}

// --------------------------------------------
//            External Loop Energy
// --------------------------------------------
int RNAstructureFunctions::external_energy(const std::vector<LoopNode*>& children,
                                           const ProcessedRNAEntry& pRNA,
                                           const rnastructure_param& rsp) {
    std::vector<LoopElement> elements;
    const int helix_ends = add_loop_elements(elements, children, 0, pRNA.get_sequence().size(),
                                             "external loop", pRNA, rsp);
    return RNAstructureUtils::from_rnastructure(
        helix_ends + lowest_stacking_energy(elements, rnastructure::Table4D::Tstack, rsp.data));
}

}  // namespace knotergy
