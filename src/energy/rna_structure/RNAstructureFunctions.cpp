#include "energy/rna_structure/RNAstructureFunctions.hpp"

#include <RNAstructure/eval/exterior.hpp>
#include <RNAstructure/eval/hairpin.hpp>
#include <RNAstructure/eval/internal.hpp>
#include <RNAstructure/eval/loop_element.hpp>
#include <RNAstructure/eval/multibranch.hpp>

#include <iostream>

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

// The loop that runs from position `next` up to `end` (not included) as ThermoRNA takes it: its
// unpaired nucleotides and its branches, 5' to 3'. A pseudoknotted child is not a helix
// RNAstructure has stacking rules for: ThermoRNA gives it no helix end and lets nothing stack
// on it.
std::vector<rnastructure::LoopElement> loop_elements(const std::vector<LoopNode*>& children,
                                                     std::size_t next, std::size_t end,
                                                     const char* loop_name,
                                                     const ProcessedRNAEntry& pRNA,
                                                     const rnastructure_param& rsp) {
    const std::string& sequence = pRNA.get_sequence();
    const std::vector<std::string_view>& mod_sequence = pRNA.get_modified_sequence();
    const auto code = [&](std::size_t k) {
        return RNAstructureUtils::rnastructure_code(k, sequence, mod_sequence, rsp);
    };

    std::vector<rnastructure::LoopElement> elements;
    const auto add_unpaired_up_to = [&](std::size_t position) {
        for (; next < position; ++next) {
            elements.push_back(rnastructure::LoopElement::unpaired(code(next)));
        }
    };

    for (const LoopNode* child : children) {
        if (child->begin < next || child->end <= child->begin || child->end >= end) {
            THROW_ERROR(std::string("Invalid child of the ") + loop_name + ": begin = " +
                        std::to_string(child->begin) + ", end = " + std::to_string(child->end));
        }
        add_unpaired_up_to(child->begin);
        next = child->end + 1;

        // A helix is its first pair and the nucleotides after and before it, inside the helix.
        const std::size_t p = child->begin;
        const std::size_t q = child->end;
        elements.push_back(
            child->loop_type == LoopType::Pseudoknot
                ? rnastructure::LoopElement::other(code(p), code(q))
                : rnastructure::LoopElement::helix(code(p), code(q), code(p + 1), code(q - 1)));
    }
    add_unpaired_up_to(end);
    return elements;
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

    // loop size
    const std::size_t size = j - i - 1;

    if (size < 3) {
        std::cerr
            << "Warning: Hairpin loop size is too small (usually < 3). Infinite Energy. Pairing: "
            << sequence[i] << "-" << sequence[j] << " i: " << i << ", j: " << j << std::endl;
        is_inf = true;
        return INF;
    }

    const auto code = [&](std::size_t k) {
        return RNAstructureUtils::rnastructure_code(k, sequence, mod_sequence, rsp);
    };
    const rnastructure::Code none = rnastructure::kNoNucleotide;

    // The loop from i to j: special hairpins are looked up by the whole loop, and a loop of
    // nothing but C has a penalty.
    std::vector<rnastructure::Code> loop(size + 2);
    for (std::size_t k = 0; k < loop.size(); ++k) loop[k] = code(i + k);
    // The nucleotides outside the closing pair that efn2 looks at, where the sequence has
    // them: i - 1 and j + 1 for the helix end, i - 1 and i - 2 for the GGG-U closure bonus.
    const rnastructure::Code beyond_i = i >= 1 ? code(i - 1) : none;
    const rnastructure::Code beyond_i2 = i >= 2 ? code(i - 2) : none;
    const rnastructure::Code beyond_j = j + 1 < sequence.size() ? code(j + 1) : none;

    // A closing pair RNAstructure does not allow or a forbidden table entry (N beside a closing
    // helix of one pair, for example) makes the hairpin infinite.
    const int energy = RNAstructureUtils::from_rnastructure(
        rnastructure::hairpin_energy(size, loop.data(), beyond_i, beyond_i2, beyond_j, rsp.data));
    if (energy == INF) is_inf = true;
    return energy;
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
    if (n1 + n2 == 0) return stack_energy(i, j, ci, cj, sequence, mod_sequence, rsp);

    const auto code = [&](std::size_t k) {
        return RNAstructureUtils::rnastructure_code(k, sequence, mod_sequence, rsp);
    };
    const rnastructure::Code none = rnastructure::kNoNucleotide;

    // A single bulged nucleotide can sit at any of the identical nucleotides next to it, paired
    // or not: count the positions, 5' of the bulge and 3' of it.
    std::size_t bulge_run = 1;
    if (n1 + n2 == 1) {
        const std::size_t bulge = n1 == 1 ? i + 1 : j - 1;
        const rnastructure::Code bulged = code(bulge);
        for (std::size_t k = bulge; k > 0 && code(k - 1) == bulged; --k) ++bulge_run;
        for (std::size_t k = bulge + 1; k < sequence.size() && code(k) == bulged; ++k) ++bulge_run;
    }

    return RNAstructureUtils::from_rnastructure(rnastructure::internal_loop_energy(
        n1, n2, code(i), code(j), code(ci), code(cj),
        // The nucleotides next to the pairs, inside the loop
        code(i + 1), code(j - 1), code(ci - 1), code(cj + 1),
        // Outside the outer pair, where the sequence has them, and inside the inner pair: the
        // helix ends depend on them.
        i >= 1 ? code(i - 1) : none, j + 1 < sequence.size() ? code(j + 1) : none, code(ci + 1),
        code(cj - 1), bulge_run, rsp.data));
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
                                              const rnastructure_param& rsp,
                                              rnastructure::MultibranchRule rule,
                                              rnastructure::CoaxialStacking coaxial) {
    const std::string& sequence = pRNA.get_sequence();
    const std::vector<std::string_view>& mod_sequence = pRNA.get_modified_sequence();
    if (node.loop_type != LoopType::Multibranch || node.children.empty() ||
        node.end <= node.begin || node.end >= sequence.size()) {
        THROW_ERROR("multibranch_energy should only be called for multibranch loops: begin = " +
                    std::to_string(node.begin) + ", end = " + std::to_string(node.end));
    }

    const auto code = [&](std::size_t k) {
        return RNAstructureUtils::rnastructure_code(k, sequence, mod_sequence, rsp);
    };
    const rnastructure::Code none = rnastructure::kNoNucleotide;
    const std::size_t i = node.begin;
    const std::size_t j = node.end;

    // The closing pair, the nucleotides outside it where the sequence has them (its helix end
    // depends on them), and what is inside it.
    return RNAstructureUtils::from_rnastructure(rnastructure::multibranch_energy(
        code(i), code(j), i >= 1 ? code(i - 1) : none, j + 1 < sequence.size() ? code(j + 1) : none,
        loop_elements(node.children, i + 1, j, "multibranch loop", pRNA, rsp), rsp.data, rule,
        coaxial));
}

// --------------------------------------------
//            External Loop Energy
// --------------------------------------------
int RNAstructureFunctions::external_energy(const std::vector<LoopNode*>& children,
                                           const ProcessedRNAEntry& pRNA,
                                           const rnastructure_param& rsp,
                                           rnastructure::CoaxialStacking coaxial) {
    return RNAstructureUtils::from_rnastructure(rnastructure::exterior_energy(
        loop_elements(children, 0, pRNA.get_sequence().size(), "external loop", pRNA, rsp),
        rsp.data, coaxial));
}

}  // namespace knotergy
