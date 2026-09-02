
#include "ViennaFunctions.hpp"

#include <ViennaRNA/eval/exterior.hpp>
#include <ViennaRNA/eval/hairpin.hpp>
#include <ViennaRNA/eval/internal.hpp>
#include <ViennaRNA/eval/multibranch.hpp>
#include <ViennaRNA/model.hpp>
#include <ViennaRNA/sequences/alphabet.hpp>
#include <ViennaRNA/utils/basic.hpp>

namespace viennarna = thermorna::viennarna;

namespace knotergy {

// --------------------------------------------
//               Stack Energy
// --------------------------------------------
int ViennaFunctions::stack_energy(size_t i, size_t j, size_t ci, size_t cj,
                                  const std::string& sequence, vrna_md_param& vp) {
    bool stacked = i + 1 == ci && j == cj + 1 && ci < cj && j < sequence.size();
    if (!stacked && vp.md.dangles != 3) {
        THROW_ERROR("Invalid indices for stack energy calculation. Received i: " +
                    std::to_string(i) + ", j: " + std::to_string(j) +
                    ", ci: " + std::to_string(ci) + ", cj: " + std::to_string(cj));
    }

    // c = child or nested base pair
    unsigned int type = ViennaUtils::get_pair_type(sequence[i], sequence[j], vp.md);
    unsigned int child_type = ViennaUtils::reverse_pair_type(sequence[ci], sequence[cj], vp.md);

    int salt_stack_correction = vp.p->SaltStack;
    return vp.p->stack[type][child_type] + salt_stack_correction;
}

int ViennaFunctions::stack_energy(const LoopNode& node, const std::string& sequence,
                                  vrna_md_param& vp) {
    return stack_energy(node.begin, node.end, node.children[0]->begin, node.children[0]->end,
                        sequence, vp);
}

// Helper function for stack_energy that takes PKBasePair objects instead of indices
int ViennaFunctions::stack_energy(PKBasePair pair, PKBasePair child, const std::string& sequence,
                                  vrna_md_param& vp) {
    return stack_energy(pair.i, pair.j, child.i, child.j, sequence, vp);
}

// --------------------------------------------
//             Hairpin Energy
// --------------------------------------------
int ViennaFunctions::hairpin_energy(size_t i, size_t j, const std::string& sequence, bool& is_inf,
                                    vrna_md_param& vp) {
    if (j <= i || j >= sequence.size()) {
        std::cerr << "Invalid indices for hairpin energy calculation." << std::endl;
        return 0;
    }

    // loop size
    unsigned int size = static_cast<unsigned int>(j - i - 1);

    // Max size for hairpin loop is 30 in ViennaRNA
    if (size <= 30 && vp.p->hairpin[size] == INF) {
        std::cerr
            << "Warning: Hairpin loop size is too small (usually < 3). Infinite Energy. Pairing: "
            << sequence[i] << "-" << sequence[j] << " i: " << i << ", j: " << j << std::endl;
        is_inf = true;
        return INF;
    }

    unsigned int pair_type = ViennaUtils::get_pair_type(sequence[i], sequence[j], vp.md);
    int si1 = ViennaUtils::fast_nucleotide_encode(sequence[i + 1]);
    int sj1 = ViennaUtils::fast_nucleotide_encode(sequence[j - 1]);

    // If loop size < 7, you MUST pass the loop sequence substring
    // https://github.com/ViennaRNA/ViennaRNA/blob/219394580aec203a9d6f0d5450021e22642d5a83/src/ViennaRNA/eval/hairpin.h#L78C1-L81C94
    const char* loop_seq = nullptr;
    std::string loop_subseq;
    if (size < 7) {
        loop_subseq = sequence.substr(i, size + 2);  // length == size
        loop_seq = loop_subseq.c_str();              // temporary c-string
    }

    return viennarna::vrna_E_hairpin(size, pair_type, si1, sj1, loop_seq, vp.p);
}

int ViennaFunctions::hairpin_energy(const LoopNode& node, const std::string& sequence, bool& is_inf,
                                    vrna_md_param& vp) {
    return hairpin_energy(node.begin, node.end, sequence, is_inf, vp);
}

int ViennaFunctions::hairpin_energy(const PKBasePair& pair, const std::string& sequence,
                                    bool& is_inf, vrna_md_param& vp) {
    return hairpin_energy(pair.i, pair.j, sequence, is_inf, vp);
}

// --------------------------------------------
//            Internal Loop Energy
// --------------------------------------------
int ViennaFunctions::internal_loop_energy(size_t i, size_t j, size_t ci, size_t cj,
                                          const std::string& sequence, vrna_md_param& vp) {
    // c = child or nested bp

    if (j <= i || cj <= ci || ci <= i || j <= cj || j >= sequence.size()) {
        THROW_ERROR("Invalid indices for internal loop energy calculation.");  // Add i, j, ci, cj
                                                                               // to message
    }

    unsigned int n1 = static_cast<unsigned int>(ci - i - 1);  // unpaired bases 5' side
    unsigned int n2 = static_cast<unsigned int>(j - cj - 1);  // unpaired bases 3' side
    unsigned int type1 = ViennaUtils::get_pair_type(sequence[i], sequence[j], vp.md);
    unsigned int type2 = ViennaUtils::reverse_pair_type(sequence[ci], sequence[cj], vp.md);
    int si1 = ViennaUtils::fast_nucleotide_encode(sequence[i + 1]);   // 5' nt of closing
    int sj1 = ViennaUtils::fast_nucleotide_encode(sequence[j - 1]);   // 3' nt of closing
    int sp1 = ViennaUtils::fast_nucleotide_encode(sequence[ci - 1]);  // 5' nt of child
    int sq1 = ViennaUtils::fast_nucleotide_encode(sequence[cj + 1]);  // 3' nt of child
    return viennarna::vrna_E_internal(n1, n2, type1, type2, si1, sj1, sp1, sq1, vp.p);
}

int ViennaFunctions::internal_loop_energy(const LoopNode& node, const std::string& sequence,
                                          vrna_md_param& vp) {
    return internal_loop_energy(node.begin, node.end, node.children[0]->begin,
                                node.children[0]->end, sequence, vp);
}

int ViennaFunctions::internal_loop_energy(PKBasePair pair, PKBasePair child,
                                          const std::string& sequence, vrna_md_param& vp) {
    return internal_loop_energy(pair.i, pair.j, child.i, child.j, sequence, vp);
}

// --------------------------------------------
//            Multibranch Energy
// --------------------------------------------

int ViennaFunctions::multibranch_energy(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                        vrna_md_param& vp) {
    // ------------------- Penalties -------------------
    // closing penalty + unpaired bases penalty
    int penalties = vp.p->MLclosing + node.exclusive_unpaired_bases_count * vp.p->MLbase;

    // ------------------ Dangle Energies ------------------
    if (vp.md.dangles == 1) {
        return Dangle1::get_multibranch_dangle_1(node, pRNA, vp) + penalties;
    }

    if (vp.md.dangles == 3) {
        return CoaxialStacking::get_multibranch_dangle_3(node, pRNA, vp) + penalties;
    }

    int energy = penalties;

    const std::string& sequence = pRNA.get_sequence();
    // ------------------ Closing Pair Energy ------------------
    unsigned int pair_type =
        ViennaUtils::reverse_pair_type(sequence[node.begin], sequence[node.end], vp.md);
    auto [n5d, n3d] = ViennaUtils::encode_inner_dangles(node.begin, node.end, pRNA, vp.md);

    energy += viennarna::vrna_E_multibranch_stem(pair_type, n3d, n5d, vp.p);

    // ------------------ Child Stems Energy ------------------
    for (const LoopNode* child : node.children) {
        if (child->loop_type == LoopType::Pseudoknot) continue;

        size_t ci = child->begin;
        size_t cj = child->end;
        unsigned int c_pair_type =
            ViennaUtils::get_pair_type(sequence[child->begin], sequence[child->end], vp.md);
        auto [n5d_outer, n3d_outer] = ViennaUtils::encode_outer_dangles(ci, cj, pRNA, vp.md);

        n5d_outer = vp.md.dangles != 0 ? n5d_outer : -1;
        n3d_outer = vp.md.dangles != 0 ? n3d_outer : -1;

        energy += viennarna::vrna_E_multibranch_stem(c_pair_type, n5d_outer, n3d_outer, vp.p);
    }

    return energy;
}

// --------------------------------------------
//            External Loop Energy
// --------------------------------------------

int ViennaFunctions::external_energy(const std::vector<LoopNode*>& children,
                                     const ProcessedRNAEntry& pRNA, vrna_md_param& vp) {
    // ------------------ Dangle 1 Energy ------------------
    if (vp.md.dangles == 1 || vp.md.dangles == 3) {
        return Dangle1::get_external_dangle_1(children, pRNA, vp);
    }

    const std::string& sequence = pRNA.get_sequence();

    // ------------------ No dangles or dangle type 2 ------------------
    int energy = 0;
    for (const LoopNode* child : children) {
        if (child->loop_type == LoopType::Pseudoknot) {
            continue;
        }
        auto [n5d_outer, n3d_outer] =
            ViennaUtils::encode_outer_dangles(child->begin, child->end, pRNA, vp.md);
        unsigned int pair_type =
            ViennaUtils::get_pair_type(sequence[child->begin], sequence[child->end], vp.md);
        n5d_outer = vp.md.dangles != 0 ? n5d_outer : -1;
        n3d_outer = vp.md.dangles != 0 ? n3d_outer : -1;
        energy += viennarna::vrna_E_exterior_stem(pair_type, n5d_outer, n3d_outer, vp.p);
    }
    return energy;
}
}  // namespace knotergy