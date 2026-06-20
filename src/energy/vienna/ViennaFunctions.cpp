
#include "ViennaFunctions.hpp"

namespace knotergy {

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
    return vp.p->stack[type][child_type];
}

int ViennaFunctions::stack_energy(const LoopNode& node, vrna_md_param& vp) {
    return vp.p->stack[node.pair_type][node.children[0]->r_pair_type];
}

// Helper function for stack_energy that takes BasePair objects instead of indices
int ViennaFunctions::stack_energy(BasePair pair, BasePair child, const std::string& sequence,
                                  vrna_md_param& vp) {
    return stack_energy(pair.i, pair.j, child.i, child.j, sequence, vp);
}

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
    int si1 = vrna_nucleotide_encode(sequence[i + 1], &vp.md);
    int sj1 = vrna_nucleotide_encode(sequence[j - 1], &vp.md);

    // If loop size < 7, you MUST pass the loop sequence substring
    // https://github.com/ViennaRNA/ViennaRNA/blob/219394580aec203a9d6f0d5450021e22642d5a83/src/ViennaRNA/eval/hairpin.h#L78C1-L81C94
    const char* loop_seq = nullptr;
    std::string loop_subseq;
    if (size < 7) {
        loop_subseq = sequence.substr(i, size + 2);  // length == size
        loop_seq = loop_subseq.c_str();              // temporary c-string
    }

    return vrna_E_hairpin(size, pair_type, si1, sj1, loop_seq, vp.p);
}

int ViennaFunctions::hairpin_energy(const LoopNode& node, const ProcessedRNAEntry& pRNA,
                                    bool& is_inf, vrna_md_param& vp) {
    const std::string& sequence = pRNA.get_sequence();
    size_t i = node.begin;
    size_t j = node.end;

    // loop size
    unsigned int size = static_cast<unsigned int>(node.end - node.begin - 1);

    // Max size for hairpin loop is 30 in ViennaRNA
    if (size <= 30 && vp.p->hairpin[size] == INF) {
        std::cerr
            << "Warning: Hairpin loop size is too small (usually < 3). Infinite Energy. Pairing: "
            << sequence[i] << "-" << sequence[j] << " i: " << i << ", j: " << j << std::endl;
        is_inf = true;
        return INF;
    }

    // If loop size < 7, you MUST pass the loop sequence substring
    // https://github.com/ViennaRNA/ViennaRNA/blob/219394580aec203a9d6f0d5450021e22642d5a83/src/ViennaRNA/eval/hairpin.h#L78C1-L81C94
    const char* loop_seq = nullptr;
    std::string loop_subseq;
    if (size < 7) {
        loop_subseq = sequence.substr(i, size + 2);  // length == size
        loop_seq = loop_subseq.c_str();              // temporary c-string
    }

    return vrna_E_hairpin(size, node.pair_type, node.n5d_inner, node.n3d_inner, loop_seq, vp.p);
}

int ViennaFunctions::hairpin_energy(const BasePair& pair, const std::string& sequence, bool& is_inf,
                                    vrna_md_param& vp) {
    return hairpin_energy(pair.i, pair.j, sequence, is_inf, vp);
}

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
    int si1 = vrna_nucleotide_encode(sequence[i + 1], &vp.md);   // 5' mismatch nt of closing pair
    int sj1 = vrna_nucleotide_encode(sequence[j - 1], &vp.md);   // 3' mismatch nt of closing pair
    int sp1 = vrna_nucleotide_encode(sequence[ci - 1], &vp.md);  // 5' mismatch nt of enclosed pair
    int sq1 = vrna_nucleotide_encode(sequence[cj + 1], &vp.md);  // 3' mismatch nt of enclosed pair
    return vrna_E_internal(n1, n2, type1, type2, si1, sj1, sp1, sq1, vp.p);
}

int ViennaFunctions::internal_loop_energy(const LoopNode& node, vrna_md_param& vp) {
    size_t i = node.begin;
    size_t j = node.end;
    LoopNode* child = node.children[0].get();
    size_t ci = child->begin;
    size_t cj = child->end;
    unsigned int type = node.pair_type;
    unsigned int rc_type = child->r_pair_type;  // reverse child pair type
    int si1 = node.n5d_inner;
    int sj1 = node.n3d_inner;
    int sp1 = child->n5d_outer;
    int sq1 = child->n3d_outer;

    unsigned int n1 = static_cast<unsigned int>(ci - i - 1);  // unpaired bases 5' side
    unsigned int n2 = static_cast<unsigned int>(j - cj - 1);  // unpaired bases 3' side
    return vrna_E_internal(n1, n2, type, rc_type, si1, sj1, sp1, sq1, vp.p);
}

int ViennaFunctions::internal_loop_energy(BasePair pair, BasePair child,
                                          const std::string& sequence, vrna_md_param& vp) {
    return internal_loop_energy(pair.i, pair.j, child.i, child.j, sequence, vp);
}

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

    // ------------------ Closing Pair Energy ------------------
    energy += vrna_E_multibranch_stem(node.r_pair_type, node.n3d_inner, node.n5d_inner, vp.p);

    // ------------------ Child Stems Energy ------------------
    for (const std::unique_ptr<LoopNode>& child : node.children) {
        if (child->loop_type == LoopType::Pseudoknot) continue;
        int n5d_outer = vp.md.dangles != 0 ? child->n5d_outer : -1;
        int n3d_outer = vp.md.dangles != 0 ? child->n3d_outer : -1;
        energy += vrna_E_multibranch_stem(child->pair_type, n5d_outer, n3d_outer, vp.p);
    }

    return energy;
}

int ViennaFunctions::external_energy(const std::vector<std::unique_ptr<LoopNode>>& children,
                                     const ProcessedRNAEntry& pRNA, vrna_md_param& vp) {
    // ------------------ Dangle 1 Energy ------------------
    if (vp.md.dangles == 1 || vp.md.dangles == 3) {
        return Dangle1::get_external_dangle_1(children, pRNA, vp);
    }

    // ------------------ No dangles or dangle type 2 ------------------
    int energy = 0;
    for (const std::unique_ptr<LoopNode>& child : children) {
        if (child->loop_type == LoopType::Pseudoknot) {
            continue;
        }
        int n5d_outer = vp.md.dangles != 0 ? child->n5d_outer : -1;
        int n3d_outer = vp.md.dangles != 0 ? child->n3d_outer : -1;
        energy += vrna_E_exterior_stem(child->pair_type, n5d_outer, n3d_outer, vp.p);
    }
    return energy;
}
}  // namespace knotergy