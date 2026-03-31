
#include "ViennaFunctions.hpp"

namespace knotergy {

int ViennaFunctions::stack_energy(size_t i, size_t j, size_t ci, size_t cj,
                                  const std::string& sequence, vrna_md_param& vp) {
    if (j <= i || cj <= ci || ci <= i || j <= cj || j >= sequence.size()) {
        std::cerr << "Invalid indices for stack energy calculation." << std::endl;
        return 0;
    }
    // c = child or nested base pair
    unsigned int type1 = ViennaUtils::get_pair_type(sequence[i], sequence[j], vp.md);
    unsigned int type2 = ViennaUtils::reverse_pair_type(sequence[ci], sequence[cj], vp.md);
    return vp.p->stack[type1][type2];
}

// Helper function for stack_energy that takes BasePair objects instead of indices
int ViennaFunctions::stack_energy(BasePair pair, BasePair child, const std::string& sequence, 
                                  vrna_md_param& vp) {
    return stack_energy(pair.i, pair.j, child.i, child.j, sequence, vp);
}

int ViennaFunctions::hairpin_energy(size_t i, size_t j, const std::string& sequence, 
                                    bool& is_inf, vrna_md_param& vp) {
    if (j <= i || j >= sequence.size()) {
        std::cerr << "Invalid indices for hairpin energy calculation." << std::endl;
        return 0;
    }

    // loop size
    unsigned int size = static_cast<unsigned int>(j - i - 1);
    if (vp.p->hairpin[size] == INF) {
        std::cerr << "Warning: Hairpin loop size is too small (usually < 3). Infinite Energy. Sequence: "
                  << sequence << " i: " << i << ", j: " << j << std::endl;
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

int ViennaFunctions::hairpin_energy(const BasePair& pair, const std::string& sequence, 
                                    bool& is_inf, vrna_md_param& vp) {
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
    int si1 = vrna_nucleotide_encode(sequence[i + 1], &vp.md);  // 5' mismatch nt of closing pair
    int sj1 = vrna_nucleotide_encode(sequence[j - 1], &vp.md);  // 3' mismatch nt of closing pair
    int sp1 = vrna_nucleotide_encode(sequence[ci - 1], &vp.md);  // 5' mismatch nt of enclosed pair
    int sq1 = vrna_nucleotide_encode(sequence[cj + 1], &vp.md);  // 3' mismatch nt of enclosed pair
    return vrna_E_internal(n1, n2, type1, type2, si1, sj1, sp1, sq1, vp.p);
}

int ViennaFunctions::internal_loop_energy(BasePair pair, BasePair child,
                                          const std::string& sequence, vrna_md_param& vp) {
    return internal_loop_energy(pair.i, pair.j, child.i, child.j, sequence, vp);
}

int ViennaFunctions::multibranch_energy(const LoopNode& node, const ProcessedRNAEntry& pRNA, vrna_md_param& vp) {
    const std::string& sequence = pRNA.get_sequence();

    // ------------------- Penalties -------------------
    int energy = vp.p->MLclosing;  // closing penalty
    energy += node.exclusive_unpaired_bases_count * vp.p->MLbase;  // unpaired bases penalty

    // ------------------ Dangle 1 Energy ------------------
    if (vp.md.dangles == 1) {
        return ViennaDangles::get_multibranch_dangle_1(node, pRNA, vp) + energy;
    }

    size_t i = node.begin;
    size_t j = node.end;

    // ------------------ Closing Pair Energy ------------------
    unsigned int pair_type = ViennaUtils::reverse_pair_type(sequence[i], sequence[j], vp.md);
    auto [n5d, n3d] = ViennaUtils::encode_inner_dangles(i, j, pRNA, vp.md);
    if (vp.md.dangles == 0) {
        n5d = -1;
        n3d = -1;
    }

    energy += vrna_E_multibranch_stem(pair_type, n3d, n5d, vp.p);

    // ------------------ Child Stems Energy ------------------
    for (const std::shared_ptr<LoopNode>& child : node.children) {
        size_t ci = child->begin;
        size_t cj = child->end;
        unsigned int c_pair_type = ViennaUtils::get_pair_type(sequence[ci], sequence[cj], vp.md);
        auto [child_n5d, child_n3d] = ViennaUtils::encode_outer_dangles(ci, cj, pRNA, vp.md);

        if (vp.md.dangles == 0) {
            child_n5d = -1;
            child_n3d = -1;
        }

        energy += vrna_E_multibranch_stem(c_pair_type, child_n5d, child_n3d, vp.p);
    }

    return energy;
}

int ViennaFunctions::external_energy(const std::vector<std::shared_ptr<LoopNode>>& children,
                                     const ProcessedRNAEntry& pRNA, vrna_md_param& vp) {
    const std::string& sequence = pRNA.get_sequence();

    // ------------------ Dangle 1 Energy ------------------
    if (vp.md.dangles == 1) {
        return ViennaDangles::get_external_dangle_1(children, pRNA, vp);
    }

    // ------------------ No dangles or dangle type 2 ------------------
    int energy = 0;
    for (std::shared_ptr<LoopNode> c : children) {
        if (c->loop_type != LoopType::Pseudoknot) {
            unsigned int pair_type = 
                           ViennaUtils::get_pair_type(sequence[c->begin], sequence[c->end], vp.md);

            // Check for dangling ends at sequence boundaries (-1 indicates no dangle)
            auto [n5d, n3d] = ViennaUtils::encode_outer_dangles(c->begin, c->end, pRNA, vp.md);

            if (vp.md.dangles == 0) {
                n5d = -1;
                n3d = -1;
            }
            energy += vrna_E_exterior_stem(pair_type, n5d, n3d, vp.p);
        }
    }
    return energy;
}
}  // namespace knotergy