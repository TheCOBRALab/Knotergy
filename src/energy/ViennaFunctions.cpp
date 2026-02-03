
#include "ViennaFunctions.hpp"

namespace knotergy {

int ViennaFunctions::stack_energy(size_t i, size_t j, size_t ci, size_t cj,
                                  const std::string& sequence) {
    if (j <= i || cj <= ci || ci <= i || j <= cj || j >= sequence.size()) {
        std::cerr << "Invalid indices for stack energy calculation." << std::endl;
        return 0;
    }
    // c = child or nested base pair
    unsigned int type1 = ViennaUtils::get_pair_type(sequence[i], sequence[j]);
    unsigned int type2 = ViennaUtils::reverse_pair_type(sequence[ci], sequence[cj]);
    return ViennaParams::p->stack[type1][type2];
}

int ViennaFunctions::stack_energy(BasePair pair, BasePair child, const std::string& sequence) {
    return stack_energy(pair.i, pair.j, child.i, child.j, sequence);
}

int ViennaFunctions::hairpin_energy(size_t i, size_t j, const std::string& sequence) {
    if (j <= i || j >= sequence.size()) {
        std::cerr << "Invalid indices for hairpin energy calculation." << std::endl;
        return 0;
    }

    // loop size
    unsigned int size = static_cast<unsigned int>(j - i - 1);
    if (size < 3) {
        std::cerr << "Warning: Hairpin loop size is less than 3. Infinite Energy. Sequence: "
                  << sequence << " i: " << i << ", j: " << j << std::endl;
        return INF;
    }

    unsigned int pair_type = ViennaUtils::get_pair_type(sequence[i], sequence[j]);
    int si1 = vrna_nucleotide_encode(sequence[i + 1], &ViennaParams::md);
    int sj1 = vrna_nucleotide_encode(sequence[j - 1], &ViennaParams::md);

    // If loop size < 7, you MUST pass the loop sequence substring
    // https://github.com/ViennaRNA/ViennaRNA/blob/219394580aec203a9d6f0d5450021e22642d5a83/src/ViennaRNA/eval/hairpin.h#L78C1-L81C94
    const char* loop_seq = nullptr;
    std::string loop_subseq;
    if (size < 7) {
        loop_subseq = sequence.substr(i, size + 2);  // length == size
        loop_seq = loop_subseq.c_str();              // temporary c-string
    }

    return vrna_E_hairpin(size, pair_type, si1, sj1, loop_seq, ViennaParams::p);
}

int ViennaFunctions::hairpin_energy(const BasePair& pair, const std::string& sequence) {
    return hairpin_energy(pair.i, pair.j, sequence);
}

int ViennaFunctions::internal_loop_energy(size_t i, size_t j, size_t ci, size_t cj,
                                          const std::string& sequence) {
    // c = child or nested bp

    if (j <= i || cj <= ci || ci <= i || j <= cj || j >= sequence.size()) {
        THROW_ERROR("Invalid indices for internal loop energy calculation.");  // Add i, j, ci, cj
                                                                               // to message
    }

    unsigned int n1 = static_cast<unsigned int>(ci - i - 1);  // unpaired bases 5' side
    unsigned int n2 = static_cast<unsigned int>(j - cj - 1);  // unpaired bases 3' side
    unsigned int type1 = ViennaUtils::get_pair_type(sequence[i], sequence[j]);
    unsigned int type2 = ViennaUtils::reverse_pair_type(sequence[ci], sequence[cj]);
    int si1 = vrna_nucleotide_encode(sequence[i + 1],
                                     &ViennaParams::md);  // 5' mismatch nt of closing pair
    int sj1 = vrna_nucleotide_encode(sequence[j - 1],
                                     &ViennaParams::md);  // 3' mismatch nt of closing pair
    int sp1 = vrna_nucleotide_encode(sequence[ci - 1],
                                     &ViennaParams::md);  // 5' mismatch nt of enclosed pair
    int sq1 = vrna_nucleotide_encode(sequence[cj + 1],
                                     &ViennaParams::md);  // 3' mismatch nt of enclosed pair
    return vrna_E_internal(n1, n2, type1, type2, si1, sj1, sp1, sq1, ViennaParams::p);
}

int ViennaFunctions::internal_loop_energy(BasePair pair, BasePair child,
                                          const std::string& sequence) {
    return internal_loop_energy(pair.i, pair.j, child.i, child.j, sequence);
}

int ViennaFunctions::multibranch_energy(const LoopNode& node, const std::string& sequence) {
    // ------------------- Penalties -------------------
    int energy = ViennaParams::p->MLclosing;  // closing penalty
    energy +=
        node.exclusive_unpaired_bases_count * ViennaParams::p->MLbase;  // unpaired bases penalty

    // ------------------ Dangle 1 Energy ------------------
    if (ViennaParams::md.dangles == 1) {
        return ViennaDangles::get_multibranch_dangle_1(node, sequence) + energy;
    }

    size_t i = node.begin;
    size_t j = node.end;

    // ------------------ Closing Pair Energy ------------------
    unsigned int pair_type = ViennaUtils::reverse_pair_type(sequence[i], sequence[j]);
    auto [n5d, n3d] = ViennaUtils::encode_inner_dangles(i, j, sequence);
    if (ViennaParams::md.dangles == 0) {
        n5d = -1;
        n3d = -1;
    }

    energy += vrna_E_multibranch_stem(pair_type, n3d, n5d, ViennaParams::p);

    // ------------------ Child Stems Energy ------------------
    for (const std::shared_ptr<LoopNode>& child : node.children) {
        // if (child->loop_type == LoopType::Pseudoknot) {
        //     energy += child->number_of_bands * P->MLintern[child_pair_type];
        //     continue;
        // }

        size_t ci = child->begin;
        size_t cj = child->end;
        unsigned int child_pair_type = ViennaUtils::get_pair_type(sequence[ci], sequence[cj]);
        auto [child_n5d, child_n3d] =
            ViennaUtils::encode_outer_dangles(ci, cj, sequence);

        if (ViennaParams::md.dangles == 0) {
            child_n5d = -1;
            child_n3d = -1;
        }

        energy += vrna_E_multibranch_stem(child_pair_type, child_n5d, child_n3d, ViennaParams::p);
    }

    return energy;
}

int ViennaFunctions::external_energy(const std::vector<std::shared_ptr<LoopNode>>& children,
                                     const std::string& sequence) {

    // ------------------ Dangle 1 Energy ------------------
    if (ViennaParams::md.dangles == 1) {
        return ViennaDangles::get_external_dangle_1(children, sequence);
    }

    // ------------------ No dangles or dangle type 2 ------------------
    int energy = 0;
    for (std::shared_ptr<LoopNode> c : children) {
        if (c->loop_type != LoopType::Pseudoknot) {
            unsigned int pair_type = ViennaUtils::get_pair_type(sequence[c->begin], sequence[c->end]);

            // Check for dangling ends at sequence boundaries (-1 indicates no dangle)
            auto [n5d, n3d] = ViennaUtils::encode_outer_dangles(c->begin, c->end, sequence);

            if (ViennaParams::md.dangles == 0) {
                n5d = -1;
                n3d = -1;
            }
            energy += vrna_E_exterior_stem(pair_type, n5d, n3d, ViennaParams::p);
        }
    }
    return energy;
}
}  // namespace knotergy