
#include "ViennaFunctions.hpp"

namespace knotergy {

const vrna_param_t* ViennaFunctions::get_parameters() const {return P;};

int ViennaFunctions::stack_energy(size_t i, size_t j, size_t ci, size_t cj, const std::string& sequence){
    if (j <= i || cj <= ci || ci <= i || j <= cj || j >= sequence.size()) {
        std::cerr << "Invalid indices for stack energy calculation." << std::endl;
        return 0;
    }
    // c = child or nested base pair
    unsigned int type1 = ViennaUtils::get_pair_type(sequence[i], sequence[j], md);
    unsigned int type2 = ViennaUtils::get_pair_type(sequence[ci], sequence[cj], md);
    return P->stack[type1][ViennaUtils::reverse_pair_type(type2, md)];
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
    }

    unsigned int pair_type = ViennaUtils::get_pair_type(sequence[i], sequence[j], md);
    int si1 = vrna_nucleotide_encode(sequence[i + 1], &md);
    int sj1 = vrna_nucleotide_encode(sequence[j - 1], &md);

    // If loop size < 7, you MUST pass the loop sequence substring
    // https://github.com/ViennaRNA/ViennaRNA/blob/219394580aec203a9d6f0d5450021e22642d5a83/src/ViennaRNA/eval/hairpin.h#L78C1-L81C94
    const char* loop_seq = nullptr;
    std::string loop_subseq;
    if (size < 7) {
        loop_subseq = sequence.substr(i, size + 2);  // length == size
        loop_seq = loop_subseq.c_str();              // temporary c-string
    }

    return vrna_E_hairpin(size, pair_type, si1, sj1, loop_seq, P);
}

int ViennaFunctions::hairpin_energy(const BasePair& pair, const std::string& sequence) {
    return hairpin_energy(pair.i, pair.j, sequence);
}

int ViennaFunctions::internal_loop_energy(size_t i, size_t j, size_t ci, size_t cj, const std::string& sequence) {
    // c = child or nested bp

    if (j <= i || cj <= ci || ci <= i || j <= cj || j >= sequence.size()) {
        THROW_ERROR("Invalid indices for internal loop energy calculation."); // Add i, j, ci, cj to message
    }

    unsigned int n1 = static_cast<unsigned int>(ci - i - 1);
    unsigned int n2 = static_cast<unsigned int>(j - cj - 1);
    unsigned int pair_type = ViennaUtils::get_pair_type(sequence[i], sequence[j], md);
    unsigned int pair_type2 = ViennaUtils::reverse_pair_type(ViennaUtils::get_pair_type(sequence[ci], sequence[cj], md), md);
    int si1 = vrna_nucleotide_encode(sequence[i + 1], &md);   // 5' mismatch of closing pair
    int sj1 = vrna_nucleotide_encode(sequence[j - 1], &md);   // 3' mismatch of closing pair
    int sp1 = vrna_nucleotide_encode(sequence[ci - 1], &md);  // 5' mismatch of enclosed pair
    int sq1 = vrna_nucleotide_encode(sequence[cj + 1], &md);  // 3' mismatch of enclosed pair

    return vrna_E_internal(n1, n2, pair_type, pair_type2, si1, sj1, sp1, sq1, P);
}

int ViennaFunctions::internal_loop_energy(BasePair pair, BasePair child, const std::string& sequence) {
    return internal_loop_energy(pair.i, pair.j, child.i, child.j, sequence);
}

int ViennaFunctions::multibranch_energy(const LoopNode& node, const std::string& sequence) {
    size_t i = node.begin;
    size_t j = node.end;

    // penalties
    int energy = P->MLclosing; // closing penalty
    energy += node.exclusive_unpaired_bases_count * P->MLbase; // unpaired bases penalty

    int n5d, n3d;
    switch (md.dangles){
        case 0:
            n5d = -1;
            n3d = -1;
            break;
        case 1:
            return ViennaDangles::get_multi_dangle_1(node, sequence, md) + energy;
        case 2:
            n5d = vrna_nucleotide_encode(sequence[i + 1], &md);
            n3d = vrna_nucleotide_encode(sequence[j - 1], &md);
            break;
        default:
            THROW_ERROR("Dangle model `" + std::to_string(md.dangles) + "` not supported in multi-branch energy calculation.");
    }

    unsigned int pair_type = ViennaUtils::get_pair_type(sequence[i], sequence[j], md);

    // energy of the closing pair
    energy += vrna_E_multibranch_stem(ViennaUtils::reverse_pair_type(pair_type, md), n3d, n5d, P);

    const std::vector<std::shared_ptr<LoopNode>>& children = node.children;
    for (const std::shared_ptr<LoopNode>& child : children) {
        if (child->loop_type == LoopType::Pseudoknot) {
            // energy += child->number_of_bands * P->MLintern[child_pair_type];
            continue;
        }

        size_t ci = child->begin;
        size_t cj = child->end;
        unsigned int child_pair_type = ViennaUtils::get_pair_type(sequence[ci], sequence[cj], md);

        int child_n5d;
        int child_n3d;
        switch (md.dangles){
            case 0:
                child_n5d = -1;
                child_n3d = -1;
                break;
            case 2:
                child_n5d = vrna_nucleotide_encode(sequence[ci - 1], &md);
                child_n3d = vrna_nucleotide_encode(sequence[cj + 1], &md);
                break;
            default:
                // Dangle 1 is handled separately above
                THROW_ERROR("Dangle model `" + std::to_string(md.dangles) + "` not supported in multi-branch energy calculation. (HOW DID WE GET HERE?)");
        }
        energy += vrna_E_multibranch_stem(child_pair_type, child_n5d, child_n3d, P);
    }

    return energy;
}

int ViennaFunctions::external_energy(const std::vector<std::shared_ptr<LoopNode>>& children,
                    const std::string& sequence) {

    if (md.dangles == 1) {
        return ViennaDangles::get_external_dangle_1(children, sequence, md);
    }

    int energy = 0;
    for (std::shared_ptr<LoopNode> c : children) {
        if (c->loop_type != LoopType::Pseudoknot) {
            unsigned int pair_type = ViennaUtils::get_pair_type(sequence[c->begin], sequence[c->end], md);
            int n5d = c->begin > 0 ? vrna_nucleotide_encode(sequence[c->begin - 1], &md) : -1;
            int n3d = c->end < sequence.size() - 1 ? vrna_nucleotide_encode(sequence[c->end + 1], &md) : -1;
            if (md.dangles == 0) {
                n5d = -1;
                n3d = -1;
            }
            energy += vrna_E_exterior_stem(pair_type, n5d, n3d, P);
        }
    }
    return energy;
}
} // namespace knotergy