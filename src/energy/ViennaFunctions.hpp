#pragma once

#include "../loops/LoopNode.hpp"
#include "../rna_regions/RNAEntry.hpp"

extern "C" {
#include <ViennaRNA/eval/exterior.h>
#include <ViennaRNA/eval/hairpin.h>
#include <ViennaRNA/eval/internal.h>
#include <ViennaRNA/eval/multibranch.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/sequences/alphabet.h>
#include <ViennaRNA/utils/basic.h>
}

namespace knotergy {
class ViennaFunctions {
   public:
    ViennaFunctions() {
        vrna_md_set_default(&md);
        P = vrna_params(&md);
    }
    ~ViennaFunctions() { free(P); }

    int stack_energy(size_t i, size_t j, size_t ci, size_t cj, const std::string& sequence) {
        if (j <= i || cj <= ci || ci <= i || j <= cj || j >= sequence.size()) {
            std::cerr << "Invalid indices for stack energy calculation." << std::endl;
            return 0;
        }
        // c = child or nested base pair
        unsigned int type1 = get_pair_type(sequence[i], sequence[j]);
        unsigned int type2 = get_pair_type(sequence[ci], sequence[cj]);
        return P->stack[type1][md.rtype[type2]];
    }

    int hairpin_energy(size_t i, size_t j, const std::string& sequence) {
        if (j <= i || j >= sequence.size()) {
            std::cerr << "Invalid indices for hairpin energy calculation." << std::endl;
            return 0;
        }

        unsigned int size = static_cast<unsigned int>(j - i - 1);
        if (size < 3) {
            std::cerr << "Hairpin loop size is less than 3. Infinite Energy. Sequence: " << sequence
                      << " i: " << i << ", j: " << j << std::endl;
        }

        unsigned int pair_type = get_pair_type(sequence[i], sequence[j]);
        int si1 = vrna_nucleotide_encode(sequence[i + 1], &md);
        int sj1 = vrna_nucleotide_encode(sequence[j - 1], &md);

        // For loop size < 7, you MUST pass the loop sequence substring
        // https://github.com/ViennaRNA/ViennaRNA/blob/219394580aec203a9d6f0d5450021e22642d5a83/src/ViennaRNA/eval/hairpin.h#L78C1-L81C94
        const char* loop_seq = nullptr;
        std::string loop_subseq;
        if (size < 7) {
            loop_subseq = sequence.substr(i, size + 2);  // length == size
            loop_seq = loop_subseq.c_str();              // temporary c-string
        }

        return vrna_E_hairpin(size, pair_type, si1, sj1, loop_seq, P);
    }

    int internal_loop_energy(size_t i, size_t j, size_t ci, size_t cj,
                             const std::string& sequence) {
        // c = child or nested bp

        if (j <= i || cj <= ci || ci <= i || j <= cj || j >= sequence.size()) {
            std::cerr << "Invalid indices for internal loop energy calculation." << std::endl;
            return 0;
        }

        unsigned int n1 = static_cast<unsigned int>(ci - i - 1);
        unsigned int n2 = static_cast<unsigned int>(j - cj - 1);
        unsigned int pair_type = get_pair_type(sequence[i], sequence[j]);
        unsigned int pair_type2 = get_pair_type(sequence[ci], sequence[cj]);
        int si1 = vrna_nucleotide_encode(sequence[i + 1], &md);   // 5' mismatch of closing pair
        int sj1 = vrna_nucleotide_encode(sequence[j - 1], &md);   // 3' mismatch of closing pair
        int sp1 = vrna_nucleotide_encode(sequence[ci - 1], &md);  // 5' mismatch of enclosed pair
        int sq1 = vrna_nucleotide_encode(sequence[cj + 1], &md);  // 3' mismatch of enclosed pair

        return vrna_E_internal(n1, n2, pair_type, pair_type2, si1, sj1, sp1, sq1, P);
    }

    int multibranch_energy(const LoopNode& node, const std::string& sequence) {
        int energy = P->MLclosing;  // closing penalty
        size_t i = node.begin;
        size_t j = node.end;

        unsigned int pair_type = get_pair_type(sequence[i], sequence[j]);
        int n5d = vrna_nucleotide_encode(sequence[i + 1], &md);
        int n3d = vrna_nucleotide_encode(sequence[j - 1], &md);

        energy += vrna_E_multibranch_stem(md.rtype[pair_type], n3d, n5d, P);

        int unpaired = node.number_of_unpaired_bases;
        const std::vector<std::shared_ptr<LoopNode>>& children = node.children;
        for (const std::shared_ptr<LoopNode>& child : children) {
            size_t ci = child->begin;
            size_t cj = child->end;

            unsigned int pair_type = get_pair_type(sequence[ci], sequence[cj]);
            int n5d = ci > 0 ? vrna_nucleotide_encode(sequence[ci - 1], &md) : -1;
            int n3d = cj < sequence.size() - 1 ? vrna_nucleotide_encode(sequence[cj + 1], &md) : -1;

            energy += vrna_E_multibranch_stem(pair_type, n5d, n3d, P);
            unpaired -= child->number_of_unpaired_bases;
        }

        energy += unpaired * P->MLbase;

        return energy;
    }

    int pseudoknot_energy(size_t i, size_t j, const std::string& sequence) {
        // Placeholder for actual energy calculation logic
        return 0;  // Replace with actual energy calculation
    }

    int external_energy(const std::vector<std::shared_ptr<LoopNode>>& children,
                        const std::string& sequence) {
        int energy = 0;

        for (std::shared_ptr<LoopNode> c : children) {
            if (c->loop_type != LoopType::Pseudoknot) {
                unsigned int pair_type = get_pair_type(sequence[c->begin], sequence[c->end]);
                int n5d = c->begin > 0 ? vrna_nucleotide_encode(sequence[c->begin - 1], &md) : -1;
                int n3d = c->end < sequence.size() - 1
                              ? vrna_nucleotide_encode(sequence[c->end + 1], &md)
                              : -1;
                energy += vrna_E_exterior_stem(pair_type, n5d, n3d, P);
            }
        }
        return energy;
    }

   private:
    vrna_md_t md;
    vrna_param_t* P;

    unsigned int get_pair_type(const char& i, const char& j) {
        int encoded_i = vrna_nucleotide_encode(i, &md);
        int encoded_j = vrna_nucleotide_encode(j, &md);
        return vrna_get_ptype_md(encoded_i, encoded_j, &md);
    }
};
}  // namespace knotergy