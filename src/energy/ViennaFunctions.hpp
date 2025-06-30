#pragma once

#include "../rna_regions/RNAEntry.hpp"

extern "C" {
#include <ViennaRNA/eval/hairpin.h>
#include <ViennaRNA/mfe/exterior.h>
#include <ViennaRNA/mfe/internal.h>
#include <ViennaRNA/mfe/multibranch.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/sequences/alphabet.h>
#include <ViennaRNA/utils/basic.h>
}

namespace knotergy {
class ViennaFunctions {
   public:
    ViennaFunctions() { vrna_md_set_default(&md); P = vrna_params(&md);}
    ~ViennaFunctions() { free(P); }

    int stack_energy(size_t i, size_t j, std::string& sequence) {
        // Placeholder for actual energy calculation logic
        return 0.0f;  // Replace with actual energy calculation
    }

    int hairpin_energy(size_t i, size_t j, const std::string& sequence) {
        unsigned int size = static_cast<int>(j - i - 1);
        if (size < 3) {
            std::cerr << "Hairpin loop size is less than 3. Infinite Energy. Sequence: " << sequence
                      << " i: " << i << ", j: " << j << std::endl;
        }

        unsigned int type = get_pair_type(sequence[i], sequence[j]);
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

        int energy = vrna_E_hairpin(size, type, si1, sj1, loop_seq, P);
        return energy;
    }

    int internal_loop_energy(size_t i, size_t j, std::string& sequence) {
        // Placeholder for actual energy calculation logic
        return 0.0f;  // Replace with actual energy calculation
    }

    int multi_energy(size_t i, std::string& sequence) {
        // Placeholder for actual energy calculation logic
        return 0.0f;  // Replace with actual energy calculation
    }

    int pseudoknot_energy(size_t i, size_t j, std::string& sequence) {
        // Placeholder for actual energy calculation logic
        return 0.0f;  // Replace with actual energy calculation
    }

    int external_energy(size_t i, size_t j, std::string& sequence) {
        // Placeholder for actual energy calculation logic
        return 0.0f;  // Replace with actual energy calculation
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