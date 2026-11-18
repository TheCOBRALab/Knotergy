#pragma once

#include "io/parameters/RNAstructureParams.hpp"
#include "utils/common.hpp"

#include <ViennaRNA/params/constants.hpp>

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace knotergy {

class RNAstructureUtils {
   public:
    [[nodiscard]] static rnastructure::Code rnastructure_code(
        std::size_t k, const std::string& sequence,
        const std::vector<std::string_view>& mod_sequence, const rnastructure_param& rsp) {
        const rnastructure::Datatable& d = rsp.data;
        // mod_sequence is empty when the sequence has no modified bases.
        const std::string_view symbol =
            mod_sequence.empty() ? std::string_view(sequence.data() + k, 1) : mod_sequence[k];
        const rnastructure::Code code = d.code(symbol);
        if (code != rnastructure::kInvalidCode && !d.is_linker(code)) return code;

        // RNAstructure does not have a fallback for modified nucleotides, so check json first.
        THROW_ERROR("Unknown or ambiguous symbol: " + std::string(symbol) + " at position " +
                    std::to_string(k));
    }

    // RNAstructure marks a forbidden table entry with INFINITE_ENERGY (1400 kcal/mol) and adds it
    // up like any other term, so a forbidden loop comes out near a multiple of that (efn2 prints
    // 1399.9, 1400.5, 2800.5, ...). No allowed loop is anywhere near half of it.
    [[nodiscard]] static int from_rnastructure(int energy) {
        return energy >= rnastructure::kInfiniteEnergy / 2 ? INF : 10 * energy;
    }
};
}  // namespace knotergy
