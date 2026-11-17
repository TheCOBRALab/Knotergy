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

    [[nodiscard]] static int from_rnastructure(int energy) {
        return energy >= rnastructure::kInfiniteEnergy ? INF : 10 * energy;
    }
};
}  // namespace knotergy
