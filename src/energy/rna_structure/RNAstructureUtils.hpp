#include "io/parameters/RNAstructureParams.hpp"
#include "utils/common.hpp"

#include <ViennaRNA/params/constants.hpp>
namespace rs = thermorna::rnastructure;

namespace knotergy {

class RNAstructureUtils {
   public:
    [[nodiscard]] static rnastructure::Code rnastructure_code(
        std::size_t k, const std::string& sequence,
        const std::vector<std::string_view>& mod_sequence, const rnastructure_param& rsp) {
        const rnastructure::Datatable& d = rsp.data;
        if (!mod_sequence.empty()) {
            const rnastructure::Code code = d.code(mod_sequence[k]);
            if (code != rnastructure::kInvalidCode && !d.is_linker(code)) return code;
        } else {
            const rnastructure::Code code = d.code(std::string_view(sequence.data() + k, 1));
            if (code != rnastructure::kInvalidCode && !d.is_linker(code)) return code;
        }

        // RNAstructure does not have a fallback for modified nucleotides, so check json first.
        THROW_ERROR("Unknown or ambiguous symbol: " + std::string(mod_sequence[k]) +
                    " at position " + std::to_string(k));
    }

    [[nodiscard]] static int from_rnastructure(int energy) {
        return energy >= rnastructure::kInfiniteEnergy ? INF : 10 * energy;
    }
};
}  // namespace knotergy