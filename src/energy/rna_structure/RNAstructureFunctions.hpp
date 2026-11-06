#include "energy/rna_structure/RNAstructureUtils.hpp"
#include "io/parameters/RNAstructureParams.hpp"
#include "utils/common.hpp"
namespace rs = thermorna::rnastructure;

namespace knotergy {

class RNAstructureFunctions {
   public:
    [[nodiscard]] static int stack_energy(std::size_t i, std::size_t j, std::size_t ci,
                                          std::size_t cj, const std::string& sequence,
                                          const std::vector<std::string_view>& mod_sequence,
                                          const rnastructure_param& rsp);
};

}  // namespace knotergy