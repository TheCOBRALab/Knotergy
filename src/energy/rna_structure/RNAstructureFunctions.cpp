#include "energy/rna_structure/RNAstructureFunctions.hpp"

namespace knotergy {

int RNAstructureFunctions::stack_energy(std::size_t i, std::size_t j, std::size_t ci,
                                        std::size_t cj, const std::string& sequence,
                                        const std::vector<std::string_view>& mod_sequence,
                                        const rnastructure_param& rsp) {
    const rnastructure::Datatable& d = rsp.data;
    rnastructure::Code l_code =
        RNAstructureUtils::rnastructure_code(i, sequence, mod_sequence, rsp);
    rnastructure::Code r_code =
        RNAstructureUtils::rnastructure_code(j, sequence, mod_sequence, rsp);
    rnastructure::Code cl_code =
        RNAstructureUtils::rnastructure_code(ci, sequence, mod_sequence, rsp);
    rnastructure::Code cr_code =
        RNAstructureUtils::rnastructure_code(cj, sequence, mod_sequence, rsp);

    return RNAstructureUtils::from_rnastructure(d.stack(l_code, r_code, cl_code, cr_code));
}

}  // namespace knotergy