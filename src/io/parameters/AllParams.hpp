
#include "io/parameters/ModParams.hpp"
#include "io/parameters/PseudoknotParams.hpp"
#include "io/parameters/RNAstructureParams.hpp"
#include "io/parameters/ViennaParams.hpp"

namespace knotergy {

class AllParams {
   public:
    AllParams(const all_mod_params& mod, const pk_param& pk, const vrna_md_param& vienna,
              const rnastructure_param& rnastructure)
        : mp_(mod), pkp_(pk), vp_(vienna), rsp_(rnastructure) {}

    const all_mod_params mp_;
    const pk_param pkp_;
    vrna_md_param vp_;
    const rnastructure_param rsp_;

    bool rnastructure_params_loaded = false;
    bool vienna_params_loaded = false;
    bool pseudoknot_params_loaded = false;
    bool modified_params_loaded = false;

    [[nodiscard]] int stack(const std::size_t i, const std::size_t j, const std::size_t ci,
                            const std::size_t cj, const std::string& sequence,
                            const std::vector<std::string_view>& mod_sequence);
    [[nodiscard]] int d5();
    [[nodiscard]] int d3();
    [[nodiscard]] int mismatch();
    [[nodiscard]] int terminal();
    [[nodiscard]] int hairpin();
    [[nodiscard]] int internal();
    [[nodiscard]] int multibranch();
};

}  // namespace knotergy