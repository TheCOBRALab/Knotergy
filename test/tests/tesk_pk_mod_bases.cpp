#include "test_utils.hpp"

namespace {

const int round = 0;

std::tuple<double, double, double, double> pipeline(std::string sequence, std::string structure,
                                                    std::string param_file = turner_file,
                                                    std::string pseudoknot_param_file = pkp_file,
                                                    std::string mod_param_path = mod_folder) {
    double d0 = get_energy(sequence, structure, 0, round, param_file, pseudoknot_param_file,
                           mod_param_path);
    double d1 = get_energy(sequence, structure, 1, round, param_file, pseudoknot_param_file,
                           mod_param_path);
    double d2 = get_energy(sequence, structure, 2, round, param_file, pseudoknot_param_file,
                           mod_param_path);
    double d3 = get_energy(sequence, structure, 3, round, param_file, pseudoknot_param_file,
                           mod_param_path);

    return {d0, d1, d2, d3};
}
}  // namespace

TEST(mod_pk, pk_in_multiloop) {
    std::string sequence = "GAGGGA6AACCCAACGCAAAAGAGAAGCGUAAACUCUAC";
    std::string structure = "(.(((....))).((((...[[[[..))))...]]]].)";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 6.976, 0.000005);  // Turner 2004
    EXPECT_NEAR(d1, 4.176, 0.000005);  // Turner 2004
    EXPECT_NEAR(d2, 4.376, 0.000005);  // Turner 2004
    EXPECT_NEAR(d3, 3.676, 0.000005);  // Turner 2004
}