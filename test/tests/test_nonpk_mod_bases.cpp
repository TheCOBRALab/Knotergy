#include "test_utils.hpp"

namespace {

const bool round = false;

std::tuple<double, double, double, double> pipeline(std::string sequence, std::string structure,
                                                    std::string param_file            = turner_file,
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

TEST(mod_nonpk, single_base) {
    std::string sequence  = "P";
    std::string structure = ".";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 0.00, 0.000005);  // Turner 2004
    EXPECT_NEAR(d1, 0.00, 0.000005);  // Turner 2004
    EXPECT_NEAR(d2, 0.00, 0.000005);  // Turner 2004
    EXPECT_NEAR(d3, 0.00, 0.000005);  // Turner 2004
}

// echo -e "AAAAAAAAAA6AAAAAAAAAAUUUUUUUUUUUUUUUUUUUUUUUUUUU" | RNAfold -m
TEST(mod_nonpk, small) {
    std::string sequence  = "6U";
    std::string structure = "()";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 100000.00, 0.000005);  // Turner 2004
    EXPECT_NEAR(d1, 100000.00, 0.000005);  // Turner 2004
    EXPECT_NEAR(d2, 100000.00, 0.000005);  // Turner 2004
    EXPECT_NEAR(d3, 100000.00, 0.000005);  // Turner 2004
}

// echo -e "AAAUU6\n(...)." | RNAfold -C  --enforceConstraint  -m
TEST(mod_nonpk, dangle_3_external) {
    std::string sequence  = "AAAUU6";
    std::string structure = "(...).";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    // RNAfold gives different values due to a likely bug in their code
    // d0 has a dangle modification when there shouldn't be one
    // d1 & d2 use the wrong unmodified energy for the dangle and wrong keys for the modification
    EXPECT_NEAR(d0, 6.40, 0.000005);  // Turner 2004 // RNAfold gives 6.77
    EXPECT_NEAR(d1, 5.70, 0.000005);  // Turner 2004 // RNAfold gives 6.07
    EXPECT_NEAR(d2, 5.70, 0.000005);  // Turner 2004 // RNAfold gives 6.07
    EXPECT_NEAR(d3, 5.70, 0.000005);  // Turner 2004 // RNAfold gives 6.07
}

// echo -e "A6AAAAUUU\n(((...)))" | RNAfold -C  --enforceConstraint  -m
TEST(mod_nonpk, simple_stack) {
    std::string sequence  = "A6AAAAUUU";
    std::string structure = "(((...)))";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 4.65, 0.000005);  // Turner 2004 // RNAfold gives 4.65
    EXPECT_NEAR(d1, 4.65, 0.000005);  // Turner 2004 // RNAfold gives 4.65
    EXPECT_NEAR(d2, 4.65, 0.000005);  // Turner 2004 // RNAfold gives 4.65
    EXPECT_NEAR(d3, 4.65, 0.000005);  // Turner 2004 // RNAfold gives 4.65
}

// echo -e
// "AAAAAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUUPGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCC"
// | RNAfold -mod-file=./params/modified_bases/rna_mod_pseudouridine_parameters.json
// ./build/Knotergy -s
// "AAAAAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUUPGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCC" -r
// "(((((((((((((((((((....)))))))))))))))))))((((((((((((((((((((....))))))))))))))))))))" -p
// ./params/common/rna_turner2004.par
TEST(mod_nonpk, external1_0_dangle) {
    std::string sequence =
        "AAAAAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUUPGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCC";
    std::string structure =
        "(((((((((((((((((((....)))))))))))))))))))((((((((((((((((((((....))))))))))))))))))))";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, -70.11, 0.000005);  // Turner 2004 // RNAfold gives -70.11
    EXPECT_NEAR(d1, -70.11, 0.000005);  // Turner 2004 // RNAfold gives -70.11
    EXPECT_NEAR(d2, -70.81, 0.000005);  // Turner 2004 // RNAfold gives -70.81
    EXPECT_NEAR(d3, -70.11, 0.000005);  // Turner 2004 // RNAfold gives -70.81
}

TEST(mod_nonpk, multiloop_no_lookup) {
    std::string sequence =
        "GGUUUUUUUUAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUGGGGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCA6";
    std::string structure =
        ".(((((((((((((((..........))))))((((((...)))))))))))))).((((((((........)))))))))...";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, -2.50, 0.000005);
    EXPECT_NEAR(d1, -4.30, 0.000005);
    EXPECT_NEAR(d2, -8.90, 0.000005);
    EXPECT_NEAR(d3, -10.90, 0.000005);
}

// echo -e "AAAAAUUUUUUUUPGGGGGGGGGCCCCCCCCC" | RNAfold
// -mod-file=./params/modified_bases/rna_mod_pseudouridine_parameters.json