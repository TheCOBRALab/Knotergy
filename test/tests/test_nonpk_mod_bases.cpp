#include "utils.hpp"

#include <energy/ComputeEnergy.hpp>
#include <gtest/gtest.h>
#include <io/RNAInputManager.hpp>
#include <loop_tree/LoopFactory.hpp>
#include <preprocessing/ProcessedRNAEntry.hpp>
#include <preprocessing/RNAEntry.hpp>
#include <preprocessing/RNAProcessor.hpp>
#include <utils/common.hpp>

#include <string>
#include <vector>

namespace {

const bool round = false;

std::tuple<double, double, double> dangle_get_energy(std::string sequence, std::string structure,
                                                     std::string param_file = turner_file,
                                                     std::string pseudoknot_param_file = pkp_file,
                                                     std::string mod_param_path = mod_folder) {
    double d0 = get_energy(sequence, structure, 0, round, param_file, pseudoknot_param_file,
                           mod_param_path);
    double d1 = get_energy(sequence, structure, 1, round, param_file, pseudoknot_param_file,
                           mod_param_path);
    double d2 = get_energy(sequence, structure, 2, round, param_file, pseudoknot_param_file,
                           mod_param_path);

    return {d0, d1, d2};
}
}  // namespace

TEST(mod_nonpk, single_base) {
    std::string sequence = "P";
    std::string structure = ".";
    auto [d0, d1, d2] = dangle_get_energy(sequence, structure);

    EXPECT_NEAR(d0, 0.00, 0.000005);  // Turner 2004
    EXPECT_NEAR(d1, 0.00, 0.000005);  // Turner 2004
    EXPECT_NEAR(d2, 0.00, 0.000005);  // Turner 2004
}

// echo -e "AAAAAAAAAA6AAAAAAAAAAUUUUUUUUUUUUUUUUUUUUUUUUUUU" | RNAfold -m
TEST(mod_nonpk, small) {
    std::string sequence = "6U";
    std::string structure = "()";
    auto [d0, d1, d2] = dangle_get_energy(sequence, structure);

    EXPECT_NEAR(d0, 100000.00, 0.000005);  // Turner 2004
    EXPECT_NEAR(d1, 100000.00, 0.000005);  // Turner 2004
    EXPECT_NEAR(d2, 100000.00, 0.000005);  // Turner 2004
}

// echo -e "AAAUU6\n(...)." | RNAfold -C  --enforceConstraint  -m
TEST(mod_nonpk, dangle_3_external) {
    std::string sequence = "AAAUU6";
    std::string structure = "(...).";
    auto [d0, d1, d2] = dangle_get_energy(sequence, structure);

    // RNAfold gives different values due to a likely bug in their code
    // d0 has a dangle modification when there shouldn't be one
    // d1 & d2 use the wrong unmodified energy for the dangle and wrong keys for the modification
    EXPECT_NEAR(d0, 6.40, 0.000005);  // Turner 2004 // RNAfold gives 6.77
    EXPECT_NEAR(d1, 5.70, 0.000005);  // Turner 2004 // RNAfold gives 6.07
    EXPECT_NEAR(d2, 5.70, 0.000005);  // Turner 2004 // RNAfold gives 6.07
}

// echo -e "A6AAAAUUU\n(((...)))" | RNAfold -C  --enforceConstraint  -m
TEST(mod_nonpk, simple_stack) {
    std::string sequence = "A6AAAAUUU";
    std::string structure = "(((...)))";
    auto [d0, d1, d2] = dangle_get_energy(sequence, structure);

    // ViennaRNA's RNAfold gives 4.65 but it's likely a bug in their code
    // (Vienna checks 2 different keys for the stack)
    EXPECT_NEAR(d0, 4.58, 0.000005);  // Turner 2004 // RNAfold gives 4.65
    EXPECT_NEAR(d1, 4.58, 0.000005);  // Turner 2004 // RNAfold gives 4.65
    EXPECT_NEAR(d2, 4.58, 0.000005);  // Turner 2004 // RNAfold gives 4.65
}

// // echo -e "GUUUUUAAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUUUUUUUU66AAAC" | RNAfold -m
// TEST(mod_nonpk, stack2_0_dangle) {
//     std::string sequence = "GUUUUUAAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUUUUUUUU66AAAC";
//     std::string structure= "(((((((((((((((((((((((((......)))))))))))))))))))))))))";
//     auto [d0, d1, d2] = dangle_get_energy(sequence, structure);

//     // EXPECT_NEAR(d0, -18.55, 0.000005); // Turner 2004
//     // EXPECT_NEAR(d1, -18.55, 0.000005); // Turner 2004
//     EXPECT_NEAR(d2, -18.55, 0.000005); // Turner 2004
// }

// echo -e "AAAAAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUUPGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCC"
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
    auto [d0, d1, d2] = dangle_get_energy(sequence, structure);

    // RNAfold gives different values due to a likely bug in their code
    // where they check 2 different keys for the stack
    EXPECT_NEAR(d0, -69.39, 0.000005);  // Turner 2004 // RNAfold gives -70.11
    EXPECT_NEAR(d1, -69.39, 0.000005);  // Turner 2004 // RNAfold gives -70.11
    EXPECT_NEAR(d2, -70.09, 0.000005);  // Turner 2004 // RNAfold gives -70.81
}

// // echo -e
// "AAAAAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUUU6GGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCC\n(((((((((((((((((((....))))))))))))))))))).((((((((((((((((((((....))))))))))))))))))))"
// | RNAfold -m
// // ./build/Knotergy -p ./params/common/rna_turner2004.par -s
// "AAAAAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUUU6GGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCC" -r
// "(((((((((((((((((((....))))))))))))))))))).((((((((((((((((((((....))))))))))))))))))))"
// TEST(mod_nonpk, external2_1_dangle) {
//     std::string sequence =
//     "AAAAAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUUU6GGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCC";
//     std::string structure=
//     "(((((((((((((((((((....))))))))))))))))))).((((((((((((((((((((....))))))))))))))))))))";
//     auto [d0, d1, d2] = dangle_get_energy(sequence, structure);

//     EXPECT_NEAR(d0, -68.83, 0.000005); // Turner 2004
//     EXPECT_NEAR(d1, -69.53, 0.000005); // Turner 2004
//     EXPECT_NEAR(d2, -69.73, 0.000005); // Turner 2004
// }

// echo -e "AAAAAUUUUUUUUPGGGGGGGGGCCCCCCCCC" | RNAfold
// -mod-file=./params/modified_bases/rna_mod_pseudouridine_parameters.json