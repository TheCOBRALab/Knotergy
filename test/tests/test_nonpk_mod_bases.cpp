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

TEST(mod_nonpk, single_base) {
    std::string sequence = "P";
    std::string structure = ".";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 0.00, 0.000005);  // Turner 2004
    EXPECT_NEAR(d1, 0.00, 0.000005);  // Turner 2004
    EXPECT_NEAR(d2, 0.00, 0.000005);  // Turner 2004
    EXPECT_NEAR(d3, 0.00, 0.000005);  // Turner 2004
}

// echo -e "AAAAAAAAAA6AAAAAAAAAAUUUUUUUUUUUUUUUUUUUUUUUUUUU" | RNAfold -m
TEST(mod_nonpk, small) {
    std::string sequence = "6U";
    std::string structure = "()";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 100000.00, 0.000005);  // Turner 2004
    EXPECT_NEAR(d1, 100000.00, 0.000005);  // Turner 2004
    EXPECT_NEAR(d2, 100000.00, 0.000005);  // Turner 2004
    EXPECT_NEAR(d3, 100000.00, 0.000005);  // Turner 2004
}

// echo -e "AAAUU6\n(...)." | RNAfold -C  --enforceConstraint  -m
TEST(mod_nonpk, dangle_3_external) {
    std::string sequence = "AAAUU6";
    std::string structure = "(...).";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 6.40, 0.000005);  // Turner 2004
    EXPECT_NEAR(d1, 5.97, 0.000005);  // Turner 2004
    EXPECT_NEAR(d2, 5.97, 0.000005);  // Turner 2004
    EXPECT_NEAR(d3, 5.97, 0.000005);  // Turner 2004
}

// echo -e "A6AAAAUUU\n(((...)))" | RNAfold -C  --enforceConstraint  -m
TEST(mod_nonpk, simple_stack) {
    std::string sequence = "A6AAAAUUU";
    std::string structure = "(((...)))";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 4.65, 0.000005);  // Turner 2004
    EXPECT_NEAR(d1, 4.65, 0.000005);  // Turner 2004
    EXPECT_NEAR(d2, 4.65, 0.000005);  // Turner 2004
    EXPECT_NEAR(d3, 4.65, 0.000005);  // Turner 2004
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

    EXPECT_NEAR(d0, -70.11, 0.000005);  // Turner 2004
    EXPECT_NEAR(d1, -70.11, 0.000005);  // Turner 2004
    EXPECT_NEAR(d2, -70.81, 0.000005);  // Turner 2004
    EXPECT_NEAR(d3, -70.11, 0.000005);  // Turner 2004
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

TEST(mod_nonpk, non_au_terminal_penalty) {
    std::string sequence = "GGGIGAAACCCC";
    std::string structure = "((((....))))";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, -6.20, 0.000005);
    EXPECT_NEAR(d1, -6.20, 0.000005);
    EXPECT_NEAR(d2, -6.20, 0.000005);
    EXPECT_NEAR(d3, -6.20, 0.000005);
}

TEST(mod_nonpk, non_au_terminal_penalty2) {
    std::string sequence = "GGGIGAAAUCCC";
    std::string structure = "((((....))))";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, -5.17, 0.000005);
    EXPECT_NEAR(d1, -5.17, 0.000005);
    EXPECT_NEAR(d2, -5.17, 0.000005);
    EXPECT_NEAR(d3, -5.17, 0.000005);
}

TEST(mod_nonpk, pseudoknot_in_external_loop) {
    std::string sequence = "6CGCAAAAGAGAAGCGUAAACUCU";
    std::string structure = "((((...[[[[..))))...]]]]";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, -8.6091, 0.000005);
    EXPECT_NEAR(d1, -8.6091, 0.000005);
    EXPECT_NEAR(d2, -8.6091, 0.000005);
    EXPECT_NEAR(d3, -8.6091, 0.000005);
}

TEST(mod_nonpk, multiloop_with_unpaired_both_sides) {
    std::string sequence = "GACGGAAACCGUGGAAACC6AC";
    std::string structure = "(.(((...)))(((...))).)";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 7.16, 0.000005);
    EXPECT_NEAR(d1, 5.26, 0.000005);
    EXPECT_NEAR(d2, 4.16, 0.000005);
    EXPECT_NEAR(d3, 3.77, 0.000005);
}

TEST(mod_nonpk, multiloop_closing_dangle5) {
    std::string sequence = "UAGGGAAACCCAAGGGAAACCCA6A";
    std::string structure = "(.(((...)))..(((...)))..)";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 5.20, 0.000005);
    EXPECT_NEAR(d1, 0.76, 0.000005);
    EXPECT_NEAR(d2, 1.20, 0.000005);
    EXPECT_NEAR(d3, -0.64, 0.000005);
}

TEST(mod_nonpk, multiloop_closing_dangle3) {
    std::string sequence = "C6GGGAAACCCAAGGGAAACCCAG";
    std::string structure = "(.(((...)))..(((...))).)";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 4.20, 0.000005);
    EXPECT_NEAR(d1, -1.34, 0.000005);
    EXPECT_NEAR(d2, -0.30, 0.000005);
    EXPECT_NEAR(d3, -1.54, 0.000005);
}

TEST(mod_nonpk, multiloop_one_unpaired_between_each_closing_pair) {
    std::string sequence = "GGAIGCAAAAGCCAGCGAAAACGCACC";
    std::string structure = "((.(((....))).(((....))).))";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 0.39, 0.000005);
    EXPECT_NEAR(d1, -4.11, 0.000005);
    EXPECT_NEAR(d2, -3.71, 0.000005);
    EXPECT_NEAR(d3, -4.03, 0.000005);
}

// Terminal mismatches of hairpin and internal loops, and helix ends of internal loops with
// unpaired bases on both sides. Expected values from ViennaRNA's development version (test_dev,
// 2d5a5c3), e.g. echo -e "GGC6AACGCC\n(((....)))" | RNAfold -C --enforceConstraint -m -d2
// Under dangle model 0, ViennaRNA zeroes its multiloop mismatch table, which is the baseline of the
// mismatch corrections, and so adds the whole modified mismatch energy (-4.45 here). Knotergy keeps
// the baseline of the other dangle models.
TEST(mod_nonpk, hairpin_mismatch) {
    std::string sequence = "GGC6AACGCC";
    std::string structure = "(((....)))";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, -2.95, 0.000005);
    EXPECT_NEAR(d1, -2.95, 0.000005);
    EXPECT_NEAR(d2, -2.95, 0.000005);
    EXPECT_NEAR(d3, -2.95, 0.000005);
}

// m6A in the terminal mismatch of the closing pair of a 2x2 internal loop
TEST(mod_nonpk, internal_loop_closing_mismatch) {
    std::string sequence = "GGGC6AGCGCAAAAGCGCACGCCC";
    std::string structure = "((((..((((....))))..))))";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, -13.85, 0.000005);
    EXPECT_NEAR(d1, -13.85, 0.000005);
    EXPECT_NEAR(d2, -13.85, 0.000005);
    EXPECT_NEAR(d3, -13.85, 0.000005);
}

// m6A in the terminal mismatch of the enclosed pair of a 2x2 internal loop, outside both pairs
TEST(mod_nonpk, internal_loop_enclosed_mismatch) {
    std::string sequence = "GGGCACGCGCAAAAGCGC6AGCCC";
    std::string structure = "((((..((((....))))..))))";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, -13.85, 0.000005);
    EXPECT_NEAR(d1, -13.85, 0.000005);
    EXPECT_NEAR(d2, -13.85, 0.000005);
    EXPECT_NEAR(d3, -13.85, 0.000005);
}

// A 1x1 internal loop gets mismatch corrections
TEST(mod_nonpk, internal_loop_1x1_mismatch) {
    std::string sequence = "GGGC6GCGCAAAAGCGCCGCCC";
    std::string structure = "((((.((((....)))).))))";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, -14.95, 0.000005);
    EXPECT_NEAR(d1, -14.95, 0.000005);
    EXPECT_NEAR(d2, -14.95, 0.000005);
    EXPECT_NEAR(d3, -14.95, 0.000005);
}

// Inosine in the closing pair of a 2x2 internal loop: helix-end correction (and stacking)
TEST(mod_nonpk, internal_loop_helix_end) {
    std::string sequence = "GGGIAAGCGCAAAAGCGCAACCCC";
    std::string structure = "((((..((((....))))..))))";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, -13.10, 0.000005);
    EXPECT_NEAR(d1, -13.10, 0.000005);
    EXPECT_NEAR(d2, -13.10, 0.000005);
    EXPECT_NEAR(d3, -13.10, 0.000005);
}

// A 1x1 internal loop gets no helix-end correction (only the stacking correction applies)
TEST(mod_nonpk, internal_loop_1x1_no_helix_end) {
    std::string sequence = "GGGIAGCGCAAAAGCGCACCCC";
    std::string structure = "((((.((((....)))).))))";
    auto [d0, d1, d2, d3] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, -13.42, 0.000005);
    EXPECT_NEAR(d1, -13.42, 0.000005);
    EXPECT_NEAR(d2, -13.42, 0.000005);
    EXPECT_NEAR(d3, -13.42, 0.000005);
}

// echo -e "AAAAAUUUUUUUUPGGGGGGGGGCCCCCCCCC" | RNAfold
// -mod-file=./params/modified_bases/rna_mod_pseudouridine_parameters.json