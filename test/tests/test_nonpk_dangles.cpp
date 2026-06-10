#include "test_utils.hpp"

namespace {

const std::string DP_file =
    std::string(KNOTERGY_SOURCE_DIR) + "/params/common/rna_DirksPierce09.par";
const std::string turner_file =
    std::string(KNOTERGY_SOURCE_DIR) + "/params/common/rna_turner2004.par";
const std::string pkp_file =
    std::string(KNOTERGY_SOURCE_DIR) + "/params/pseudo/rna_pk_DirksPierce09_HotKnotsV2.json";
const bool round = false;

std::tuple<double, double, double> pipeline(std::string sequence, std::string structure,
                                            std::string param_file = turner_file) {
    double d0 = get_energy(sequence, structure, 0, round, param_file);
    double d1 = get_energy(sequence, structure, 1, round, param_file);
    double d2 = get_energy(sequence, structure, 2, round, param_file);

    return {d0, d1, d2};
}
}  // namespace

TEST(Dangles, empty) {
    std::string sequence = "U";
    std::string structure = ".";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 0, 0.000005);
    EXPECT_NEAR(d1, 0, 0.000005);
    EXPECT_NEAR(d2, 0, 0.000005);
}

TEST(Dangles, external_simple) {
    std::string sequence = "AAAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUU";
    std::string structure = ".(((((.........................))))).";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 3.90, 0.000005);
    EXPECT_NEAR(d1, 3.2, 0.000005);
    EXPECT_NEAR(d2, 3.2, 0.000005);
}

TEST(Dangles, external_simple_NoEnds) {
    std::string sequence = "AAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUU";
    std::string structure = "(((((.........................)))))";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 3.90, 0.000005);
    EXPECT_NEAR(d1, 3.90, 0.000005);
    EXPECT_NEAR(d2, 3.90, 0.000005);
}

TEST(Dangles, external_simple_adjacent) {
    std::string sequence = "AAAAAAUUUUUGGGGGCCCCCCAAAAAUUUUUUAAAAUUUUGGGCCAAAAAAAUUUUU";
    std::string structure = ".((((((((((((((...)))))))))))))).((((((((........)))))))).";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, -16.40, 0.000005);
    EXPECT_NEAR(d1, -17.20, 0.000005);
    EXPECT_NEAR(d2, -17.60, 0.000005);
}

TEST(Dangles, multiloop_left_dangle) {
    std::string sequence = "AAAAAAAAAUUUUUUUAAAAAUUUUUUU";
    std::string structure = "(((.(((....)))..((...))..)))";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 15.10, 0.000005);
    EXPECT_NEAR(d1, 13.5, 0.000005);
    EXPECT_NEAR(d2, 12.9, 0.000005);
}

TEST(Dangles, multiloop_left_dangle_chained) {
    std::string sequence = "AAAAAAAAAUUUUUUAAAAAUUUUUUU";
    std::string structure = "(((.(((....))).((...))..)))";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 15.10, 0.000005);
    EXPECT_NEAR(d1, 13.6, 0.000005);
    EXPECT_NEAR(d2, 12.9, 0.000005);
}

TEST(Dangles, multiloop_left_touch) {
    std::string sequence = "AAAAAAAAUUUUUUAAAAAUUUUUUU";
    std::string structure = "((((((....))).((...))..)))";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 15.10, 0.000005);
    EXPECT_NEAR(d1, 14.4, 0.000005);
    EXPECT_NEAR(d2, 12.9, 0.000005);
}

TEST(Dangles, multiloop_right_dangle) {
    std::string sequence = "AAAAAAAAAUUUUUUAAAAAUUUUUUU";
    std::string structure = "(((..(((....))).((...)).)))";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 15.10, 0.000005);
    EXPECT_NEAR(d1, 13.1, 0.000005);
    EXPECT_NEAR(d2, 12.4, 0.000005);
}

TEST(Dangles, multiloop_right_touch) {
    std::string sequence = "AAAAAAAAAUUUUUUAAAAAUUUUUU";
    std::string structure = "(((..(((....))).((...)))))";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 15.10, 0.000005);
    EXPECT_NEAR(d1, 13.3, 0.000005);
    EXPECT_NEAR(d2, 12.4, 0.000005);
}

TEST(Dangles, multiloop_both_touch) {
    std::string sequence = "AAAAAAAAAUUUUUUAAAAAUUUUUU";
    std::string structure = "((((((....)))...((...)))))";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 15.10, 0.000005);
    EXPECT_NEAR(d1, 14.70, 0.000005);
    EXPECT_NEAR(d2, 12.70, 0.000005);
}

TEST(Dangles, multiloop_both_RDangle) {
    std::string sequence = "AAAAAAAAAUUUUUUAAAAAUUUUUUU";
    std::string structure = "((((((....)))...((...)).)))";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 15.10, 0.000005);
    EXPECT_NEAR(d1, 14.30, 0.000005);
    EXPECT_NEAR(d2, 12.70, 0.000005);
}

TEST(Dangles, multiloop_both_LDangle) {
    std::string sequence = "AAAAAAAAAUUUUUUAAAAAUUUUUUU";
    std::string structure = "(((.(((....)))...((...)))))";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 15.10, 0.000005);
    EXPECT_NEAR(d1, 13.90, 0.000005);
    EXPECT_NEAR(d2, 12.70, 0.000005);
}

TEST(Dangles, multiloop_both_BDangle) {
    std::string sequence = "AAAAAAAAAUUUUUUAAAAAUUUUUUUU";
    std::string structure = "(((.(((....)))...((...)).)))";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 15.10, 0.000005);
    EXPECT_NEAR(d1, 13.50, 0.000005);
    EXPECT_NEAR(d2, 12.70, 0.000005);
}

TEST(Dangles, multiloop_loop_touch) {
    std::string sequence = "AAAAAAAAAUUUUUAAAAAUUUUU";
    std::string structure = "((((((....))).((...)))))";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 15.10, 0.000005);
    EXPECT_NEAR(d1, 14.90, 0.000005);
    EXPECT_NEAR(d2, 12.90, 0.000005);
}

TEST(Dangles, multiloop_loop_RDangle) {
    std::string sequence = "AAAAAAAAAUUUUAAAAAUUUUUUU";
    std::string structure = "((((((....))).((...)).)))";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 15.10, 0.000005);
    EXPECT_NEAR(d1, 14.20, 0.000005);
    EXPECT_NEAR(d2, 12.40, 0.000005);
}

TEST(Dangles, multiloop_loop_LDangle) {
    std::string sequence = "AAAAAAAAAUUUUUAAAAUUUUUUU";
    std::string structure = "(((.(((....))).((...)))))";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 15.10, 0.000005);
    EXPECT_NEAR(d1, 13.60, 0.000005);
    EXPECT_NEAR(d2, 12.40, 0.000005);
}

TEST(Dangles, multiloop_loop_BDangle) {
    std::string sequence = "AAAAAAAAAUUUUUAAAAUUUUUUUU";
    std::string structure = "(((.(((....))).((...)).)))";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 15.10, 0.000005);
    EXPECT_NEAR(d1, 13.40, 0.000005);
    EXPECT_NEAR(d2, 12.40, 0.000005);
}

// echo -e
// "GGUUUUUUUUAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUGGGGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCC\n.(((((((((((((((..........))))))((((((...)))))))))))))).((((((((........)))))))))."
// | RNAeval -d
TEST(Dangles, external_multiloop) {
    std::string sequence =
        "GGUUUUUUUUAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUGGGGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCC";
    std::string structure =
        ".(((((((((((((((..........))))))((((((...)))))))))))))).((((((((........))))))))).";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, -2.50, 0.000005);
    EXPECT_NEAR(d1, -4.30, 0.000005);
    EXPECT_NEAR(d2, -8.90, 0.000005);
}