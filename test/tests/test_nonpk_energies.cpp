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

// echo -e "AU\n()" | RNAeval
TEST(NonPseudoKnottedEnergies, nonPK_small_Turner) {
    std::string sequence = "AU";
    std::string structure = "()";
    double turner_energy = get_energy(sequence, structure, turner_file);

    EXPECT_NEAR(turner_energy, 100000.50, 0.000005);  // Turner 2004
}

// echo -e "AU\n()" | RNAeval
TEST(NonPseudoKnottedEnergies, nonPK_small_DP) {
    std::string sequence = "AU";
    std::string structure = "()";
    double dp_energy = get_energy(sequence, structure, DP_file);

    EXPECT_NEAR(dp_energy, 100000.50, 0.000005);  // Dirks & Pierce 2009
}

// echo -e "AGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGU\n(........................................)"
// | RNAeval
TEST(NonPseudoKnottedEnergies, nonPK_hairpin_long_Turner) {
    std::string sequence = "AGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGU";
    std::string structure = "(........................................)";
    double turner_energy = get_energy(sequence, structure, turner_file);

    EXPECT_NEAR(turner_energy, 7.41, 0.000005);  // Turner 2004
}

// echo -e
// "AGGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC\n.(((((.........................)))))................"
// | RNAeval
TEST(NonPseudoKnottedEnergies, BasicStack_Turner) {
    std::string sequence = "AGGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    std::string structure = ".(((((.........................)))))................";

    double turner_energy = get_energy(sequence, structure, turner_file);
    EXPECT_NEAR(turner_energy, -8.30, 0.000005);  // Turner 2004
}

// echo -e
// "AGGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC\n.(((((.........................)))))................"
// | RNAeval
TEST(NonPseudoKnottedEnergies, BasicStack_DP) {
    std::string sequence = "AGGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    std::string structure = ".(((((.........................)))))................";

    double dp_energy = get_energy(sequence, structure, DP_file);
    EXPECT_NEAR(dp_energy, -5.91, 0.000005);  // Dirks & Pierce 2009
}

// echo -e "GGGGAGAAAAAAAAAUUUUUU\n((((((.........))))))" | RNAeval
TEST(NonPseudoKnottedEnergies, StacksFullStructure_Turner) {
    std::string sequence = "GGGGAGAAAAAAAAAUUUUUU";
    std::string structure = "((((((.........))))))";

    double turner_energy = get_energy(sequence, structure, turner_file);
    EXPECT_NEAR(turner_energy, 3.70, 0.000005);
}

// echo -e "GGGGAGAAAAAAAAAUUUUUU\n((((((.........))))))" | RNAeval
TEST(NonPseudoKnottedEnergies, StacksFullStructure_DP) {
    std::string sequence = "GGGGAGAAAAAAAAAUUUUUU";
    std::string structure = "((((((.........))))))";

    double dp_energy = get_energy(sequence, structure, DP_file);
    EXPECT_NEAR(dp_energy, 0.87, 0.000005);  // Dirks & Pierce 2009
}

// echo -e
// "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCC\n((((((.....................))))))................((((((......))))))...(((...)))......"
// | RNAeval
TEST(NonPseudoKnottedEnergies, MultipleStacks_Turner) {
    std::string sequence =
        "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCC";
    std::string structure =
        "((((((.....................))))))................((((((......))))))...(((...)))......";

    double turner_energy = get_energy(sequence, structure, turner_file);
    EXPECT_NEAR(turner_energy, -28.10, 0.000005);
}

// echo -e
// "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCC\n((((((.....................))))))................((((((......))))))...(((...)))......"
// | RNAeval
TEST(NonPseudoKnottedEnergies, MultipleStacks_DP) {
    std::string sequence =
        "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCC";
    std::string structure =
        "((((((.....................))))))................((((((......))))))...(((...)))......";

    double dp_energy = get_energy(sequence, structure, DP_file);
    EXPECT_NEAR(dp_energy, -19.21, 0.000005);  // Dirks & Pierce 2009
}

// echo -e
// "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCCA\n...........((((((....((((((............))))))....((((((......))))))...((....)).))))))."
// | RNAeval
TEST(NonPseudoKnottedEnergies, MultiLoopWithStacks_Turner) {
    std::string sequence =
        "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCCA";
    std::string structure =
        "...........((((((....((((((............))))))....((((((......))))))...((....)).)))))).";

    double turner_energy = get_energy(sequence, structure, turner_file);
    EXPECT_NEAR(turner_energy, -41.10, 0.000005);
}

// echo -e
// "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCCA\n...........((((((....((((((............))))))....((((((......))))))...((....)).))))))."
// | RNAeval
TEST(NonPseudoKnottedEnergies, MultiLoopWithStacks_DP) {
    std::string sequence =
        "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCCA";
    std::string structure =
        "...........((((((....((((((............))))))....((((((......))))))...((....)).)))))).";

    double dp_energy = get_energy(sequence, structure, DP_file);
    EXPECT_NEAR(dp_energy, -29.38, 0.000005);  // Dirks & Pierce 2009
}

// echo -e
// "AAAAAAGGGGCCCCCCAAAAAAAAAAAAGGGGGGAAAAGGGGGGUUUUUUCCCCCCAAAUUUUUUAAGUUUUUU\n((((((....((((((............))))))....((((((......))))))...((....)).))))))"
// | RNAeval
TEST(NonPseudoKnottedEnergies, MultiLoopWithStacksFull_Turner) {
    std::string sequence =
        "AAAAAAGGGGCCCCCCAAAAAAAAAAAAGGGGGGAAAAGGGGGGUUUUUUCCCCCCAAAUUUUUUAAGUUUUUU";
    std::string structure =
        "((((((....((((((............))))))....((((((......))))))...((....)).))))))";

    double turner_energy = get_energy(sequence, structure, turner_file);
    EXPECT_NEAR(turner_energy, -21.9, 0.000005);
}

// echo -e
// "AAAAAAGGGGCCCCCCAAAAAAAAAAAAGGGGGGAAAAGGGGGGUUUUUUCCCCCCAAAUUUUUUAAGUUUUUU\n((((((....((((((............))))))....((((((......))))))...((....)).))))))"
// | RNAeval
TEST(NonPseudoKnottedEnergies, MultiLoopWithStacksFull_DP) {
    std::string sequence =
        "AAAAAAGGGGCCCCCCAAAAAAAAAAAAGGGGGGAAAAGGGGGGUUUUUUCCCCCCAAAUUUUUUAAGUUUUUU";
    std::string structure =
        "((((((....((((((............))))))....((((((......))))))...((....)).))))))";

    double dp_energy = get_energy(sequence, structure, DP_file);
    EXPECT_NEAR(dp_energy, -16.47, 0.000005);  // Dirks & Pierce 2009
}

// echo -e
// "GGGGUUAUUUUAUUAAAAAUAACCCUGGUUUUUAAGGCGGGGUCGUGCGGUAAGGGAACCC\n((((..(...).((...)))..(((.(.((...))..(...).)..(...)..)))..)))"
// | RNAeval
TEST(NonPseudoKnottedEnergies, MultiWithMultiLoop_Turner) {
    std::string sequence = "GGGGUUAUUUUAUUAAAAAUAACCCUGGUUUUUAAGGCGGGGUCGUGCGGUAAGGGAACCC";
    std::string structure = "((((..(...).((...)))..(((.(.((...))..(...).)..(...)..)))..)))";

    double turner_energy = get_energy(sequence, structure, turner_file);
    EXPECT_NEAR(turner_energy, 30.90, 0.000005);
}

// echo -e
// "GGGGUUAUUUUAUUAAAAAUAACCCUGGUUUUUAAGGCGGGGUCGUGCGGUAAGGGAACCC\n((((..(...).((...)))..(((.(.((...))..(...).)..(...)..)))..)))"
// | RNAeval
TEST(NonPseudoKnottedEnergies, MultiWithMultiLoop_DP) {
    std::string sequence = "GGGGUUAUUUUAUUAAAAAUAACCCUGGUUUUUAAGGCGGGGUCGUGCGGUAAGGGAACCC";
    std::string structure = "((((..(...).((...)))..(((.(.((...))..(...).)..(...)..)))..)))";

    double dp_energy = get_energy(sequence, structure, DP_file);
    EXPECT_NEAR(dp_energy, 19.81, 0.000005);  // Dirks & Pierce 2009
}

TEST(NonPseudoKnottedEnergies, PositiveInternalLoops_Turner) {
    std::string sequence =
        "GGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAA"
        "CCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUU"
        "GGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAAUUUCGAGGAUUC"
        "GAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGA"
        "GGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGG"
        "AUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGCGUACGUACGUACGUACGU"
        "ACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUAC"
        "GUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGU"
        "ACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACG";
    std::string structure =
        "(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((..."
        "(((...(((...(((...(((...(((...(((...)))...)))...)))...)))...)))...)))...)))...)))...)))..."
        ")))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))..........((((((..((("
        "(((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((("
        "(..((((((..((((((..((((((..........))))))..))))))..))))))..))))))..))))))..))))))..))))))."
        ".))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))((((((((((((((((((("
        "(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((("
        "((((((((((((((((((((((((....))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))"
        ")))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))";

    double turner_energy = get_energy(sequence, structure, turner_file);
    EXPECT_NEAR(turner_energy, -471.90, 0.000005);
}

TEST(NonPseudoKnottedEnergies, PositiveInternalLoops_DP) {
    std::string sequence =
        "GGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAA"
        "CCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUU"
        "GGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAAUUUCGAGGAUUC"
        "GAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGA"
        "GGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGG"
        "AUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGCGUACGUACGUACGUACGU"
        "ACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUAC"
        "GUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGU"
        "ACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACG";
    std::string structure =
        "(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((..."
        "(((...(((...(((...(((...(((...(((...)))...)))...)))...)))...)))...)))...)))...)))...)))..."
        ")))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))..........((((((..((("
        "(((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((("
        "(..((((((..((((((..((((((..........))))))..))))))..))))))..))))))..))))))..))))))..))))))."
        ".))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))((((((((((((((((((("
        "(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((("
        "((((((((((((((((((((((((....))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))"
        ")))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))";

    double dp_energy = get_energy(sequence, structure, DP_file);
    EXPECT_NEAR(dp_energy, -317.37, 0.000005);  // Dirks & Pierce 2009
}

TEST(NonPseudoKnottedEnergies, MoreLongSequences_Turner) {
    std::string sequence =
        "UUAAAAGGGAUGCCUCUCCUGUUCAUCUUGUGGAGAAGCAUUCGAUAAGGUCAUCAUAAUGGGUCCAGCUUUGCGACCUGGCGAGAUUAG"
        "UCAGGAAAAUGUGAAGUGGGUCUUCGCUUUCCAGGUACAGGAGGCUCGCCCCGCUCAUCCAGUUCGUCCCCUAACCACUUGUUUUCUCAG"
        "GAUAGUUUGUUUUGUACACCCGUGUACAUACACAUGUAUCACACCCCAGAUUGCCGAAUGUUUCGUUCGGUCGAGCCUGACUAUGCAUAA"
        "ACCCUACCUCUGAAACCUUGGGCAACUCACUACUUCCGAGCUAAAUCCCUCUGUUUGAGCUAGCCUGAGAUUUCAACUGGCUUCGGCCUU"
        "GUUUAUACCAUCGUUUGCUGAUCCAUUGAAGAAAUAAGUUACCGAUGGCCCCAAACUGACGAUCACUAUUCUUUCCAUAGGAGUUAUGGG"
        "UAUACUGCCCGUAGACGGAAAGAUGAAUGCCUGUAUCCGGGAGUCAGAUG";
    std::string structure =
        ".......((((((.(((((............))))).)))))).....(((((((..((((((((.(((...((((((((((.......)"
        ")))))...(((.(..((((((((((..............))))))))))....).)))....................(((...(((((("
        "(.((((((.............((((((((....))))).)))....((((..((((((((...))))))))...((((((.....((..."
        "..........)).....))))))..(((.........))).........))))...)))))).)))))))...)))..(((....))).."
        "..........))))..)))))))))))..(....)........))))))).....(((((...((....(((((((...(...((((((("
        ".......))))))).))))))))....))((((....)))).)))))...";

    double turner_energy = get_energy(sequence, structure, turner_file);
    EXPECT_NEAR(turner_energy, -121.2, 0.000005);
}

TEST(NonPseudoKnottedEnergies, MoreLongSequences_DP) {
    std::string sequence =
        "UUAAAAGGGAUGCCUCUCCUGUUCAUCUUGUGGAGAAGCAUUCGAUAAGGUCAUCAUAAUGGGUCCAGCUUUGCGACCUGGCGAGAUUAG"
        "UCAGGAAAAUGUGAAGUGGGUCUUCGCUUUCCAGGUACAGGAGGCUCGCCCCGCUCAUCCAGUUCGUCCCCUAACCACUUGUUUUCUCAG"
        "GAUAGUUUGUUUUGUACACCCGUGUACAUACACAUGUAUCACACCCCAGAUUGCCGAAUGUUUCGUUCGGUCGAGCCUGACUAUGCAUAA"
        "ACCCUACCUCUGAAACCUUGGGCAACUCACUACUUCCGAGCUAAAUCCCUCUGUUUGAGCUAGCCUGAGAUUUCAACUGGCUUCGGCCUU"
        "GUUUAUACCAUCGUUUGCUGAUCCAUUGAAGAAAUAAGUUACCGAUGGCCCCAAACUGACGAUCACUAUUCUUUCCAUAGGAGUUAUGGG"
        "UAUACUGCCCGUAGACGGAAAGAUGAAUGCCUGUAUCCGGGAGUCAGAUG";
    std::string structure =
        ".......((((((.(((((............))))).)))))).....(((((((..((((((((.(((...((((((((((.......)"
        ")))))...(((.(..((((((((((..............))))))))))....).)))....................(((...(((((("
        "(.((((((.............((((((((....))))).)))....((((..((((((((...))))))))...((((((.....((..."
        "..........)).....))))))..(((.........))).........))))...)))))).)))))))...)))..(((....))).."
        "..........))))..)))))))))))..(....)........))))))).....(((((...((....(((((((...(...((((((("
        ".......))))))).))))))))....))((((....)))).)))))...";

    double dp_energy = get_energy(sequence, structure, DP_file);
    EXPECT_NEAR(dp_energy, -90.42, 0.000005);  // Dirks & Pierce 2009
}