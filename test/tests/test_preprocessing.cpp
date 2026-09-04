#include "test_utils.hpp"

#include <gtest/gtest.h>

TEST(RNAEntry, ConstructsWithNameSequenceStructure) {
    std::string name = "Test Input";
    std::string sequence = "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    std::string structure = "(((((.........................)))))................";
    knotergy::RNAEntry rna(name, sequence, structure);

    EXPECT_EQ(rna.name, "Test Input");
    EXPECT_EQ(rna.sequence, "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC");
    EXPECT_EQ(rna.structure, "(((((.........................)))))................");
}

// -------------- Getpair_table tests --------------

TEST(ProcessedRNAEntry, Getpair_table_ReturnsEmptyForUnpaired) {
    std::string sequence = "AAAAA";
    std::string structure = ".....";
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    std::vector<std::size_t> expected_pair_table = {SIZE_MAX, SIZE_MAX, SIZE_MAX, SIZE_MAX,
                                                    SIZE_MAX};
    EXPECT_EQ(processed_rna.get_pair_table(), expected_pair_table);
}

TEST(ProcessedRNAEntry, Getpair_table_ReturnsExpectedForSimpleHairpin) {
    std::string sequence = "AAUU";
    std::string structure = "(..)";
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    std::vector<std::size_t> expected_pair_table = {3, SIZE_MAX, SIZE_MAX, 0};

    EXPECT_EQ(processed_rna.get_pair_table(), expected_pair_table);
}

TEST(ProcessedRNAEntry, Getpair_table_HandlesMultiplePairingNotations) {
    std::string sequence = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU";
    std::string structure = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ)]}>abcdefghijklmnopqrstuvwxyz";
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    std::vector<std::size_t> expected_pair_table = {
        30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
        50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
        10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
    EXPECT_EQ(processed_rna.get_pair_table(), expected_pair_table);
}

// ------------ ClosedRegion tests ------------
TEST(ProcessedRNAEntry, GetClosedRegions_ReturnsEmptyForUnpaired) {
    std::string sequence = "AAAAA";
    std::string structure = ".....";
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    std::vector<knotergy::ClosedRegion> expected_closed_regions = {};
    EXPECT_EQ(processed_rna.get_closed_regions(), expected_closed_regions);
}

TEST(ProcessedRNAEntry, GetClosedRegions_ReturnsExpectedForSimpleHairpin) {
    std::string sequence = "AAUU";
    std::string structure = "(..)";
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    std::vector<knotergy::ClosedRegion> expected_closed_regions = {{0, 3}};
    EXPECT_EQ(processed_rna.get_closed_regions(), expected_closed_regions);
}

TEST(ProcessedRNAEntry, GetClosedRegions_IsDeterministicForPseudoknotNotation) {
    std::string sequence = "AAUU";
    std::string structure = "[(])";
    std::string structure2 = "([)]";
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::RNAEntry rna2(sequence, structure2);
    knotergy::ProcessedRNAEntry processed_rna1(knotergy::RNAProcessor::process_rna(rna));
    knotergy::ProcessedRNAEntry processed_rna2(knotergy::RNAProcessor::process_rna(rna2));
    std::vector<knotergy::ClosedRegion> expected_closed_regions = {{0, 3}};
    EXPECT_EQ(processed_rna1.get_closed_regions(), processed_rna2.get_closed_regions());
    EXPECT_EQ(processed_rna1.get_closed_regions(), expected_closed_regions);
}

TEST(ProcessedRNAEntry, GetClosedRegions_HandlesPseudoknotWithNestedPKBasePair) {
    std::string sequence = "AAGGUAGUU";
    std::string structure = "[(..().])";
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    std::vector<knotergy::ClosedRegion> expected_pair_table = {{0, 8}, {4, 5}};
    EXPECT_EQ(processed_rna.get_closed_regions(), expected_pair_table);
}

TEST(ProcessedRNAEntry, GetClosedRegions_HandlesPseudoknotWithNestedPseudoknot) {
    std::string sequence = "AAGGTTGGUUU";
    std::string structure = "[(..([)].])";
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    std::vector<knotergy::ClosedRegion> expected_closed_regions = {{0, 10}, {4, 7}};
    EXPECT_EQ(processed_rna.get_closed_regions(), expected_closed_regions);
}

TEST(ProcessedRNAEntry, GetClosedRegions_HandlesFauxPseudoknots) {
    std::string sequence = "AAAUAUUAAUUU";
    std::string structure = "(([][])[()])";
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    std::vector<knotergy::ClosedRegion> expected_closed_regions = {{0, 11}, {1, 6},  {2, 3},
                                                                   {4, 5},  {7, 10}, {8, 9}};
    EXPECT_EQ(processed_rna.get_closed_regions(), expected_closed_regions);
}

TEST(ProcessedRNAEntry, GetClosedRegions_HandlesMultiplePairingNotations) {
    std::string sequence = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU";
    std::string structure = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ)]}>abcdefghijklmnopqrstuvwxyz";
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    std::vector<knotergy::ClosedRegion> expected_closed_regions = {{0, 59}};
    EXPECT_EQ(processed_rna.get_closed_regions(), expected_closed_regions);
}

TEST(ProcessedRNAEntry, GetClosedRegions_HandlesComplexPseudoknot) {
    std::string sequence = "AGUAGUAUTTAAAAGGGGUUUUATTTTUGAUAATU";
    std::string structure = "([)([)()].(((([[[[))))(]]]]).().(])";
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    std::vector<knotergy::ClosedRegion> expected_closed_regions = {
        {0, 34}, {3, 8}, {6, 7}, {10, 27}, {29, 30}};
    EXPECT_EQ(processed_rna.get_closed_regions(), expected_closed_regions);
}
