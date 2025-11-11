#include <gtest/gtest.h>

#include <pipeline/shared.hpp>
#include <preprocessing/RNAEntry.hpp>
#include <preprocessing/ProcessedRNAEntry.hpp>
#include <preprocessing/RNAProcessor.hpp>
#include <string>
#include <vector>

TEST(RNAEntryBasics, populate) {
    std::string name = "Test Input";
    std::string sequence = "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    std::string structure = "(((((.........................)))))................";
    knotergy::RNAEntry rna(name, sequence, structure);

    EXPECT_EQ(rna.name, "Test Input");
    EXPECT_EQ(rna.sequence, "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC");
    EXPECT_EQ(rna.structure, "(((((.........................)))))................");
}

TEST(ProcessedRNAEntryPairings, SimplePairing) {
    std::string sequence = "AAUU";
    std::string structure = "(..)";
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    std::vector<size_t> expected_pairings = {3, SIZE_MAX, SIZE_MAX, 0};

    EXPECT_EQ(processed_rna.get_pairings(), expected_pairings);
}

TEST(RNAEntryClosedRegions, SimpleRegion) {
    std::string sequence = "AAUU";
    std::string structure = "(..)";
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    std::vector<knotergy::ClosedRegion> expected_pairings = {{0, 3}};
    EXPECT_EQ(processed_rna.get_closed_regions(), expected_pairings);
}

TEST(RNAEntryClosedRegions, SimplePseudoKnot) {
    std::string sequence = "AAUU";
    std::string structure = "[(])";
    std::string structure2 = "[(])";
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::RNAEntry rna2(sequence, structure2);
    knotergy::ProcessedRNAEntry processed_rna1(knotergy::RNAProcessor::process_rna(rna));
    knotergy::ProcessedRNAEntry processed_rna2(knotergy::RNAProcessor::process_rna(rna));
    std::vector<knotergy::ClosedRegion> expected_pairings = {{0, 3}};
    EXPECT_EQ(processed_rna1.get_closed_regions(), processed_rna2.get_closed_regions());
    EXPECT_EQ(processed_rna1.get_closed_regions(), expected_pairings);
}

TEST(RNAEntryClosedRegions, PseudoknotWithNestedBP) {
    std::string sequence = "AAGGUAGUU";
    std::string structure = "[(..().])";
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    std::vector<knotergy::ClosedRegion> expected_pairings = {{4, 5}, {0, 8}};
    EXPECT_EQ(processed_rna.get_closed_regions(), expected_pairings);
}

TEST(RNAEntryClosedRegions, PseudoknotWithNestedPseudoknot) {
    std::string sequence = "AAGGTTGGUUU";
    std::string structure = "[(..([)].])";
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    std::vector<knotergy::ClosedRegion> expected_pairings = {{4, 7}, {0, 10}};
    EXPECT_EQ(processed_rna.get_closed_regions(), expected_pairings);
}

TEST(RNAEntryClosedRegions, ComplexPseudoKnot) {
    std::string sequence = "AGUAGUAUTTAAAAGGGGUUUUATTTTUGAUAATU";
    std::string structure = "([)([)()].(((([[[[))))(]]]]).().(])";
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    std::vector<knotergy::ClosedRegion> expected_pairings = {
        {6, 7}, {3, 8}, {10, 27}, {29, 30}, {0, 34}};
    EXPECT_EQ(processed_rna.get_closed_regions(), expected_pairings);
}

TEST(RNAEntryClosedRegions, FauxPseudoKnots) {
    std::string sequence = "AAAUAUUAAUUU";
    std::string structure = "(([][])[()])";
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    std::vector<knotergy::ClosedRegion> expected_pairings = {
        {2, 3}, {4, 5}, {1, 6}, {8, 9}, {7, 10}, {0, 11}};
    EXPECT_EQ(processed_rna.get_closed_regions(), expected_pairings);
}