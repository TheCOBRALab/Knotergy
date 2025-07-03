#include <gtest/gtest.h>

#include <helpers/common.hpp>
#include <preprocessing/RNAEntry.hpp>
#include <string>
#include <vector>

std::string name = "Test Input";

TEST(RNAEntryBasics, populate) {
    std::string sequence = "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    std::string structure = "(((((.........................)))))................";
    knotergy::RNAEntry rna(name, sequence, structure);

    EXPECT_EQ(rna.get_name(), "Test Input");
    EXPECT_EQ(rna.get_sequence(), "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC");
    EXPECT_EQ(rna.get_structure(), "(((((.........................)))))................");
}

TEST(RNAEntryBasics, GettersAndSetters) {
    std::string sequence = "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    std::string structure = "(((((.........................)))))................";
    knotergy::RNAEntry rna;

    rna.set_name(name);
    rna.set_sequence(sequence);
    rna.set_structure(structure);

    EXPECT_EQ(rna.get_name(), "Test Input");
    EXPECT_EQ(rna.get_sequence(), "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC");
    EXPECT_EQ(rna.get_structure(), "(((((.........................)))))................");
}

TEST(RNAEntryPairings, SimplePairing) {
    std::string sequence = "AAUU";
    std::string structure = "(..)";
    knotergy::RNAEntry rna(name, sequence, structure);
    std::vector<size_t> expected_pairings = {3, knotergy::NULL_INDEX,
                                             knotergy::NULL_INDEX, 0};

    EXPECT_EQ(rna.get_pairings(), expected_pairings);
}

TEST(RNAEntryClosedRegions, SimpleRegion) {
    std::string sequence = "AAUU";
    std::string structure = "(..)";
    knotergy::RNAEntry rna(name, sequence, structure);
    std::vector<knotergy::ClosedRegion> expected_pairings = {{0, 3}};
    EXPECT_EQ(rna.get_closed_regions(), expected_pairings);
}

TEST(RNAEntryClosedRegions, SimplePseudoKnot) {
    std::string sequence = "AAUU";
    std::string structure = "[(])";
    std::string structure2 = "[(])";
    knotergy::RNAEntry rna(name, sequence, structure);
    knotergy::RNAEntry rna2(name, sequence, structure2);
    std::vector<knotergy::ClosedRegion> expected_pairings = {{0, 3}};
    EXPECT_EQ(rna.get_closed_regions(), rna2.get_closed_regions());
    EXPECT_EQ(rna.get_closed_regions(), expected_pairings);
}

TEST(RNAEntryClosedRegions, PseudoknotWithNestedBP) {
    std::string sequence = "AAGGUAGUU";
    std::string structure = "[(..().])";
    knotergy::RNAEntry rna(name, sequence, structure);
    std::vector<knotergy::ClosedRegion> expected_pairings = {{4, 5}, {0, 8}};
    EXPECT_EQ(rna.get_closed_regions(), expected_pairings);
}

TEST(RNAEntryClosedRegions, PseudoknotWithNestedPseudoknot) {
    std::string sequence = "AAGGTTGGUUU";
    std::string structure = "[(..([)].])";
    knotergy::RNAEntry rna(name, sequence, structure);
    std::vector<knotergy::ClosedRegion> expected_pairings = {{4, 7}, {0, 10}};
    EXPECT_EQ(rna.get_closed_regions(), expected_pairings);
}

TEST(RNAEntryClosedRegions, ComplexPseudoKnot) {
    std::string sequence = "AGUAGUAUTTAAAAGGGGUUUUATTTTUGAUAATU";
    std::string structure = "([)([)()].(((([[[[))))(]]]]).().(])";
    knotergy::RNAEntry rna(name, sequence, structure);
    std::vector<knotergy::ClosedRegion> expected_pairings = {
        {6, 7}, {3, 8}, {10, 27}, {29, 30}, {0, 34}};
    EXPECT_EQ(rna.get_closed_regions(), expected_pairings);
}

TEST(RNAEntryClosedRegions, FauxPseudoKnots) {
    std::string sequence = "AAAUAUUAAUUU";
    std::string structure = "(([][])[()])";
    knotergy::RNAEntry rna(name, sequence, structure);
    std::vector<knotergy::ClosedRegion> expected_pairings = {
        {2, 3}, {4, 5}, {1, 6}, {8, 9}, {7, 10}, {0, 11}};
    EXPECT_EQ(rna.get_closed_regions(), expected_pairings);
}