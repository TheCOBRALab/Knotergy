#include <gtest/gtest.h>
#include <rna_regions/RNAEntry.hpp>
#include <string>
#include <vector>
#include <helpers/common.hpp>

std::string name = "Test Input";

TEST(RNAEntryBasics, populate){
    std::string sequence = "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    std::string structure = "(((((.........................)))))................";
    compute_energy::RNAEntry rna(name, sequence, structure);
    
    EXPECT_EQ(rna.get_name(), "Test Input");
    EXPECT_EQ(rna.get_sequence(), "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC");
    EXPECT_EQ(rna.get_structure(), "(((((.........................)))))................");
}

TEST(RNAEntryBasics, GettersAndSetters){
    std::string sequence = "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    std::string structure = "(((((.........................)))))................";
    compute_energy::RNAEntry rna;

    rna.set_name(name);
    rna.set_sequence(sequence);
    rna.set_structure(structure);

    EXPECT_EQ(rna.get_name(), "Test Input");
    EXPECT_EQ(rna.get_sequence(), "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC");
    EXPECT_EQ(rna.get_structure(), "(((((.........................)))))................");
}

TEST(RNAEntryPairings, SimplePairing){
    std::string sequence = "AAUU";
    std::string structure = "(..)";
    compute_energy::RNAEntry rna(name, sequence, structure);
    std::vector<size_t> expected_pairings = {3,compute_energy::NULL_INDEX, compute_energy::NULL_INDEX,0};
    
    EXPECT_EQ(rna.get_pairings(), expected_pairings);
}

TEST(RNAEntryClosedRegions, SimpleRegion){
    std::string sequence = "AAUU";
    std::string structure = "(..)";
    compute_energy::RNAEntry rna(name, sequence, structure);
    std::vector<compute_energy::Region> expected_pairings;
    expected_pairings.emplace_back(0, 3);
    EXPECT_EQ(rna.get_closed_regions(), expected_pairings);
}

TEST(RNAEntryClosedRegions, SimplePseudoKnot){
    std::string sequence = "AAUU";
    std::string structure = "[(])";
    std::string structure2 = "[(])";
    compute_energy::RNAEntry rna(name, sequence, structure);
    compute_energy::RNAEntry rna2(name, sequence, structure2);
    std::vector<compute_energy::Region> expected_pairings;
    expected_pairings.emplace_back(0, 3);

    EXPECT_EQ(rna.get_closed_regions(), rna2.get_closed_regions());
    EXPECT_EQ(rna.get_closed_regions(), expected_pairings);
}


TEST(RNAEntryClosedRegions, PseudoknotWithNestedBP){
    std::string sequence = "AAGGUAGUU";
    std::string structure = "[(..().])";
    compute_energy::RNAEntry rna(name, sequence, structure);
    std::vector<compute_energy::Region> expected_pairings;
    expected_pairings.emplace_back(4, 5);
    expected_pairings.emplace_back(0, 8);
    EXPECT_EQ(rna.get_closed_regions(), expected_pairings);
}

TEST(RNAEntryClosedRegions, PseudoknotWithNestedPseudoknot){
    std::string sequence = "AAGGTTGGUUU";
    std::string structure = "[(..([)].])";
    compute_energy::RNAEntry rna(name, sequence, structure);
    std::vector<compute_energy::Region> expected_pairings;
    expected_pairings.emplace_back(4, 7);
    expected_pairings.emplace_back(0, 10);
    EXPECT_EQ(rna.get_closed_regions(), expected_pairings);
}

TEST(RNAEntryClosedRegions, ComplexPseudoKnot){
    std::string sequence = "AGUAGUAUTTAAAAGGGGUUUUATTTTUGAUAATU";
    std::string structure = "([)([)()].(((([[[[))))(]]]]).().(])";
    compute_energy::RNAEntry rna(name, sequence, structure);
    std::vector<compute_energy::Region> expected_pairings;
    expected_pairings.emplace_back(6, 7);
    expected_pairings.emplace_back(3, 8);
    expected_pairings.emplace_back(10, 27);
    expected_pairings.emplace_back(29, 30);
    expected_pairings.emplace_back(0, 34);
    EXPECT_EQ(rna.get_closed_regions(), expected_pairings);
}

TEST(RNAEntryClosedRegions, FauxPseudoKnots){
    std::string sequence = "AAAUAUUAAUUU";
    std::string structure = "(([][])[()])";
    compute_energy::RNAEntry rna(name, sequence, structure);
    std::vector<compute_energy::Region> expected_pairings;
    expected_pairings.emplace_back(2, 3);
    expected_pairings.emplace_back(4, 5);
    expected_pairings.emplace_back(1, 6);
    expected_pairings.emplace_back(8, 9);
    expected_pairings.emplace_back(7, 10);
    expected_pairings.emplace_back(0, 11);
    EXPECT_EQ(rna.get_closed_regions(), expected_pairings);
}