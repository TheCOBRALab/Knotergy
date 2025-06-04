#include <gtest/gtest.h>
#include <rna_regions/RNAEntry.hpp>
#include <string>
#include <vector>


TEST(RNAEntryBasics, populate){
    std::string name = "Test Input";
    std::string sequence = "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    std::string structure = "(((((.........................)))))................";
    compute_energy::RNAEntry rna(name, sequence, structure);
    
    EXPECT_EQ(rna.get_name(), "Test Input");
    EXPECT_EQ(rna.get_sequence(), "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC");
    EXPECT_EQ(rna.get_structure(), "(((((.........................)))))................");
}

TEST(RNAEntryBasics, GettersAndSetters){
    std::string name = "Test Input";
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
    std::string name = "Test Input";
    std::string sequence = "AAUU";
    std::string structure = "(..)";
    compute_energy::RNAEntry rna(name, sequence, structure);
    std::vector<int> expected_pairings = {3,-1,-1,0};
    
    EXPECT_EQ(rna.get_pairings(), expected_pairings);
}