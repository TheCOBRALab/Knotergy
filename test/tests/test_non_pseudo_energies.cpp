#include <gtest/gtest.h>

#include <helpers/common.hpp>
#include <helpers/main_helpers.hpp>
#include <rna_regions/RNAEntry.hpp>
#include <loops/LoopFactory.hpp>
#include <energy/ComputeEnergy.hpp>

#include <string>
#include <vector>


// echo -e "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC\n(((((.........................)))))................" | RNAeval -P ./rna_langdon2018.par
TEST(NonPseudoKnottedEnergies, BasicStack) {
    std::string sequence = "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    std::string structure = "(((((.........................)))))................";
    knotergy::load_energy_parameters();
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::LoopFactory factory(rna);
    knotergy::ComputeEnergy energy(factory.get_root_node(), sequence);
    
    EXPECT_NEAR(energy.getEnergy(), -11.30, 0.001);
}

// echo -e "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCC\n...........((((((....((((((............))))))....((((((......))))))...((....)).))))))" | RNAeval -P ./rna_langdon2018.par
TEST(NonPseudoKnottedEnergies, MultiLoopWithStacks) {
    std::string sequence = "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCC";
    std::string structure = "...........((((((....((((((............))))))....((((((......))))))...((....)).))))))";
    knotergy::load_energy_parameters();
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::LoopFactory factory(rna);
    knotergy::ComputeEnergy energy(factory.get_root_node(), sequence);

    EXPECT_NEAR(energy.getEnergy(), -44.00, 0.001);
}