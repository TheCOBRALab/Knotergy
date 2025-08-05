#include <gtest/gtest.h>

#include <pipeline/shared.hpp>
#include <pipeline/input_pipeline.hpp>
#include <preprocessing/RNAEntry.hpp>
#include <preprocessing/RNAProcessedEntry.hpp>
#include <preprocessing/RNAProcessor.hpp>
#include <loop_tree/LoopFactory.hpp>
#include <energy/ComputeEnergy.hpp>

#include <string>
#include <vector>


float pipeline(std::string sequence, std::string structure, bool round = false){
    knotergy::load_energy_parameters();
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::RNAProcessedEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    knotergy::LoopFactory factory(processed_rna);
    knotergy::ComputeEnergy energy(factory.get_root_node(), sequence,processed_rna, round);

    return energy.getEnergy();
}

TEST(PseudoknottedEnergies, SimplePseudoknot) {
    std::string sequence  = "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    std::string structure = "[[[[[.......((((((((((........]]]]]......))))))))))";
    float result = pipeline(sequence, structure);
    EXPECT_NEAR(result, -33.38, 0.009);
}

TEST(PseudoknottedEnergies, SimplePseudoknotRounded) {
    std::string sequence  = "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    std::string structure = "[[[[[.......((((((((((........]]]]]......))))))))))";
    bool round = true;
    float result = pipeline(sequence, structure, round);
    EXPECT_NEAR(result, -33.29, 0.001);
}
