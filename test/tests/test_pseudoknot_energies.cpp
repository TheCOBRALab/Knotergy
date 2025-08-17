#include <gtest/gtest.h>

#include <pipeline/shared.hpp>
#include <pipeline/input_pipeline.hpp>
#include <preprocessing/RNAEntry.hpp>
#include <preprocessing/ProcessedRNAEntry.hpp>
#include <preprocessing/RNAProcessor.hpp>
#include <loop_tree/LoopFactory.hpp>
#include <energy/ComputeEnergy.hpp>

#include <string>
#include <vector>


float pipeline(std::string sequence, std::string structure, bool round = false){
    knotergy::ViennaParams::load_energy_parameters();
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
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
    EXPECT_NEAR(result, -33.42, 0.001);
}

TEST(PseudoknottedEnergies, PseudoknotWithInband) {
    std::string sequence  = "GGGGGGGAGGGGGAAAACCCCCAGGGGGGGGACCCCCCCAAACCCCCCCC";
    std::string structure = "(((((((.(((((....))))).[[[[[[[[.)))))))...]]]]]]]]";
    float result = pipeline(sequence, structure);
    EXPECT_NEAR(result, -42.02, 0.009);
}

TEST(PseudoknottedEnergies, PseudoknotWithInbandRounded) {
    std::string sequence  = "GGGGGGGAGGGGGAAAACCCCCAGGGGGGGGACCCCCCCAAACCCCCCCC";
    std::string structure = "(((((((.(((((....))))).[[[[[[[[.)))))))...]]]]]]]]";
    bool round = true;
    std::cout << sequence << structure<<std::endl;
    float result = pipeline(sequence, structure, round);
    EXPECT_NEAR(result, -42.06, 0.001);
}

TEST(PseudoknottedEnergies, ExtendedPseudoknotWithInband) {
    std::string sequence  = "AUCCAUGCGAAGAACUAUGGAUCUCUGAAUGUUUUCGGUACAUUUCGGUGGUCCUUUAACGCCUUCCUUUGUGACACCAC";
    std::string structure = ".[[[[...[[[[......[[[[....((((((.......)))))).((((]]]]........]]]]...]].]]))))..";
    float result = pipeline(sequence, structure);
    EXPECT_NEAR(result, -15.08, 0.009);
}

TEST(PseudoknottedEnergies, ExtendedPseudoknotWithInbandRounded) {
    std::string sequence  = "AUCCAUGCGAAGAACUAUGGAUCUCUGAAUGUUUUCGGUACAUUUCGGUGGUCCUUUAACGCCUUCCUUUGUGACACCAC";
    std::string structure = ".[[[[...[[[[......[[[[....((((((.......)))))).((((]]]]........]]]]...]].]]))))..";
    bool round = true;
    float result = pipeline(sequence, structure, round);
    EXPECT_NEAR(result, -15.11, 0.001);
}
