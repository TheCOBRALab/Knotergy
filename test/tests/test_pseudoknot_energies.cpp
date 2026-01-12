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

namespace {
float pipeline(std::string sequence, std::string structure, std::string param_file = "../../params/common/rna_turner2004.par", bool round = false){
    knotergy::ViennaParams::load_energy_parameters(param_file);
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    knotergy::LoopFactory factory(processed_rna);
    std::vector<knotergy::modified_base_params> mod_params;  // empty for unmodified bases
    int dangle = 2;
    knotergy::ComputeEnergy energy(factory.get_root_node(), processed_rna, mod_params, dangle, round);

    return energy.getEnergy();
}
}
TEST(PseudoknottedEnergies, SimplePseudoknot) {
    std::string sequence  = "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    std::string structure = "[[[[[.......((((((((((........]]]]]......))))))))))";
    float result = pipeline(sequence, structure);
    float dp_result = pipeline(sequence, structure, "../../params/common/rna_DirksPierce09.par");

    EXPECT_NEAR(result, -33.38, 0.009);
    EXPECT_NEAR(dp_result, -20.1912, 0.009);
}

TEST(PseudoknottedEnergies, SimplePseudoknotRounded) {
    std::string sequence  = "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    std::string structure = "[[[[[.......((((((((((........]]]]]......))))))))))";
    bool round = true;
    float result = pipeline(sequence, structure, "", round);
    float dp_result = pipeline(sequence, structure, "../../params/common/rna_DirksPierce09.par", round);

    EXPECT_NEAR(result, -33.42, 0.001);
    EXPECT_NEAR(dp_result, -20.16, 0.001);
}

TEST(PseudoknottedEnergies, PseudoknotWithInband) {
    std::string sequence  = "GGGGGGGAGGGGGAAAACCCCCAGGGGGGGGACCCCCCCAAACCCCCCCC";
    std::string structure = "(((((((.(((((....))))).[[[[[[[[.)))))))...]]]]]]]]";
    float result = pipeline(sequence, structure);
    float dp_result = pipeline(sequence, structure, "../../params/common/rna_DirksPierce09.par");

    EXPECT_NEAR(result, -42.021, 0.009);
    EXPECT_NEAR(dp_result, -25.8912, 0.009);
}

TEST(PseudoknottedEnergies, PseudoknotWithInbandRounded) {
    std::string sequence  = "GGGGGGGAGGGGGAAAACCCCCAGGGGGGGGACCCCCCCAAACCCCCCCC";
    std::string structure = "(((((((.(((((....))))).[[[[[[[[.)))))))...]]]]]]]]";
    bool round = true;
    
    float result = pipeline(sequence, structure, "", round);
    float dp_result = pipeline(sequence, structure, "../../params/common/rna_DirksPierce09.par", round);

    EXPECT_NEAR(result, -42.06, 0.001);
    EXPECT_NEAR(dp_result, -25.86, 0.001);
}

TEST(PseudoknottedEnergies, ExtendedPseudoknotWithInband) {
    std::string sequence  = "AUCCAUGCGAAGAACUAUGGAUCUCUGAAUGUUUUCGGUACAUUUCGGUGGUCCUUUAACGCCUUCCUUUGUGACACCAC";
    std::string structure = ".[[[[...[[[[......[[[[....((((((.......)))))).((((]]]]........]]]]...]].]]))))..";
    float result = pipeline(sequence, structure);
    float dp_result = pipeline(sequence, structure, "../../params/common/rna_DirksPierce09.par");

    EXPECT_NEAR(result, -15.08, 0.009);
    EXPECT_NEAR(dp_result, -8.9765, 0.009);
}

TEST(PseudoknottedEnergies, ExtendedPseudoknotWithInbandRounded) {
    std::string sequence  = "AUCCAUGCGAAGAACUAUGGAUCUCUGAAUGUUUUCGGUACAUUUCGGUGGUCCUUUAACGCCUUCCUUUGUGACACCAC";
    std::string structure = ".[[[[...[[[[......[[[[....((((((.......)))))).((((]]]]........]]]]...]].]]))))..";
    bool round = true;
    float result = pipeline(sequence, structure, "", round);
    float dp_result = pipeline(sequence, structure, "../../params/common/rna_DirksPierce09.par", round);

    EXPECT_NEAR(result, -15.11, 0.001);
    EXPECT_NEAR(dp_result, -8.98, 0.001);
}

TEST(PseudoknottedEnergies, MultiLoopThatSpansABand) {
    std::string sequence  = "AAAAAAAGGGAAAGGGUUUGGGAAAGGGGGGGGGGGGUUUGGGUUUGGGUUUUGGGCCCCCCC";
    std::string structure = "(((((((...(((...)))...(((..[[[[[[[...)))...)))...))))...]]]]]]]";
    float result = pipeline(sequence, structure);
    float dp_result = pipeline(sequence, structure, "../../params/common/rna_DirksPierce09.par");

    EXPECT_NEAR(result, -5.831, 0.009);
    EXPECT_NEAR(dp_result, -0.7629, 0.009);
}

TEST(PseudoknottedEnergies, MultiLoopThatSpansABandRounded) {
    std::string sequence  = "AAAAAAAGGGAAAGGGUUUGGGAAAGGGGGGGGGGGGUUUGGGUUUGGGUUUUGGGCCCCCCC";
    std::string structure = "(((((((...(((...)))...(((..[[[[[[[...)))...)))...))))...]]]]]]]";
    bool round = true;
    float result = pipeline(sequence, structure, "", round);
    float dp_result = pipeline(sequence, structure, "../../params/common/rna_DirksPierce09.par", round);

    EXPECT_NEAR(result, -5.84, 0.001);
    EXPECT_NEAR(dp_result, -0.74, 0.001);
}

TEST(PseudoknottedEnergies, MultiLoopThatSpansABand2) {
    std::string sequence  = "AAAGGAAAGGGUUUGGGGGGGGGGGAAAGGGUUUGGGGGGGGGGGUUUGGGGGCCCCCCCCCCCCCCCCCC";
    std::string structure = "(((..(((...)))...[[[[[...(((...)))...[[[[[...)))......]]]]]]]]]].......";
    float result = pipeline(sequence, structure);
    float dp_result = pipeline(sequence, structure, "../../params/common/rna_DirksPierce09.par");

    EXPECT_NEAR(result, -5.748, 0.009);
    EXPECT_NEAR(dp_result, 0.047, 0.009);
}

TEST(PseudoknottedEnergies, MultiLoopThatSpansABandRounded2) {
    std::string sequence  = "AAAGGAAAGGGUUUGGGGGGGGGGGAAAGGGUUUGGGGGGGGGGGUUUGGGGGCCCCCCCCCCCCCCCCCC";
    std::string structure = "(((..(((...)))...[[[[[...(((...)))...[[[[[...)))......]]]]]]]]]].......";
    bool round = true;
    float result = pipeline(sequence, structure, "", round);
    float dp_result = pipeline(sequence, structure, "../../params/common/rna_DirksPierce09.par", round);

    EXPECT_NEAR(result, -5.77, 0.001);
    EXPECT_NEAR(dp_result, 0.07, 0.001);
}

TEST(PseudoknottedEnergies, TwoMultiLoopsThatSpansABand) {
    std::string sequence  = "AAAAGGAAAGGGGUUUGGGAAAGGGAAAGGGUUUGGGAAAGGGGGGGGGGGGUUUGGGUUUGGGUUUUGGGCCCCCCC";
    std::string structure = "((((..(((....)))...(((...(((...)))...(((..[[[[[[[...)))...)))...))))...]]]]]]]";
    float result = pipeline(sequence, structure);
    float dp_result = pipeline(sequence, structure, "../../params/common/rna_DirksPierce09.par");
    
    EXPECT_NEAR(result, -0.189, 0.009);
    EXPECT_NEAR(dp_result, 4.6223, 0.009);
}

TEST(PseudoknottedEnergies, TwoMultiLoopsThatSpansABandRounded) {
    std::string sequence  = "AAAAGGAAAGGGGUUUGGGAAAGGGAAAGGGUUUGGGAAAGGGGGGGGGGGGUUUGGGUUUGGGUUUUGGGCCCCCCC";
    std::string structure = "((((..(((....)))...(((...(((...)))...(((..[[[[[[[...)))...)))...))))...]]]]]]]";
    bool round = true;
    float result = pipeline(sequence, structure, "", round);
    float dp_result = pipeline(sequence, structure, "../../params/common/rna_DirksPierce09.par", round);

    EXPECT_NEAR(result, -0.2, 0.001);
    EXPECT_NEAR(dp_result, 4.65, 0.001);
}

TEST(PseudoknottedEnergies, TwoMultiLoopsThatSpansABand2) {
    std::string sequence  = "AAAAGGGAAAGGGUUUGGGAAAAGGGGGGGGGGGUUUUGGUUUUGGGGGGGGGGGGGGGGGAAAAGGGAAAACCCCCCCCCCCUUUUGGGGGAAAGGGUUUGGUUUU";
    std::string structure = "((((...(((...)))...(((([[[[[[[[[[[))))..)))).................((((...((((]]]]]]]]]]])))).....(((...)))..))))";
    float result = pipeline(sequence, structure);
    float dp_result = pipeline(sequence, structure, "../../params/common/rna_DirksPierce09.par");

    EXPECT_NEAR(result, -11.422, 0.009);
    EXPECT_NEAR(dp_result, -1.9068, 0.009);
}

TEST(PseudoknottedEnergies, TwoMultiLoopsThatSpansABandRounded2) {
    std::string sequence  = "AAAAGGGAAAGGGUUUGGGAAAAGGGGGGGGGGGUUUUGGUUUUGGGGGGGGGGGGGGGGGAAAAGGGAAAACCCCCCCCCCCUUUUGGGGGAAAGGGUUUGGUUUU";
    std::string structure = "((((...(((...)))...(((([[[[[[[[[[[))))..)))).................((((...((((]]]]]]]]]]])))).....(((...)))..))))";
    bool round = true;
    float result = pipeline(sequence, structure, "", round);
    float dp_result = pipeline(sequence, structure, "../../params/common/rna_DirksPierce09.par", round);

    EXPECT_NEAR(result, -11.44, 0.001);
    EXPECT_NEAR(dp_result, -1.8600, 0.001);
}