#include <gtest/gtest.h>

#include <io/common.hpp>
#include <io/RNAInputManager.hpp>
#include <preprocessing/RNAEntry.hpp>
#include <preprocessing/ProcessedRNAEntry.hpp>
#include <preprocessing/RNAProcessor.hpp>
#include <loop_tree/LoopFactory.hpp>
#include <energy/ComputeEnergy.hpp>

#include <string>
#include <vector>

namespace {
float pipeline(std::string sequence, std::string structure, std::string param_file, bool round = false){
    int dangle = 2;
    knotergy::vrna_md_param vp = knotergy::ViennaParams::load_energy_parameters(param_file, dangle, sequence);
    knotergy::pk_param pkp = knotergy::PseudoknotParams::load_pk_param("../../params/pseudo/rna_pk_DirksPierce09_HotKnotsV2.json");
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    knotergy::LoopFactory factory(processed_rna);
    std::vector<knotergy::modified_base_param> mp;  // empty for unmodified bases
    knotergy::ComputeEnergy energy(factory.get_root_node(), processed_rna, vp, pkp, mp, round);

    return energy.getEnergy();
}

std::pair<float, float> get_turner_results(const std::string& sequence, const std::string& structure){
    bool round = true;
    float turner_result = pipeline(sequence, structure, "../../params/common/rna_turner2004.par");
    float turner_rounded = pipeline(sequence, structure, "../../params/common/rna_turner2004.par", round);
    return {turner_result, turner_rounded};
}
std::pair<float, float> get_dp_results(const std::string& sequence, const std::string& structure){
    bool round = true;
    float dp_result = pipeline(sequence, structure, "../../params/common/rna_DirksPierce09.par");
    float dp_rounded = pipeline(sequence, structure, "../../params/common/rna_DirksPierce09.par", round);
    return {dp_result, dp_rounded};
}

TEST(PK_energies, InfiniteEnergy_DP) {
    std::string sequence  = "GGCC";
    std::string structure = "[(])";
    float dp_result = pipeline(sequence, structure, "../../params/common/rna_DirksPierce09.par");

    EXPECT_NEAR(dp_result, 200003.5400, 0.009);
}

TEST(PK_energies, InfiniteEnergy_Turner) {
    std::string sequence  = "GGCC";
    std::string structure = "[(])";
    float turner_result = pipeline(sequence, structure, "../../params/common/rna_turner2004.par");

    EXPECT_NEAR(turner_result, 200003.5400, 0.009);
}

TEST(PK_energies, HType_DP) {
    std::string sequence  = "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    std::string structure = "[[[[[.......((((((((((........]]]]]......))))))))))";
    auto [dp_result, dp_rounded] = get_dp_results(sequence, structure);

    EXPECT_NEAR(dp_result,  -20.1912, 0.009);
    EXPECT_NEAR(dp_rounded, -20.1600, 0.001);
}

TEST(PK_energies, HType_Turner) {
    std::string sequence  = "GGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    std::string structure = "[[[[[.......((((((((((........]]]]]......))))))))))";
    auto [turner_result, turner_rounded] = get_turner_results(sequence, structure);

    EXPECT_NEAR(turner_result, -33.38, 0.009);
    EXPECT_NEAR(turner_rounded, -33.42, 0.009);
}

TEST(PK_energies, LLType_DP) {
    std::string sequence  = "GGGGGGGAGGGGGAAAACCCCCAGGGGGGGGACCCCCCCAAACCCCCCCC";
    std::string structure = "(((((((.(((((....))))).[[[[[[[[.)))))))...]]]]]]]]";
    auto [dp_result, dp_rounded] = get_dp_results(sequence, structure);

    EXPECT_NEAR(dp_result, -25.8912, 0.009);
    EXPECT_NEAR(dp_rounded, -25.86, 0.009);
}

TEST(PK_energies, LLType_Turner) {
    std::string sequence  = "GGGGGGGAGGGGGAAAACCCCCAGGGGGGGGACCCCCCCAAACCCCCCCC";
    std::string structure = "(((((((.(((((....))))).[[[[[[[[.)))))))...]]]]]]]]";
    auto [turner_result, turner_rounded] = get_turner_results(sequence, structure);

    EXPECT_NEAR(turner_result, -42.021, 0.009);
    EXPECT_NEAR(turner_rounded, -42.06, 0.001);
}

TEST(PK_energies, LLType_Long_DP) {
    std::string sequence  = "AUCCAUGCGAAGAACUAUGGAUCUCUGAAUGUUUUCGGUACAUUUCGGUGGUCCUUUAACGCCUUCCUUUGUGACACCAC";
    std::string structure = ".[[[[...[[[[......[[[[....((((((.......)))))).((((]]]]........]]]]...]].]]))))..";
    auto [dp_result, dp_rounded] = get_dp_results(sequence, structure);

    EXPECT_NEAR(dp_result, -8.9765, 0.009);
    EXPECT_NEAR(dp_rounded, -8.98, 0.001);
}


TEST(PK_energies, LLType_Long_Turner) {
    std::string sequence  = "AUCCAUGCGAAGAACUAUGGAUCUCUGAAUGUUUUCGGUACAUUUCGGUGGUCCUUUAACGCCUUCCUUUGUGACACCAC";
    std::string structure = ".[[[[...[[[[......[[[[....((((((.......)))))).((((]]]]........]]]]...]].]]))))..";
    auto [turner_result, turner_rounded] = get_turner_results(sequence, structure);

    EXPECT_NEAR(turner_result, -15.08, 0.009);
    EXPECT_NEAR(turner_rounded, -15.11, 0.001);
}

TEST(PK_energies, HLinType_DP) {
    std::string sequence  = "AAAAAAAGGGAAAGGGUUUGGGAAAGGGGGGGGGGGGUUUGGGUUUGGGUUUUGGGCCCCCCC";
    std::string structure = "(((((((...(((...)))...(((..[[[[[[[...)))...)))...))))...]]]]]]]";
    auto [dp_result, dp_rounded] = get_dp_results(sequence, structure);

    EXPECT_NEAR(dp_result, -0.7629, 0.009);
    EXPECT_NEAR(dp_rounded, -0.74, 0.001);
}


TEST(PK_energies, HLinType_Turner) {
    std::string sequence  = "AAAAAAAGGGAAAGGGUUUGGGAAAGGGGGGGGGGGGUUUGGGUUUGGGUUUUGGGCCCCCCC";
    std::string structure = "(((((((...(((...)))...(((..[[[[[[[...)))...)))...))))...]]]]]]]";
    auto [turner_result, turner_rounded] = get_turner_results(sequence, structure);

    EXPECT_NEAR(turner_result, -5.831, 0.009);
    EXPECT_NEAR(turner_rounded, -5.84, 0.001);
}


TEST(PK_energies, multi_nested_DP) {
    std::string sequence  = "AAAGGAAAGGGUUUGGGGGGGGGGGAAAGGGUUUGGGGGGGGGGGUUUGGGGGCCCCCCCCCCCCCCCCCC";
    std::string structure = "(((..(((...)))...[[[[[...(((...)))...[[[[[...)))......]]]]]]]]]].......";
    auto [dp_result, dp_rounded] = get_dp_results(sequence, structure);

    EXPECT_NEAR(dp_result, 0.047, 0.009);
    EXPECT_NEAR(dp_rounded, 0.07, 0.001);
}

TEST(PK_energies, multi_nested_Turner) {
    std::string sequence  = "AAAGGAAAGGGUUUGGGGGGGGGGGAAAGGGUUUGGGGGGGGGGGUUUGGGGGCCCCCCCCCCCCCCCCCC";
    std::string structure = "(((..(((...)))...[[[[[...(((...)))...[[[[[...)))......]]]]]]]]]].......";
    auto [turner_result, turner_rounded] = get_turner_results(sequence, structure);

    EXPECT_NEAR(turner_result, -5.748, 0.009);
    EXPECT_NEAR(turner_rounded, -5.77, 0.009);
}

TEST(PK_energies, two_multiloops_DP) {
    std::string sequence  = "AAAAGGAAAGGGGUUUGGGAAAGGGAAAGGGUUUGGGAAAGGGGGGGGGGGGUUUGGGUUUGGGUUUUGGGCCCCCCC";
    std::string structure = "((((..(((....)))...(((...(((...)))...(((..[[[[[[[...)))...)))...))))...]]]]]]]";
    auto [dp_result, dp_rounded] = get_dp_results(sequence, structure);
    
    EXPECT_NEAR(dp_result, 4.6223, 0.009);
    EXPECT_NEAR(dp_rounded, 4.65, 0.001);
}

TEST(PK_energies, two_multiloops_Turner) {
    std::string sequence  = "AAAAGGAAAGGGGUUUGGGAAAGGGAAAGGGUUUGGGAAAGGGGGGGGGGGGUUUGGGUUUGGGUUUUGGGCCCCCCC";
    std::string structure = "((((..(((....)))...(((...(((...)))...(((..[[[[[[[...)))...)))...))))...]]]]]]]";
    auto [turner_result, turner_rounded] = get_turner_results(sequence, structure);

    EXPECT_NEAR(turner_result, -0.189, 0.009);
    EXPECT_NEAR(turner_rounded, -0.2, 0.001);
}


TEST(PK_energies, kissing_muliloop_DP) {
    std::string sequence  = "AAAAGGGAAAGGGUUUGGGAAAAGGGGGGGGGGGUUUUGGUUUUGGGGGGGGGGGGGGGGGAAAAGGGAAAACCCCCCCCCCCUUUUGGGGGAAAGGGUUUGGUUUU";
    std::string structure = "((((...(((...)))...(((([[[[[[[[[[[))))..)))).................((((...((((]]]]]]]]]]])))).....(((...)))..))))";
    auto [dp_result, dp_rounded] = get_dp_results(sequence, structure);

    EXPECT_NEAR(dp_result, -1.9068, 0.009);
    EXPECT_NEAR(dp_rounded, -1.8600, 0.001);
}

TEST(PK_energies, kissing_muliloop_Turner) {
    std::string sequence  = "AAAAGGGAAAGGGUUUGGGAAAAGGGGGGGGGGGUUUUGGUUUUGGGGGGGGGGGGGGGGGAAAAGGGAAAACCCCCCCCCCCUUUUGGGGGAAAGGGUUUGGUUUU";
    std::string structure = "((((...(((...)))...(((([[[[[[[[[[[))))..)))).................((((...((((]]]]]]]]]]])))).....(((...)))..))))";
    auto [turner_result, turner_rounded] = get_turner_results(sequence, structure);

    EXPECT_NEAR(turner_result, -11.422, 0.009);
    EXPECT_NEAR(turner_rounded, -11.44, 0.001);
}

TEST(PK_energies, nestedPK_in_band_DP) {
    std::string sequence  = "AAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGUUUUUUUUUUUCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAGGGGGGGGGGGGGGGGUUUUUUUUUUUUUUUUUUUUUCCCCCCCCCCCCCCCC";
    std::string structure = "((((((((((((((((((((([[[[[[[[[[[[[[[[[[[[)))))))))))]]]]]]]]]]]]]]]]]]]]((((((((((([[[[[[[[[[[[[[[[)))))))))))))))))))))]]]]]]]]]]]]]]]]";
    auto [dp_result, dp_rounded] = get_dp_results(sequence, structure);

    EXPECT_NEAR(dp_result, -60.0667, 0.009);
    EXPECT_NEAR(dp_rounded, -59.93, 0.009);
}

TEST(PK_energies, nestedPK_in_band_Turner) {
    std::string sequence  = "AAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGUUUUUUUUUUUCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAGGGGGGGGGGGGGGGGUUUUUUUUUUUUUUUUUUUUUCCCCCCCCCCCCCCCC";
    std::string structure = "((((((((((((((((((((([[[[[[[[[[[[[[[[[[[[)))))))))))]]]]]]]]]]]]]]]]]]]]((((((((((([[[[[[[[[[[[[[[[)))))))))))))))))))))]]]]]]]]]]]]]]]]";
    auto [turner_result, turner_rounded] = get_turner_results(sequence, structure);

    EXPECT_NEAR(turner_result, -99.4670, 0.009);
    EXPECT_NEAR(turner_rounded, -99.54, 0.009);
}

}
