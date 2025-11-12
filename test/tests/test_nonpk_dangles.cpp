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
#include <tuple>

namespace {
float get_energy(std::string sequence, std::string structure, int dangle, std::string param_file = "../../params/common/rna_turner2004.par") {
    knotergy::ViennaParams::load_energy_parameters(param_file);
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    knotergy::LoopFactory factory(processed_rna);
    knotergy::ComputeEnergy energy(factory.get_root_node(), sequence, processed_rna, false, dangle);

    return energy.getEnergy();
}

std::tuple<float, float, float> pipeline(
    std::string sequence,
    std::string structure,
    std::string param_file = "../../params/common/rna_turner2004.par"
) {
    float d0 = get_energy(sequence, structure, 0, param_file);
    float d1 = get_energy(sequence, structure, 1, param_file);
    float d2 = get_energy(sequence, structure, 2, param_file);

    return {d0, d1, d2};
}
} // namespace




TEST(Dangles, D0_external_simple) {
    std::string sequence  = "AAAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUU";
    std::string structure = ".(((((.........................))))).";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 3.90, 0.009);
    EXPECT_NEAR(d1, 3.2, 0.009);
    EXPECT_NEAR(d2, 3.2, 0.009);
}

TEST(Dangles, external_simple_NoEnds) {
    std::string sequence  = "AAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUU";
    std::string structure = "(((((.........................)))))";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 3.90, 0.009);
    EXPECT_NEAR(d1, 3.90, 0.009);
    EXPECT_NEAR(d2, 3.90, 0.009);
}

TEST(Dangles, external_simple_adjacent) {
    std::string sequence  = "AAAAAAUUUUUGGGGGCCCCCCAAAAAUUUUUUAAAAUUUUGGGCCAAAAAAAUUUUU";
    std::string structure = ".((((((((((((((...)))))))))))))).((((((((........)))))))).";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, -16.40, 0.009);
    EXPECT_NEAR(d1, -17.20, 0.009);
    EXPECT_NEAR(d2, -17.60, 0.009);
}

TEST(Dangles, multiloop_left_dangle) {
    std::string sequence  = "AAAAAAAAAUUUUUUUAAAAAUUUUUUU";
    std::string structure = "(((.(((....)))..((...))..)))";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 15.10, 0.009);
    EXPECT_NEAR(d1, 13.5, 0.009);
    EXPECT_NEAR(d2, 12.9, 0.009);
}

TEST(Dangles, multiloop_left_dangle_chained) {
    std::string sequence  = "AAAAAAAAAUUUUUUAAAAAUUUUUUU";
    std::string structure = "(((.(((....))).((...))..)))";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 15.10, 0.009);
    EXPECT_NEAR(d1, 13.6, 0.009);
    EXPECT_NEAR(d2, 12.9, 0.009);
}

TEST(Dangles, multiloop_left_touch) {
    std::string sequence  = "AAAAAAAAUUUUUUAAAAAUUUUUUU";
    std::string structure = "((((((....))).((...))..)))";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 15.10, 0.009);
    EXPECT_NEAR(d1, 14.4, 0.009);
    EXPECT_NEAR(d2, 12.9, 0.009);
}


TEST(Dangles, multiloop_right_dangle) {
    std::string sequence  = "AAAAAAAAAUUUUUUAAAAAUUUUUUU";
    std::string structure = "(((..(((....))).((...)).)))";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 15.10, 0.009);
    EXPECT_NEAR(d1, 13.1, 0.009);
    EXPECT_NEAR(d2, 12.4, 0.009);
}

TEST(Dangles, multiloop_right_touch) {
    std::string sequence  = "AAAAAAAAAUUUUUUAAAAAUUUUUU";
    std::string structure = "(((..(((....))).((...)))))";
    auto [d0, d1, d2] = pipeline(sequence, structure);

    EXPECT_NEAR(d0, 15.10, 0.009);
    EXPECT_NEAR(d1, 13.3, 0.009);
    EXPECT_NEAR(d2, 12.4, 0.009);
}

// // echo -e "GGUUUUUUUUAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUGGGGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCC\n.(((((((((((((((..........))))))((((((...)))))))))))))).((((((((........)))))))))." | RNAeval -d 
// TEST(Dangles, external_multiloop) {
//     std::string sequence  = "GGUUUUUUUUAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUGGGGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCC";
//     std::string structure = ".(((((((((((((((..........))))))((((((...)))))))))))))).((((((((........))))))))).";
//     auto [d0, d1, d2] = pipeline(sequence, structure);

//     EXPECT_NEAR(d0, -2.50, 0.009);
//     EXPECT_NEAR(d1, -4.30, 0.009);
//     EXPECT_NEAR(d2, -8.90, 0.009);
// }