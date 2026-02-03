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
float pipeline(std::string sequence, std::string structure, int dangle = 2,
               std::string param_file = "../../params/common/rna_turner2004.par",
               std::string pseudoknot_param_file = "../../params/pseudo/rna_pk_DirksPierce09_HotKnotsV2.json",
               std::string mod_param_path = "../../params/modified_bases") {

    knotergy::ViennaParams::load_energy_parameters(param_file, dangle, sequence);
    knotergy::PseudoknotParams::load_pk_param(pseudoknot_param_file);
    knotergy::RNAEntry rna(sequence, structure);
    std::vector<knotergy::modified_base_params> modified_params = knotergy::ViennaParams::load_modified_energy_parameters(mod_param_path);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna), modified_params));
    knotergy::LoopFactory factory(processed_rna);
    
    knotergy::ComputeEnergy energy(factory.get_root_node(), processed_rna, modified_params);

    return energy.getEnergy();
}

std::tuple<float, float, float> dangle_pipeline(
    std::string sequence,
    std::string structure,
    std::string param_file = "../../params/common/rna_turner2004.par",
    std::string pseudoknot_param_file = "../../params/pseudo/rna_pk_DirksPierce09_HotKnotsV2.json",
    std::string mod_param_path = "../../params/modified_bases"
) {
    float d0 = pipeline(sequence, structure, 0, param_file, pseudoknot_param_file, mod_param_path);
    float d1 = pipeline(sequence, structure, 1, param_file, pseudoknot_param_file, mod_param_path);
    float d2 = pipeline(sequence, structure, 2, param_file, pseudoknot_param_file, mod_param_path);

    return {d0, d1, d2};
}
}

// echo -e "AAAAAAAAAA6AAAAAAAAAAUUUUUUUUUUUUUUUUUUUUUUUUUUU" | RNAfold --mod-file=./params/modified_bases/rna_mod_m6A_parameters.json
TEST(mod_nonpk, small) {
    std::string sequence = "6U";
    std::string structure= "()";
    auto [d0, d1, d2] = dangle_pipeline(sequence, structure);

    EXPECT_NEAR(d0, 100000.0000, 0.001); // Turner 2004
    EXPECT_NEAR(d1, 100000.0000, 0.001); // Turner 2004
    EXPECT_NEAR(d2, 100000.0000, 0.001); // Turner 2004
}

// echo -e "AAAAAAAAAA6AAAAAAAAAAUUUUUUUUUUUUUUUUUUUUUUUUUUU" | RNAfold --mod-file=./params/modified_bases/rna_mod_m6A_parameters.json
TEST(mod_nonpk, stack_0_dangles) {
    std::string sequence = "AAAAAAAAAA6AAAAAAAAAAUUUUUUUUUUUUUUUUUUUUUUUUUUU";
    std::string structure= "(((((((((((((((((((((......)))))))))))))))))))))";
    auto [d0, d1, d2] = dangle_pipeline(sequence, structure);

    EXPECT_NEAR(d0, -13.25, 0.001); // Turner 2004
    EXPECT_NEAR(d1, -13.25, 0.001); // Turner 2004
    EXPECT_NEAR(d2, -13.25, 0.001); // Turner 2004
}

// // echo -e "GUUUUUAAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUUUUUUUU66AAAC" | RNAfold --mod-file=./params/modified_bases/rna_mod_m6A_parameters.json
// TEST(mod_nonpk, stack2_0_dangle) {
//     std::string sequence = "GUUUUUAAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUUUUUUUU66AAAC";
//     std::string structure= "(((((((((((((((((((((((((......)))))))))))))))))))))))))";
//     auto [d0, d1, d2] = dangle_pipeline(sequence, structure);

//     // EXPECT_NEAR(d0, -13.63, 0.001); // Turner 2004
//     // EXPECT_NEAR(d1, -14.23, 0.001); // Turner 2004
//     EXPECT_NEAR(d2, -18.55, 0.001); // Turner 2004
// }


// echo -e "AAAAAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUUPGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCC" | RNAfold --mod-file=./params/modified_bases/rna_mod_pseudouridine_parameters.json
// ./build/Knotergy -s "AAAAAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUUPGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCC" -r "(((((((((((((((((((....)))))))))))))))))))((((((((((((((((((((....))))))))))))))))))))" -p ./params/common/rna_turner2004.par
TEST(mod_nonpk, external1_0_dangle) {
    std::string sequence = "AAAAAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUUPGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCC";
    std::string structure= "(((((((((((((((((((....)))))))))))))))))))((((((((((((((((((((....))))))))))))))))))))";
    auto [d0, d1, d2] = dangle_pipeline(sequence, structure);

    EXPECT_NEAR(d0, -70.11, 0.001); // Turner 2004
    EXPECT_NEAR(d1, -70.11, 0.001); // Turner 2004
    EXPECT_NEAR(d2, -70.81, 0.001); // Turner 2004
}

// // echo -e "AAAAAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUUU6GGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCC\n(((((((((((((((((((....))))))))))))))))))).((((((((((((((((((((....))))))))))))))))))))" | RNAfold --mod-file=./params/modified_bases/rna_mod_m6A_parameters.json
// // ./build/Knotergy -p ./params/common/rna_turner2004.par -s "AAAAAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUUU6GGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCC" -r "(((((((((((((((((((....))))))))))))))))))).((((((((((((((((((((....))))))))))))))))))))"
// TEST(mod_nonpk, external2_1_dangle) {
//     std::string sequence = "AAAAAAAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUUUUU6GGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCC";
//     std::string structure= "(((((((((((((((((((....))))))))))))))))))).((((((((((((((((((((....))))))))))))))))))))";
//     auto [d0, d1, d2] = dangle_pipeline(sequence, structure);

//     EXPECT_NEAR(d0, -68.83, 0.001); // Turner 2004
//     EXPECT_NEAR(d1, -69.53, 0.001); // Turner 2004
//     EXPECT_NEAR(d2, -69.73, 0.001); // Turner 2004
// }

// echo -e "AAAAAUUUUUUUUPGGGGGGGGGCCCCCCCCC" | RNAfold --mod-file=./params/modified_bases/rna_mod_pseudouridine_parameters.json