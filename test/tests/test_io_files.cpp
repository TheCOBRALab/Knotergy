#include "test_utils.hpp"

namespace {

double pipeline(std::string input_file, std::string param_file,
                std::string pseudoknot_param_file = pkp_file) {
    knotergy::vrna_md_param vp = knotergy::ViennaParams::load_energy_parameters(param_file);
    knotergy::pk_param pkp = knotergy::PseudoknotParams::load_pk_param(pseudoknot_param_file);
    knotergy::RNAEntry rna = knotergy::RNAInputManager::get_all_inputs(input_file, "", "").front();
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    knotergy::LoopFactory factory(processed_rna, vp);
    knotergy::all_mod_params mp;  // empty for no modified bases
    knotergy::ComputeEnergy energy(factory.get_root_node(), processed_rna, vp, pkp, mp);

    return energy.getEnergy();
}

TEST(IO_InputFile, tree_construction_overflow) {
    const std::string input_file =
        "../../test/tests/input_files/loop_tree_destructor_overflow_test.txt";

    double turner_energy = pipeline(input_file, turner_file);

    EXPECT_NEAR(turner_energy, -22492.70, 0.0005);  // Turner 2004
}

}  // namespace