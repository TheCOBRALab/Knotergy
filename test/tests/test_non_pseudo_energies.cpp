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

float pipeline(std::string sequence, std::string structure){
    knotergy::load_energy_parameters();
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    knotergy::LoopFactory factory(processed_rna);
    knotergy::ComputeEnergy energy(factory.get_root_node(), sequence, processed_rna);

    return energy.getEnergy();
}

// echo -e "AGGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC\n.(((((.........................)))))................" | RNAeval
TEST(NonPseudoKnottedEnergies, BasicStack) {
    std::string sequence = "AGGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    std::string structure = ".(((((.........................)))))................";
    float energy = pipeline(sequence, structure);

    EXPECT_NEAR(energy, -8.30, 0.001);
}

// echo -e "GGGGAGAAAAAAAAAUUUUUU\n((((((.........))))))" | RNAeval
TEST(NonPseudoKnottedEnergies, StacksFullStructure) {
    std::string sequence = "GGGGAGAAAAAAAAAUUUUUU";
    std::string structure = "((((((.........))))))";
    float energy = pipeline(sequence, structure);

    EXPECT_NEAR(energy, 3.70, 0.001);
}

// echo -e "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCC\n((((((.....................))))))................((((((......))))))...(((...)))......" | RNAeval
TEST(NonPseudoKnottedEnergies, MultipleStacks) {
    std::string sequence = "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCC";
    std::string structure = "((((((.....................))))))................((((((......))))))...(((...)))......";
    float energy = pipeline(sequence, structure);

    EXPECT_NEAR(energy, -28.10, 0.001);
}


// echo -e "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCCA\n...........((((((....((((((............))))))....((((((......))))))...((....)).))))))." | RNAeval
TEST(NonPseudoKnottedEnergies, MultiLoopWithStacks) {
    std::string sequence = "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCCA";
    std::string structure = "...........((((((....((((((............))))))....((((((......))))))...((....)).)))))).";
    float energy = pipeline(sequence, structure);

    EXPECT_NEAR(energy, -41.10, 0.001);
}

// echo -e "AAAAAAGGGGCCCCCCAAAAAAAAAAAAGGGGGGAAAAGGGGGGUUUUUUCCCCCCAAAUUUUUUAAGUUUUUU\n((((((....((((((............))))))....((((((......))))))...((....)).))))))" | RNAeval
TEST(NonPseudoKnottedEnergies, MultiLoopWithStacksFull) {
    std::string sequence = "AAAAAAGGGGCCCCCCAAAAAAAAAAAAGGGGGGAAAAGGGGGGUUUUUUCCCCCCAAAUUUUUUAAGUUUUUU";
    std::string structure = "((((((....((((((............))))))....((((((......))))))...((....)).))))))";
    float energy = pipeline(sequence, structure);

    EXPECT_NEAR(energy, -21.9, 0.001);
}

// echo -e "GGGGUUAUUUUAUUAAAAAUAACCCUGGUUUUUAAGGCGGGGUCGUGCGGUAAGGGAACCC\n((((..(...).((...)))..(((.(.((...))..(...).)..(...)..)))..)))" | RNAeval
TEST(NonPseudoKnottedEnergies, MultiWithMultiLoop) {
    std::string sequence = "GGGGUUAUUUUAUUAAAAAUAACCCUGGUUUUUAAGGCGGGGUCGUGCGGUAAGGGAACCC";
    std::string structure = "((((..(...).((...)))..(((.(.((...))..(...).)..(...)..)))..)))";
    float energy = pipeline(sequence, structure);

    EXPECT_NEAR(energy, 30.90, 0.001);
}

TEST(NonPseudoKnottedEnergies, PositiveInternalLoops) {
    std::string sequence  = "GGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAAUUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGCGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACG";
    std::string structure = "(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))..........((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..........))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((....)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))";
    float energy = pipeline(sequence, structure);
    EXPECT_NEAR(energy, -471.90, 0.001);
} 

TEST(NonPseudoKnottedEnergies, MoreLongSequences) {
    std::string sequence  = "UUAAAAGGGAUGCCUCUCCUGUUCAUCUUGUGGAGAAGCAUUCGAUAAGGUCAUCAUAAUGGGUCCAGCUUUGCGACCUGGCGAGAUUAGUCAGGAAAAUGUGAAGUGGGUCUUCGCUUUCCAGGUACAGGAGGCUCGCCCCGCUCAUCCAGUUCGUCCCCUAACCACUUGUUUUCUCAGGAUAGUUUGUUUUGUACACCCGUGUACAUACACAUGUAUCACACCCCAGAUUGCCGAAUGUUUCGUUCGGUCGAGCCUGACUAUGCAUAAACCCUACCUCUGAAACCUUGGGCAACUCACUACUUCCGAGCUAAAUCCCUCUGUUUGAGCUAGCCUGAGAUUUCAACUGGCUUCGGCCUUGUUUAUACCAUCGUUUGCUGAUCCAUUGAAGAAAUAAGUUACCGAUGGCCCCAAACUGACGAUCACUAUUCUUUCCAUAGGAGUUAUGGGUAUACUGCCCGUAGACGGAAAGAUGAAUGCCUGUAUCCGGGAGUCAGAUG";
    std::string structure = ".......((((((.(((((............))))).)))))).....(((((((..((((((((.(((...((((((((((.......))))))...(((.(..((((((((((..............))))))))))....).)))....................(((...(((((((.((((((.............((((((((....))))).)))....((((..((((((((...))))))))...((((((.....((.............)).....))))))..(((.........))).........))))...)))))).)))))))...)))..(((....)))............))))..)))))))))))..(....)........))))))).....(((((...((....(((((((...(...(((((((.......))))))).))))))))....))((((....)))).)))))...";
    float energy = pipeline(sequence, structure);
    EXPECT_NEAR(energy, -121.2, 0.001);
}   