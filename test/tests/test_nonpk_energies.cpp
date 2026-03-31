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

float pipeline(std::string sequence, std::string structure, std::string param_file, std::string pseudoknot_param_file = "../../params/pseudo/rna_pk_DirksPierce09_HotKnotsV2.json") {
    knotergy::vrna_md_param vp = knotergy::ViennaParams::load_energy_parameters(param_file);
    knotergy::pk_param pkp = knotergy::PseudoknotParams::load_pk_param(pseudoknot_param_file);
    knotergy::RNAEntry rna(sequence, structure);
    knotergy::ProcessedRNAEntry processed_rna(knotergy::RNAProcessor::process_rna(std::move(rna)));
    knotergy::LoopFactory factory(processed_rna);
    std::vector<knotergy::modified_base_param> mp; // empty for no modified bases
    knotergy::ComputeEnergy energy(factory.get_root_node(), processed_rna, vp, pkp, mp);

    return energy.getEnergy();
}
}

const std::string turner = "../../params/common/rna_turner2004.par";
const std::string DP = "../../params/common/rna_DirksPierce09.par";

// echo -e "AU\n()" | RNAeval
TEST(NonPseudoKnottedEnergies, nonPK_small_Turner) {
    std::string sequence = "AU";
    std::string structure = "()";
    float turner_energy = pipeline(sequence, structure, turner);

    EXPECT_NEAR(turner_energy, 100000.50, 0.001); // Turner 2004
}

// echo -e "AU\n()" | RNAeval
TEST(NonPseudoKnottedEnergies, nonPK_small_DP) {
    std::string sequence = "AU";
    std::string structure = "()";
    float dp_energy = pipeline(sequence, structure, DP);

    EXPECT_NEAR(dp_energy, 100000.50, 0.001); // Dirks & Pierce 2009
}

// echo -e "AGGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC\n.(((((.........................)))))................" | RNAeval
TEST(NonPseudoKnottedEnergies, BasicStack_Turner) {
    std::string sequence = "AGGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    std::string structure = ".(((((.........................)))))................";

    float turner_energy = pipeline(sequence, structure, turner);
    EXPECT_NEAR(turner_energy, -8.30, 0.001); // Turner 2004
}

// echo -e "AGGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC\n.(((((.........................)))))................" | RNAeval
TEST(NonPseudoKnottedEnergies, BasicStack_DP) {
    std::string sequence = "AGGGGGAAAAAAAGGGGGGGGGGAAAAAAAACCCCCAAAAAACCCCCCCCCC";
    std::string structure = ".(((((.........................)))))................";

    float dp_energy = pipeline(sequence, structure, DP);
    EXPECT_NEAR(dp_energy, -5.91, 0.001); // Dirks & Pierce 2009
}

// echo -e "GGGGAGAAAAAAAAAUUUUUU\n((((((.........))))))" | RNAeval
TEST(NonPseudoKnottedEnergies, StacksFullStructure_Turner) {
    std::string sequence = "GGGGAGAAAAAAAAAUUUUUU";
    std::string structure = "((((((.........))))))";

    float turner_energy = pipeline(sequence, structure, turner);
    EXPECT_NEAR(turner_energy, 3.70, 0.001);
}

// echo -e "GGGGAGAAAAAAAAAUUUUUU\n((((((.........))))))" | RNAeval
TEST(NonPseudoKnottedEnergies, StacksFullStructure_DP) {
    std::string sequence = "GGGGAGAAAAAAAAAUUUUUU";
    std::string structure = "((((((.........))))))";

    float dp_energy = pipeline(sequence, structure, DP);
    EXPECT_NEAR(dp_energy, 0.87, 0.001); // Dirks & Pierce 2009
}

// echo -e "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCC\n((((((.....................))))))................((((((......))))))...(((...)))......" | RNAeval
TEST(NonPseudoKnottedEnergies, MultipleStacks_Turner) {
    std::string sequence = "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCC";
    std::string structure = "((((((.....................))))))................((((((......))))))...(((...)))......";

    float turner_energy = pipeline(sequence, structure, turner);
    EXPECT_NEAR(turner_energy, -28.10, 0.001);
}

// echo -e "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCC\n((((((.....................))))))................((((((......))))))...(((...)))......" | RNAeval
TEST(NonPseudoKnottedEnergies, MultipleStacks_DP) {
    std::string sequence = "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCC";
    std::string structure = "((((((.....................))))))................((((((......))))))...(((...)))......";

    float dp_energy = pipeline(sequence, structure, DP);
    EXPECT_NEAR(dp_energy, -19.21, 0.001); // Dirks & Pierce 2009
}

// echo -e "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCCA\n...........((((((....((((((............))))))....((((((......))))))...((....)).))))))." | RNAeval
TEST(NonPseudoKnottedEnergies, MultiLoopWithStacks_Turner) {
    std::string sequence = "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCCA";
    std::string structure = "...........((((((....((((((............))))))....((((((......))))))...((....)).)))))).";

    float turner_energy = pipeline(sequence, structure, turner);
    EXPECT_NEAR(turner_energy, -41.10, 0.001);
}

// echo -e "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCCA\n...........((((((....((((((............))))))....((((((......))))))...((....)).))))))." | RNAeval
TEST(NonPseudoKnottedEnergies, MultiLoopWithStacks_DP) {
    std::string sequence = "GGGGGGAAAAAGGGGGGAAAAGGGGGGCCCCCCAAAAAACCCCCCAAAAGGGGGGAAAAAACCCCCCAAAGGGAAACCCCCCCCCA";
    std::string structure = "...........((((((....((((((............))))))....((((((......))))))...((....)).)))))).";

    float dp_energy = pipeline(sequence, structure, DP);
    EXPECT_NEAR(dp_energy, -29.38, 0.001); // Dirks & Pierce 2009
}

// echo -e "AAAAAAGGGGCCCCCCAAAAAAAAAAAAGGGGGGAAAAGGGGGGUUUUUUCCCCCCAAAUUUUUUAAGUUUUUU\n((((((....((((((............))))))....((((((......))))))...((....)).))))))" | RNAeval
TEST(NonPseudoKnottedEnergies, MultiLoopWithStacksFull_Turner) {
    std::string sequence = "AAAAAAGGGGCCCCCCAAAAAAAAAAAAGGGGGGAAAAGGGGGGUUUUUUCCCCCCAAAUUUUUUAAGUUUUUU";
    std::string structure = "((((((....((((((............))))))....((((((......))))))...((....)).))))))";

    float turner_energy = pipeline(sequence, structure, turner);
    EXPECT_NEAR(turner_energy, -21.9, 0.001);
}

// echo -e "AAAAAAGGGGCCCCCCAAAAAAAAAAAAGGGGGGAAAAGGGGGGUUUUUUCCCCCCAAAUUUUUUAAGUUUUUU\n((((((....((((((............))))))....((((((......))))))...((....)).))))))" | RNAeval
TEST(NonPseudoKnottedEnergies, MultiLoopWithStacksFull_DP) {
    std::string sequence = "AAAAAAGGGGCCCCCCAAAAAAAAAAAAGGGGGGAAAAGGGGGGUUUUUUCCCCCCAAAUUUUUUAAGUUUUUU";
    std::string structure = "((((((....((((((............))))))....((((((......))))))...((....)).))))))";

    float dp_energy = pipeline(sequence, structure, DP);
    EXPECT_NEAR(dp_energy, -16.47, 0.001); // Dirks & Pierce 2009
}

// echo -e "GGGGUUAUUUUAUUAAAAAUAACCCUGGUUUUUAAGGCGGGGUCGUGCGGUAAGGGAACCC\n((((..(...).((...)))..(((.(.((...))..(...).)..(...)..)))..)))" | RNAeval
TEST(NonPseudoKnottedEnergies, MultiWithMultiLoop_Turner) {
    std::string sequence = "GGGGUUAUUUUAUUAAAAAUAACCCUGGUUUUUAAGGCGGGGUCGUGCGGUAAGGGAACCC";
    std::string structure = "((((..(...).((...)))..(((.(.((...))..(...).)..(...)..)))..)))";

    float turner_energy = pipeline(sequence, structure, turner);
    EXPECT_NEAR(turner_energy, 30.90, 0.001);
}

// echo -e "GGGGUUAUUUUAUUAAAAAUAACCCUGGUUUUUAAGGCGGGGUCGUGCGGUAAGGGAACCC\n((((..(...).((...)))..(((.(.((...))..(...).)..(...)..)))..)))" | RNAeval
TEST(NonPseudoKnottedEnergies, MultiWithMultiLoop_DP) {
    std::string sequence = "GGGGUUAUUUUAUUAAAAAUAACCCUGGUUUUUAAGGCGGGGUCGUGCGGUAAGGGAACCC";
    std::string structure = "((((..(...).((...)))..(((.(.((...))..(...).)..(...)..)))..)))";

    float dp_energy = pipeline(sequence, structure, DP);
    EXPECT_NEAR(dp_energy, 19.81, 0.001); // Dirks & Pierce 2009
}

TEST(NonPseudoKnottedEnergies, PositiveInternalLoops_Turner) {
    std::string sequence  = "GGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAAUUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGCGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACG";
    std::string structure = "(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))..........((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..........))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((....)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))";

    float turner_energy = pipeline(sequence, structure, turner);
    EXPECT_NEAR(turner_energy, -471.90, 0.001);
}

TEST(NonPseudoKnottedEnergies, PositiveInternalLoops_DP) {
    std::string sequence  = "GGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAACCCUUUGGGAAAUUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGGAUUCGAGCGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACGUACG";
    std::string structure = "(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...(((...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))...)))..........((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..((((((..........))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))..))))))(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((....)))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))";

    float dp_energy = pipeline(sequence, structure, DP);
    EXPECT_NEAR(dp_energy, -317.37, 0.001); // Dirks & Pierce 2009
}

TEST(NonPseudoKnottedEnergies, MoreLongSequences_Turner) {
    std::string sequence  = "UUAAAAGGGAUGCCUCUCCUGUUCAUCUUGUGGAGAAGCAUUCGAUAAGGUCAUCAUAAUGGGUCCAGCUUUGCGACCUGGCGAGAUUAGUCAGGAAAAUGUGAAGUGGGUCUUCGCUUUCCAGGUACAGGAGGCUCGCCCCGCUCAUCCAGUUCGUCCCCUAACCACUUGUUUUCUCAGGAUAGUUUGUUUUGUACACCCGUGUACAUACACAUGUAUCACACCCCAGAUUGCCGAAUGUUUCGUUCGGUCGAGCCUGACUAUGCAUAAACCCUACCUCUGAAACCUUGGGCAACUCACUACUUCCGAGCUAAAUCCCUCUGUUUGAGCUAGCCUGAGAUUUCAACUGGCUUCGGCCUUGUUUAUACCAUCGUUUGCUGAUCCAUUGAAGAAAUAAGUUACCGAUGGCCCCAAACUGACGAUCACUAUUCUUUCCAUAGGAGUUAUGGGUAUACUGCCCGUAGACGGAAAGAUGAAUGCCUGUAUCCGGGAGUCAGAUG";
    std::string structure = ".......((((((.(((((............))))).)))))).....(((((((..((((((((.(((...((((((((((.......))))))...(((.(..((((((((((..............))))))))))....).)))....................(((...(((((((.((((((.............((((((((....))))).)))....((((..((((((((...))))))))...((((((.....((.............)).....))))))..(((.........))).........))))...)))))).)))))))...)))..(((....)))............))))..)))))))))))..(....)........))))))).....(((((...((....(((((((...(...(((((((.......))))))).))))))))....))((((....)))).)))))...";

    float turner_energy = pipeline(sequence, structure, turner);
    EXPECT_NEAR(turner_energy, -121.2, 0.001);
}

TEST(NonPseudoKnottedEnergies, MoreLongSequences_DP) {
    std::string sequence  = "UUAAAAGGGAUGCCUCUCCUGUUCAUCUUGUGGAGAAGCAUUCGAUAAGGUCAUCAUAAUGGGUCCAGCUUUGCGACCUGGCGAGAUUAGUCAGGAAAAUGUGAAGUGGGUCUUCGCUUUCCAGGUACAGGAGGCUCGCCCCGCUCAUCCAGUUCGUCCCCUAACCACUUGUUUUCUCAGGAUAGUUUGUUUUGUACACCCGUGUACAUACACAUGUAUCACACCCCAGAUUGCCGAAUGUUUCGUUCGGUCGAGCCUGACUAUGCAUAAACCCUACCUCUGAAACCUUGGGCAACUCACUACUUCCGAGCUAAAUCCCUCUGUUUGAGCUAGCCUGAGAUUUCAACUGGCUUCGGCCUUGUUUAUACCAUCGUUUGCUGAUCCAUUGAAGAAAUAAGUUACCGAUGGCCCCAAACUGACGAUCACUAUUCUUUCCAUAGGAGUUAUGGGUAUACUGCCCGUAGACGGAAAGAUGAAUGCCUGUAUCCGGGAGUCAGAUG";
    std::string structure = ".......((((((.(((((............))))).)))))).....(((((((..((((((((.(((...((((((((((.......))))))...(((.(..((((((((((..............))))))))))....).)))....................(((...(((((((.((((((.............((((((((....))))).)))....((((..((((((((...))))))))...((((((.....((.............)).....))))))..(((.........))).........))))...)))))).)))))))...)))..(((....)))............))))..)))))))))))..(....)........))))))).....(((((...((....(((((((...(...(((((((.......))))))).))))))))....))((((....)))).)))))...";

    float dp_energy = pipeline(sequence, structure, DP);
    EXPECT_NEAR(dp_energy, -90.42, 0.001); // Dirks & Pierce 2009
}