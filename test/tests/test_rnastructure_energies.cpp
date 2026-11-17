// RNAstructure (efn2) loop energies computed from the RNAstructure parameter sets.
//
// The expected values are from RNAstructure 6.6's `efn2 -w` (thermodynamic details), which
// prints kcal/mol with one decimal; Knotergy returns centicalories. Positions in the comments
// are 1-based like efn2's output, the arguments 0-based.

#include "test_utils.hpp"

#include <energy/modified_bases/ModHairpin.hpp>
#include <energy/rna_structure/RNAstructureFunctions.hpp>
#include <energy/rna_structure/RNAstructureUtils.hpp>
#include <energy/vienna/ViennaFunctions.hpp>
#include <io/parameters/AllParams.hpp>

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace {

using knotergy::RNAstructureFunctions;

// Knotergy's view of a sequence: the unmodified sequence plus the modified symbols.
knotergy::ProcessedRNAEntry process(const std::string& sequence, const std::string& structure) {
    std::vector<knotergy::modified_base_param> json_params =
        knotergy::ModParams::load_modified_energy_parameters(mod_folder);
    const knotergy::all_mod_params mp{json_params};
    return knotergy::RNAProcessor::process_rna(knotergy::RNAEntry(sequence, structure), mp);
}

int hairpin(const knotergy::ProcessedRNAEntry& rna, std::size_t i, std::size_t j,
            const knotergy::rnastructure_param& rsp) {
    bool is_inf = false;
    const int energy = RNAstructureFunctions::hairpin_energy(
        i, j, rna.get_sequence(), rna.get_modified_sequence(), rsp, is_inf);
    EXPECT_FALSE(is_inf);
    return energy;
}

int stack(const knotergy::ProcessedRNAEntry& rna, std::size_t i, std::size_t j,
          const knotergy::rnastructure_param& rsp) {
    return RNAstructureFunctions::stack_energy(i, j, i + 1, j - 1, rna.get_sequence(),
                                               rna.get_modified_sequence(), rsp);
}

}  // namespace

TEST(RNAstructureEnergies, HairpinsMatchEfn2) {
    const knotergy::rnastructure_param rsp = knotergy::RNAstructureParams::load();

    struct Case {
        const char* what;
        std::string sequence;
        std::string structure;
        std::size_t i;
        std::size_t j;
        int expected;
    };
    const std::vector<Case> cases = {
        {"initiation + first mismatch", "GGACAAAAGUCC", "((((....))))", 3, 8, 410},
        {"mismatch with a bonus", "GGACGAAAGUCC", "((((....))))", 3, 8, 330},
        {"five nucleotides", "GGACAGUAAGUCC", "((((.....))))", 3, 9, 420},
        {"six nucleotides", "GGACACAGUGGUCC", "((((......))))", 3, 10, 400},
        {"three nucleotides: no mismatch", "GGACGAAGUCC", "((((...))))", 3, 7, 540},
        {"special tetraloop CUUCGG", "GGACUUCGGUCC", "((((....))))", 3, 8, 370},
        {"special triloop CAACG", "GGACAACGUCC", "((((...))))", 3, 7, 680},
        {"special hexaloop ACAGUGCU", "GGCACAGUGCUGCC", "((((......))))", 3, 10, 290},
        {"all-C triloop", "GGAGCCCCUCC", "((((...))))", 3, 7, 690},
        {"all-C loop", "GGAGCCCCCCUCC", "((((.....))))", 3, 9, 810},
        {"GGG-U closure", "GGGAAAAUCC", "(((....)))", 2, 7, 360},
        {"GGG-U closure after an unpaired A", "AGGGACAAGUCCU", ".(((.....))).", 3, 9, 320},
        {"35 nucleotides", "GGAC" + std::string(35, 'A') + "GUCC",
         "((((" + std::string(35, '.') + "))))", 3, 39, 640},
        {"closing pair at both sequence ends", "CAAAAG", "(....)", 0, 5, 410},
    };
    for (const Case& c : cases) {
        const knotergy::ProcessedRNAEntry rna = process(c.sequence, c.structure);
        EXPECT_EQ(hairpin(rna, c.i, c.j, rsp), c.expected) << c.what << ": " << c.sequence;
    }
}

// efn2 -w on GGACAAAAGUCC ((((....)))): stacks -3.3, -2.4, -2.2, hairpin 4.1, exterior loop 0.0
// (only the helix end of pair 1-12), total -3.8 kcal/mol.
TEST(RNAstructureEnergies, LoopsAddUpToEfn2Total) {
    const knotergy::rnastructure_param rsp = knotergy::RNAstructureParams::load();
    const knotergy::ProcessedRNAEntry rna = process("GGACAAAAGUCC", "((((....))))");

    EXPECT_EQ(stack(rna, 0, 11, rsp), -330);
    EXPECT_EQ(stack(rna, 1, 10, rsp), -240);
    EXPECT_EQ(stack(rna, 2, 9, rsp), -220);
    const int exterior = RNAstructureFunctions::helix_end(0, 11, rna.get_sequence(),
                                                          rna.get_modified_sequence(), rsp);
    EXPECT_EQ(exterior, 0);
    EXPECT_EQ(stack(rna, 0, 11, rsp) + stack(rna, 1, 10, rsp) + stack(rna, 2, 9, rsp) +
                  hairpin(rna, 3, 8, rsp) + exterior,
              -380);

    // An A-U pair at the end of a helix costs 0.5 kcal/mol next to G-C.
    const knotergy::ProcessedRNAEntry au = process("AGAAAACU", "((....))");
    EXPECT_EQ(
        RNAstructureFunctions::helix_end(0, 7, au.get_sequence(), au.get_modified_sequence(), rsp),
        50);
}

// efn2 -a m6A -w: AGGAC6AAAGUCCU .((((....)))). has hairpin 4.1; GG6CAAAAGUCC ((((....))))
// has stacks -3.3, -1.8, -1.9 and hairpin 4.1. Loops with m6A come from the m6A set, the
// others from turner2004.
TEST(RNAstructureEnergies, ModifiedBasesUseTheirSet) {
    const knotergy::rnastructure_param rsp =
        knotergy::RNAstructureParams::load("turner2004", {"m6A", "Ψ", "m1Ψ"});

    const knotergy::ProcessedRNAEntry in_loop = process("AGGAC6AAAGUCCU", ".((((....)))).");
    ASSERT_TRUE(in_loop.has_modified_bases());
    EXPECT_EQ(hairpin(in_loop, 4, 9, rsp), 410);

    const knotergy::ProcessedRNAEntry in_stem = process("GG6CAAAAGUCC", "((((....))))");
    EXPECT_EQ(stack(in_stem, 0, 11, rsp), -330);
    EXPECT_EQ(stack(in_stem, 1, 10, rsp), -180);
    EXPECT_EQ(stack(in_stem, 2, 9, rsp), -190);
    EXPECT_EQ(hairpin(in_stem, 3, 8, rsp), 410);
}

TEST(RNAstructureEnergies, ForbiddenHairpinsAreInfinite) {
    const knotergy::rnastructure_param rsp = knotergy::RNAstructureParams::load();

    // Two loop nucleotides (efn2 prints 1400.0)
    const knotergy::ProcessedRNAEntry small = process("GGACAAGUCC", "((((..))))");
    bool is_inf = false;
    EXPECT_EQ(RNAstructureFunctions::hairpin_energy(3, 6, small.get_sequence(),
                                                    small.get_modified_sequence(), rsp, is_inf),
              INF);
    EXPECT_TRUE(is_inf);

    // A closing pair RNAstructure does not allow
    const knotergy::ProcessedRNAEntry unpaired = process("GGAAAAAAAUCC", "((((....))))");
    is_inf = false;
    EXPECT_EQ(RNAstructureFunctions::hairpin_energy(3, 8, unpaired.get_sequence(),
                                                    unpaired.get_modified_sequence(), rsp, is_inf),
              INF);
    EXPECT_TRUE(is_inf);

    EXPECT_THROW((void) RNAstructureFunctions::hairpin_energy(
                     8, 3, small.get_sequence(), small.get_modified_sequence(), rsp, is_inf),
                 knotergy::DetailedException);
}

// A symbol RNAstructure's sets do not know is an error that names the symbol, with or without
// modified bases in the sequence.
TEST(RNAstructureEnergies, UnknownSymbolsAreReported) {
    const knotergy::rnastructure_param rsp = knotergy::RNAstructureParams::load();
    const std::string sequence = "GGACAZAAGUCC";
    try {
        (void) knotergy::RNAstructureUtils::rnastructure_code(5, sequence, {}, rsp);
        FAIL();
    } catch (const knotergy::DetailedException& e) {
        EXPECT_NE(std::string(e.what()).find("symbol: Z at position 5"), std::string::npos)
            << e.what();
    }
}

TEST(AllParams, HairpinFollowsTheLoadedParameters) {
    knotergy::vrna_md_param vp = knotergy::ViennaParams::load_energy_parameters(turner_file, 2);
    const knotergy::pk_param pkp =
        knotergy::PseudoknotParams::load_pk_param(pkp_file, knotergy::RoundMethod::None);
    std::vector<knotergy::modified_base_param> json_params =
        knotergy::ModParams::load_modified_energy_parameters(mod_folder);
    const knotergy::all_mod_params mp{json_params};
    const knotergy::rnastructure_param rsp =
        knotergy::RNAstructureParams::load("turner2004", {"m6A"});

    const knotergy::ProcessedRNAEntry rna = process("GG6CAAAAGUCC", "((((....))))");
    const std::string& sequence = rna.get_sequence();
    const std::vector<std::string_view>& mod_sequence = rna.get_modified_sequence();
    bool is_inf = false;

    knotergy::AllParams params(mp, pkp, vp, rsp);

    // ViennaRNA only
    params.vienna_params_loaded = true;
    EXPECT_EQ(params.hairpin(3, 8, sequence, mod_sequence, is_inf),
              knotergy::ViennaFunctions::hairpin_energy(3, 8, sequence, is_inf, vp));

    // ViennaRNA with the JSON modified-base corrections
    params.modified_params_loaded = true;
    EXPECT_EQ(params.hairpin(3, 8, sequence, mod_sequence, is_inf),
              knotergy::ModHairpin::find_mod_hairpin_energy(3, 8, sequence, mod_sequence, vp, mp,
                                                            is_inf));

    // RNAstructure, once its parameters are loaded
    params.rnastructure_params_loaded = true;
    EXPECT_EQ(params.hairpin(3, 8, sequence, mod_sequence, is_inf), 410);
    EXPECT_FALSE(is_inf);
}
