// RNAstructure (efn2) loop energies computed from the RNAstructure parameter sets.
//
// The expected values are from RNAstructure 6.6's `efn2 -w` (thermodynamic details), which
// prints kcal/mol with one decimal; Knotergy returns centicalories. Positions in the comments
// are 1-based like efn2's output, the arguments 0-based.

#include "test_utils.hpp"

#include <energy/modified_bases/ModHairpin.hpp>
#include <energy/modified_bases/ModInternal.hpp>
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

// The two pairs around an internal loop or bulge of a dot-bracket structure: the pair that
// opens at position i and the first pair inside it.
struct LoopPairs {
    std::size_t i;
    std::size_t j;
    std::size_t ci;
    std::size_t cj;
};

LoopPairs loop_at(const std::string& structure, std::size_t i) {
    const auto partner = [&](std::size_t open) {
        int depth = 0;
        for (std::size_t k = open; k < structure.size(); ++k) {
            if (structure[k] == '(') ++depth;
            if (structure[k] == ')' && --depth == 0) return k;
        }
        return structure.size();
    };
    const std::size_t ci = structure.find('(', i + 1);
    return {i, partner(i), ci, partner(ci)};
}

// Energy of the internal loop or bulge closed by the pair that opens at position i.
int internal(const std::string& sequence, const std::string& structure, std::size_t i,
             const knotergy::rnastructure_param& rsp) {
    const knotergy::ProcessedRNAEntry rna = process(sequence, structure);
    const LoopPairs loop = loop_at(structure, i);
    return RNAstructureFunctions::internal_loop_energy(
        loop.i, loop.j, loop.ci, loop.cj, rna.get_sequence(), rna.get_modified_sequence(), rsp);
}

// One loop of a structure and its energy from efn2 -w, in centicalories.
struct LoopCase {
    const char* what;
    std::string sequence;
    std::string structure;
    std::size_t i;  // 5' position of the closing pair (efn2 prints i + 1)
    int expected;
};

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

    // N on both sides of a helix of one pair: its helix end is forbidden (efn2 prints 1404.1)
    const knotergy::ProcessedRNAEntry unknown = process("ANCAAAAGNA", "..(....)..");
    is_inf = false;
    EXPECT_EQ(RNAstructureFunctions::hairpin_energy(2, 7, unknown.get_sequence(),
                                                    unknown.get_modified_sequence(), rsp, is_inf),
              INF);
    EXPECT_TRUE(is_inf);

    EXPECT_THROW((void) RNAstructureFunctions::hairpin_energy(
                     8, 3, small.get_sequence(), small.get_modified_sequence(), rsp, is_inf),
                 knotergy::DetailedException);
}

TEST(RNAstructureEnergies, InternalLoopsMatchEfn2) {
    const knotergy::rnastructure_param rsp = knotergy::RNAstructureParams::load();

    const std::string a16(16, 'A');
    const std::string a32(32, 'A');
    const std::vector<LoopCase> cases = {
        // 1x1, 2x1 and 2x2 loops have tables of their own
        {"1x1", "GGACAGCCAAAAGGCAGUCC", "((((.(((....))).))))", 3, 90},
        {"1x1 between A-U pairs", "GGAUAAGCAAAAGCUAAUCC", "((((.(((....))).))))", 3, 190},
        {"1x1 between G-U pairs", "GGGUAGGCAAAAGCUAGCCC", "((((.(((....))).))))", 3, 190},
        {"1x2", "GGACAGCCAAAAGGCAAGUCC", "((((.(((....)))..))))", 3, 230},
        {"2x1", "GGACAAGCCAAAAGGCAGUCC", "((((..(((....))).))))", 3, 230},
        {"2x2", "GGACAAGCCAAAAGGCAAGUCC", "((((..(((....)))..))))", 3, 120},
        {"2x2 with tandem G-A mismatches", "GGACGAGCCAAAAGGCGAGUCC", "((((..(((....)))..))))", 3,
         -70},
        // initiation + mismatches + asymmetry
        {"1x3", "GGACAGCCAAAAGGCAAAGUCC", "((((.(((....)))...))))", 3, 230},
        {"3x1", "GGACAAAGCCAAAAGGCAGUCC", "((((...(((....))).))))", 3, 230},
        {"2x3", "GGACAAGCCAAAAGGCAAAGUCC", "((((..(((....)))...))))", 3, 260},
        {"3x2", "GGACAAAGCCAAAAGGCAAGUCC", "((((...(((....)))..))))", 3, 260},
        {"3x3", "GGACAAAGCCAAAAGGCAAAGUCC", "((((...(((....)))...))))", 3, 200},
        {"2x5", "GGACAAGCCAAAAGGCAAAAAGUCC", "((((..(((....))).....))))", 3, 390},
        {"4x4 with a G-G mismatch", "GGACGAAAGCCAAAAGGCAAAGGUCC", "((((....(((....)))....))))", 3,
         130},
        // more than 30 nucleotides
        {"16x16", "GGAC" + a16 + "GCCAAAAGGC" + a16 + "GUCC",
         "((((" + std::string(16, '.') + "(((....)))" + std::string(16, '.') + "))))", 3, 380},
        {"1x32", "GGACAGCCAAAAGGC" + a32 + "GUCC",
         "((((.(((....)))" + std::string(32, '.') + "))))", 3, 680},
        {"32x1", "GGAC" + a32 + "GCCAAAAGGCAGUCC",
         "((((" + std::string(32, '.') + "(((....))).))))", 3, 680},
        {"closing pair at both sequence ends", "CAGCCAAAAGGCAG", "(.(((....))).)", 0, 90},
        // N is allowed inside a loop that is not looked up as a whole
        {"2x3 with N", "GGACANGCCAAAAGGCAAAGUCC", "((((..(((....)))...))))", 3, 260},
    };
    for (const LoopCase& c : cases) {
        EXPECT_EQ(internal(c.sequence, c.structure, c.i, rsp), c.expected)
            << c.what << ": " << c.sequence;
    }
}

TEST(RNAstructureEnergies, BulgesMatchEfn2) {
    const knotergy::rnastructure_param rsp = knotergy::RNAstructureParams::load();

    const std::vector<LoopCase> cases = {
        // One nucleotide: the pairs stack across it (-2.4 or -3.3 here) and the bulge costs 3.8,
        // minus RT ln(positions it can take), minus 0.9 for a C next to another C.
        {"one nucleotide", "GGACAGCCAAAAGGCGUCC", "((((.(((....)))))))", 3, 140},
        {"one nucleotide on the 3' side", "GGACGCCAAAAGGCAGUCC", "(((((((....))).))))", 3, 140},
        {"G before a G: 2 positions", "GGACGGCCAAAAGGCGUCC", "((((.(((....)))))))", 3, 100},
        {"C after a C: 2 positions, C bonus", "GGACCGCCAAAAGGCGUCC", "((((.(((....)))))))", 3, 10},
        {"C between Gs: no bonus", "GGAGCGCCAAAAGGCCUCC", "((((.(((....)))))))", 3, 50},
        {"C between Cs on the 3' side: 3 positions", "GGAGGCCAAAAGGCCCUCC", "(((((((....))).))))",
         3, -110},
        {"Gs up to the 5' end: 3 positions", "GGGCCAAAAGGCC", "(.(((....))))", 0, -20},
        {"Cs up to the 3' end: 3 positions", "GGCCAAAAGGCCC", "((((....))).)", 0, -110},
        {"N", "GGACNGCCAAAAGGCGUCC", "((((.(((....)))))))", 3, 140},
        // Longer bulges: initiation and the helix ends of both pairs
        {"two nucleotides", "GGACAAGCCAAAAGGCGUCC", "((((..(((....)))))))", 3, 280},
        {"three nucleotides", "GGACAAAGCCAAAAGGCGUCC", "((((...(((....)))))))", 3, 320},
        {"three nucleotides on the 3' side", "GGACGCCAAAAGGCAAAGUCC", "(((((((....)))...))))", 3,
         320},
        {"35 nucleotides", "GGAC" + std::string(35, 'A') + "GCCAAAAGGCGUCC",
         "((((" + std::string(35, '.') + "(((....)))))))", 3, 630},
    };
    for (const LoopCase& c : cases) {
        EXPECT_EQ(internal(c.sequence, c.structure, c.i, rsp), c.expected)
            << c.what << ": " << c.sequence;
    }
}

// efn2 -a m6A -w, efn2 -a PSU -w and efn2 -a 1 -w
TEST(RNAstructureEnergies, ModifiedBasesInInternalLoops) {
    const knotergy::rnastructure_param rsp =
        knotergy::RNAstructureParams::load("turner2004", {"m6A", "Ψ", "m1Ψ"});

    const std::vector<LoopCase> cases = {
        {"m6A in a 1x1 loop", "GGAC6GCCAAAAGGCAGUCC", "((((.(((....))).))))", 3, 90},
        {"m6A-U pair next to the closing pair", "GG6CAGCCAAAAGGCAGUCC", "((((.(((....))).))))", 3,
         90},
        {"Ψ in a 1x1 loop", "GGACPGCCAAAAGGCAGUCC", "((((.(((....))).))))", 3, 50},
        {"m1Ψ in a 1x1 loop", "GGAC1GCCAAAAGGCAGUCC", "((((.(((....))).))))", 3, 50},
        // A bulged A after an m6A-U pair has one position: m6A is not another A.
        {"A bulge after m6A", "GG6AGCCAAAAGGCUCC", "(((.(((....))))))", 2, 220},
        {"m6A bulge after m6A: 2 positions", "GG66GCCAAAAGGCUCC", "(((.(((....))))))", 2, 180},
    };
    for (const LoopCase& c : cases) {
        EXPECT_EQ(internal(c.sequence, c.structure, c.i, rsp), c.expected)
            << c.what << ": " << c.sequence;
    }
}

TEST(RNAstructureEnergies, ForbiddenInternalLoopsAreInfinite) {
    const knotergy::rnastructure_param rsp = knotergy::RNAstructureParams::load();

    // Pairs RNAstructure does not allow (efn2 prints 2800.0)
    EXPECT_EQ(internal("GGAAAGCCAAAAGGCAAUCC", "((((.(((....))).))))", 3, rsp), INF);
    EXPECT_EQ(internal("GGACAACCAAAAGGCAGUCC", "((((.(((....))).))))", 3, rsp), INF);

    // N in a 1x1 loop (1400.0), and N next to a helix of one pair: efn2 adds the forbidden
    // helix end to the -0.7 of the 2x2 loop and prints 1399.3.
    EXPECT_EQ(internal("GGACNGCCAAAAGGCAGUCC", "((((.(((....))).))))", 3, rsp), INF);
    EXPECT_EQ(internal("NCGAGCCAAAAGGCGAGN", ".(..(((....)))..).", 1, rsp), INF);

    const knotergy::ProcessedRNAEntry rna = process("GGACAGCCAAAAGGCAGUCC", "((((.(((....))).))))");
    const std::string& sequence = rna.get_sequence();
    const std::vector<std::string_view>& mod_sequence = rna.get_modified_sequence();

    // Without unpaired nucleotides the two pairs are a stack.
    EXPECT_EQ(
        RNAstructureFunctions::internal_loop_energy(0, 19, 1, 18, sequence, mod_sequence, rsp),
        stack(rna, 0, 19, rsp));

    // The inner pair has to be inside the outer one, and both inside the sequence.
    EXPECT_THROW((void) RNAstructureFunctions::internal_loop_energy(5, 14, 3, 16, sequence,
                                                                    mod_sequence, rsp),
                 knotergy::DetailedException);
    EXPECT_THROW((void) RNAstructureFunctions::internal_loop_energy(3, 20, 5, 14, sequence,
                                                                    mod_sequence, rsp),
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

TEST(AllParams, InternalFollowsTheLoadedParameters) {
    knotergy::vrna_md_param vp = knotergy::ViennaParams::load_energy_parameters(turner_file, 2);
    const knotergy::pk_param pkp =
        knotergy::PseudoknotParams::load_pk_param(pkp_file, knotergy::RoundMethod::None);
    std::vector<knotergy::modified_base_param> json_params =
        knotergy::ModParams::load_modified_energy_parameters(mod_folder);
    const knotergy::all_mod_params mp{json_params};
    const knotergy::rnastructure_param rsp =
        knotergy::RNAstructureParams::load("turner2004", {"m6A"});

    // 1x1 loop closed by pair 4-17, with an m6A-U pair before it (efn2 -a m6A -w: 0.9)
    const knotergy::ProcessedRNAEntry rna = process("GG6CAGCCAAAAGGCAGUCC", "((((.(((....))).))))");
    const std::string& sequence = rna.get_sequence();
    const std::vector<std::string_view>& mod_sequence = rna.get_modified_sequence();

    knotergy::AllParams params(mp, pkp, vp, rsp);

    // Nothing loaded
    EXPECT_EQ(params.internal(3, 16, 5, 14, sequence, mod_sequence), knotergy::NULL_ENERGY);

    // ViennaRNA only
    params.vienna_params_loaded = true;
    EXPECT_EQ(params.internal(3, 16, 5, 14, sequence, mod_sequence),
              knotergy::ViennaFunctions::internal_loop_energy(3, 16, 5, 14, sequence, vp));

    // ViennaRNA with the JSON modified-base corrections
    params.modified_params_loaded = true;
    EXPECT_EQ(params.internal(3, 16, 5, 14, sequence, mod_sequence),
              knotergy::ModInternal::find_mod_internal_energy(3, 16, 5, 14, mod_sequence, sequence,
                                                              vp, mp));

    // RNAstructure, once its parameters are loaded
    params.rnastructure_params_loaded = true;
    EXPECT_EQ(params.internal(3, 16, 5, 14, sequence, mod_sequence), 90);
}
