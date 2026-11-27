// RNAstructure (efn2) loop energies computed from the RNAstructure parameter sets.
//
// The expected values are from RNAstructure 6.6's `efn2 -w` (thermodynamic details), which
// prints kcal/mol with one decimal; Knotergy returns centicalories. Positions in the comments
// are 1-based like efn2's output, the arguments 0-based.

#include "test_utils.hpp"

#include <RNAstructure/eval/exterior.hpp>
#include <RNAstructure/eval/hairpin.hpp>
#include <RNAstructure/eval/internal.hpp>
#include <RNAstructure/eval/loop_element.hpp>
#include <RNAstructure/eval/multibranch.hpp>
#include <energy/EnergyFunctions.hpp>
#include <energy/modified_bases/ModBaseUtils.hpp>
#include <energy/modified_bases/ModExternal.hpp>
#include <energy/modified_bases/ModHairpin.hpp>
#include <energy/modified_bases/ModInternal.hpp>
#include <energy/modified_bases/ModMultiloop.hpp>
#include <energy/rna_structure/RNAstructureFunctions.hpp>
#include <energy/rna_structure/RNAstructureUtils.hpp>
#include <energy/vienna/ViennaFunctions.hpp>
#include <io/parameters/AllParams.hpp>

#include <algorithm>
#include <cstddef>
#include <string>
#include <string_view>
#include <utility>
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

// Position paired with the one that opens at position i of a dot-bracket structure.
std::size_t partner(const std::string& structure, std::size_t i) {
    int depth = 0;
    for (std::size_t k = i; k < structure.size(); ++k) {
        if (structure[k] == '(') ++depth;
        if (structure[k] == ')' && --depth == 0) return k;
    }
    return structure.size();
}

LoopPairs loop_at(const std::string& structure, std::size_t i) {
    const std::size_t ci = structure.find('(', i + 1);
    return {i, partner(structure, i), ci, partner(structure, ci)};
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

// A dangling end on the pair that opens at position i, and its energy in centicalories.
struct DangleCase {
    const char* what;
    std::string sequence;
    std::string structure;
    std::size_t i;    // 5' position of the pair
    char side;        // '5' or '3'
    bool is_closing;  // the dangle is inside the pair
    int expected;
};

int dangle(const DangleCase& c, const knotergy::rnastructure_param& rsp) {
    const knotergy::ProcessedRNAEntry rna = process(c.sequence, c.structure);
    const std::size_t j = partner(c.structure, c.i);
    return c.side == '5'
               ? RNAstructureFunctions::dangle5_energy(
                     c.i, j, rna.get_sequence(), rna.get_modified_sequence(), rsp, c.is_closing)
               : RNAstructureFunctions::dangle3_energy(
                     c.i, j, rna.get_sequence(), rna.get_modified_sequence(), rsp, c.is_closing);
}

// A structure and one of its energies from efn2 -w, in centicalories.
struct StructureCase {
    const char* what;
    std::string sequence;
    std::string structure;
    int expected;
};

namespace rs = knotergy::rnastructure;

// Energy of the exterior loop of a structure.
int exterior(const std::string& sequence, const std::string& structure,
             const knotergy::rnastructure_param& rsp,
             rs::CoaxialStacking coaxial = rs::CoaxialStacking::On) {
    const knotergy::ProcessedRNAEntry rna = process(sequence, structure);
    knotergy::LoopFactory factory(rna);
    return RNAstructureFunctions::external_energy(factory.get_root_node().children, rna, rsp,
                                                  coaxial);
}

// Energy of the multibranch loop closed by the pair that opens at position i.
int multibranch(const std::string& sequence, const std::string& structure, std::size_t i,
                const knotergy::rnastructure_param& rsp,
                rs::MultibranchRule rule = rs::MultibranchRule::Efn2,
                rs::CoaxialStacking coaxial = rs::CoaxialStacking::On) {
    const knotergy::ProcessedRNAEntry rna = process(sequence, structure);
    knotergy::LoopFactory factory(rna);
    return RNAstructureFunctions::multibranch_energy(*factory.get_node_table()[i], rna, rsp, rule,
                                                     coaxial);
}

// Energy of a structure: its exterior loop and every loop inside it.
int total(const std::string& sequence, const std::string& structure,
          const knotergy::rnastructure_param& rsp,
          rs::MultibranchRule rule = rs::MultibranchRule::Efn2,
          rs::CoaxialStacking coaxial = rs::CoaxialStacking::On) {
    const knotergy::ProcessedRNAEntry rna = process(sequence, structure);
    knotergy::LoopFactory factory(rna);
    int energy =
        RNAstructureFunctions::external_energy(factory.get_root_node().children, rna, rsp, coaxial);
    for (std::size_t i = 0; i < structure.size(); ++i) {
        if (structure[i] != '(') continue;
        const knotergy::LoopNode& loop = *factory.get_node_table()[i];
        if (loop.loop_type == knotergy::LoopType::Hairpin) {
            energy += hairpin(rna, loop.begin, loop.end, rsp);
        } else if (loop.loop_type == knotergy::LoopType::Multibranch) {
            energy += RNAstructureFunctions::multibranch_energy(loop, rna, rsp, rule, coaxial);
        } else {  // internal loop or bulge, or a stack when nothing is unpaired
            energy += RNAstructureFunctions::internal_loop_energy(
                loop.begin, loop.end, loop.children[0]->begin, loop.children[0]->end,
                rna.get_sequence(), rna.get_modified_sequence(), rsp);
        }
    }
    return energy;
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

// efn2 -w does not print dangles, so the values are differences of its loop energies: the loop
// with the nucleotide minus the same loop with N in its place, which never stacks.
// - Outside a pair: one nucleotide next to the only helix of the exterior loop.
// - Inside a pair: one nucleotide next to the pair that closes a multibranch loop whose other
//   nucleotides are all N.
TEST(RNAstructureEnergies, DanglesMatchEfn2) {
    const knotergy::rnastructure_param rsp = knotergy::RNAstructureParams::load();

    const std::string before = ".((((....))))";
    const std::string after = "((((....)))).";
    const std::string multiloop = "((..((((....))))..((((....))))..))";
    const std::vector<DangleCase> cases = {
        // "Exterior loop" of efn2 -w, and 0.0 with N
        {"3' A after a C-G helix", "CGACAAAAGUCGA", after, 0, '3', false, -110},
        {"3' C after a C-G helix", "CGACAAAAGUCGC", after, 0, '3', false, -40},
        {"3' G after a G-C helix", "GGACAAAAGUCCG", after, 0, '3', false, -170},
        {"5' A before a C-G helix", "ACGACAAAAGUCG", before, 1, '5', false, -50},
        {"5' A before a G-C helix", "AGGACAAAAGUCC", before, 1, '5', false, -20},
        // -0.2, 0.3 and 0.4, where N gives the 0.5 of the helix end
        {"3' A after an A-U helix", "AGACAAAAGUCUA", after, 0, '3', false, -70},
        {"5' G before a U-A helix", "GUGACAAAAGUCA", before, 1, '5', false, -20},
        {"3' U after a G-U helix", "GGACAAAAGUCUU", after, 0, '3', false, -10},
        // "Multiloop" of efn2 -w: 5.8, 7.3 and 6.2 against 7.5 with N, 7.7 against 8.0
        {"3' A inside a closing C-G", "GCANGGACAAAAGUCCNNGGACAAAAGUCCNNGC", multiloop, 1, '3', true,
         -170},
        {"5' A inside a closing C-G", "GCNNGGACAAAAGUCCNNGGACAAAAGUCCNAGC", multiloop, 1, '5', true,
         -20},
        {"3' G inside a closing G-C", "GGGNGGACAAAAGUCCNNGGACAAAAGUCCNNCC", multiloop, 1, '3', true,
         -130},
        {"5' A inside a closing A-U", "GANNGGACAAAAGUCCNNGGACAAAAGUCCNAUC", multiloop, 1, '5', true,
         -30},
        // N stacks with 0 on the 3' side and is not allowed on the 5' side
        {"3' N", "CGACAAAAGUCGN", after, 0, '3', false, 0},
        {"5' N", "NCGACAAAAGUCG", before, 1, '5', false, INF},
        // A pair RNAstructure does not allow
        {"3' A after an A-A pair", "AGACAAAAGUCAA", after, 0, '3', false, INF},
    };
    for (const DangleCase& c : cases) {
        EXPECT_EQ(dangle(c, rsp), c.expected) << c.what << ": " << c.sequence;
    }

    // efn2 -w on AGACAAAAGUCUA ((((....)))). has Exterior loop = -0.2: the end of the helix and
    // the A after it.
    const knotergy::ProcessedRNAEntry rna = process("AGACAAAAGUCUA", after);
    const std::string& sequence = rna.get_sequence();
    const std::vector<std::string_view>& mod_sequence = rna.get_modified_sequence();
    EXPECT_EQ(RNAstructureFunctions::helix_end(0, 11, sequence, mod_sequence, rsp) +
                  RNAstructureFunctions::dangle3_energy(0, 11, sequence, mod_sequence, rsp),
              -20);

    // There has to be a nucleotide on that side of the pair.
    EXPECT_THROW((void) RNAstructureFunctions::dangle5_energy(0, 11, sequence, mod_sequence, rsp),
                 knotergy::DetailedException);
    EXPECT_THROW((void) RNAstructureFunctions::dangle3_energy(1, 12, sequence, mod_sequence, rsp),
                 knotergy::DetailedException);
    EXPECT_THROW(
        (void) RNAstructureFunctions::dangle3_energy(3, 4, sequence, mod_sequence, rsp, true),
        knotergy::DetailedException);
    EXPECT_THROW((void) RNAstructureFunctions::dangle5_energy(11, 0, sequence, mod_sequence, rsp),
                 knotergy::DetailedException);
}

// efn2 -a m6A -w, efn2 -a PSU -w and efn2 -a 1 -w, as in DanglesMatchEfn2
TEST(RNAstructureEnergies, DanglesWithModifiedBases) {
    const knotergy::rnastructure_param rsp =
        knotergy::RNAstructureParams::load("turner2004", {"m6A", "Ψ", "m1Ψ"});

    const std::string before = ".((((....))))";
    const std::string after = "((((....)))).";
    const std::string multiloop = "((..((((....))))..((((....))))..))";
    const std::vector<DangleCase> cases = {
        {"5' m6A before an A-U helix", "6AGACAAAAGUCU", before, 1, '5', false, -90},
        {"3' m6A after an A-U helix", "AGACAAAAGUCU6", after, 0, '3', false, -100},
        {"3' m6A after a G-C helix", "GGACAAAAGUCC6", after, 0, '3', false, -130},
        {"3' A after an m6A-U helix", "6GACAAAAGUCUA", after, 0, '3', false, -90},
        {"5' A inside a closing U-m6A", "GUNNGGACAAAAGUCCNNGGACAAAAGUCCNA6C", multiloop, 1, '5',
         true, -40},
        {"3' Ψ after a G-C helix", "GGACAAAAGUCCP", after, 0, '3', false, -150},
        {"5' Ψ before a G-C helix", "PGGACAAAAGUCC", before, 1, '5', false, -50},
        {"3' m1Ψ after a G-C helix", "GGACAAAAGUCC1", after, 0, '3', false, -120},
        {"3' A after an A-m1Ψ helix", "AGACAAAAGUC1A", after, 0, '3', false, -160},
        {"3' A inside a closing m1Ψ-G", "G1ANGGACAAAAGUCCNNGGACAAAAGUCCNNGC", multiloop, 1, '3',
         true, -160},
    };
    for (const DangleCase& c : cases) {
        EXPECT_EQ(dangle(c, rsp), c.expected) << c.what << ": " << c.sequence;
    }
}

// As for the dangles, efn2 -w shows a terminal mismatch only through its loop: the loop with
// the two nucleotides minus the same loop with N for both. In these cases the mismatch is lower
// than either dangle, so the difference is the mismatch.
TEST(RNAstructureEnergies, MismatchesMatchEfn2) {
    const knotergy::rnastructure_param rsp = knotergy::RNAstructureParams::load();

    struct MismatchCase {
        const char* what;
        std::string sequence;
        std::string structure;
        std::size_t i;  // 5' position of the pair
        bool is_external;
        bool is_closing;
        int expected;
    };
    const std::string helix = ".((((....)))).";
    const std::string multiloop = "((..((((....))))..((((....))))..))";
    const std::vector<MismatchCase> cases = {
        // Exterior loop = -1.3, and 0.0 with N
        {"G and A around a C-G helix", "GCGACAAAAGUCGA", helix, 1, true, false, -130},
        // Multiloop = 6.9 and 7.3, and 8.0 with N
        {"A and G around an A-U branch", "GCNAAGACAAAAGUCUGNGGACAAAAGUCCNNGC", multiloop, 4, false,
         false, -110},
        {"C and C around a U-A branch", "GCNCUGACAAAAGUCACNGGACAAAAGUCCNNGC", multiloop, 4, false,
         false, -70},
        // Multiloop = 6.2 against 7.5 with N, and 7.3 against 8.0
        {"A and G inside a closing G-C", "GGANGGACAAAAGUCCNNGGACAAAAGUCCNGCC", multiloop, 1, false,
         true, -130},
        {"C and A inside a closing U-G", "GUCNGGACAAAAGUCCNNGGACAAAAGUCCNAGC", multiloop, 1, false,
         true, -70},
        // A pair RNAstructure does not allow
        {"G and A around a G-A pair", "GGGACAAAAGUCAA", helix, 1, true, false, INF},
    };
    for (const MismatchCase& c : cases) {
        const knotergy::ProcessedRNAEntry rna = process(c.sequence, c.structure);
        EXPECT_EQ(RNAstructureFunctions::mismatch_energy(
                      c.i, partner(c.structure, c.i), rna.get_sequence(),
                      rna.get_modified_sequence(), rsp, c.is_external, c.is_closing),
                  c.expected)
            << c.what << ": " << c.sequence;
    }

    // N does not stack: a mismatch with N on either side is 0.
    const knotergy::ProcessedRNAEntry n5 = process("NCGACAAAAGUCGA", helix);
    EXPECT_EQ(RNAstructureFunctions::mismatch_energy(1, 12, n5.get_sequence(),
                                                     n5.get_modified_sequence(), rsp, true),
              0);
    const knotergy::ProcessedRNAEntry n3 = process("GCGACAAAAGUCGN", helix);
    EXPECT_EQ(RNAstructureFunctions::mismatch_energy(1, 12, n3.get_sequence(),
                                                     n3.get_modified_sequence(), rsp, true),
              0);

    // There has to be a nucleotide on both sides, and an exterior loop has no closing pair.
    const knotergy::ProcessedRNAEntry rna = process("GCGACAAAAGUCGA", helix);
    const std::string& sequence = rna.get_sequence();
    const std::vector<std::string_view>& mod_sequence = rna.get_modified_sequence();
    EXPECT_THROW(
        (void) RNAstructureFunctions::mismatch_energy(0, 12, sequence, mod_sequence, rsp, true),
        knotergy::DetailedException);
    EXPECT_THROW(
        (void) RNAstructureFunctions::mismatch_energy(1, 13, sequence, mod_sequence, rsp, true),
        knotergy::DetailedException);
    EXPECT_THROW((void) RNAstructureFunctions::mismatch_energy(1, 12, sequence, mod_sequence, rsp,
                                                               true, true),
                 knotergy::DetailedException);
}

// "Exterior loop" of efn2 -w
TEST(RNAstructureEnergies, ExteriorLoopsMatchEfn2) {
    const knotergy::rnastructure_param rsp = knotergy::RNAstructureParams::load();

    const std::string h = "((((....))))";
    const std::vector<StructureCase> cases = {
        {"no helix", "AAAAAA", "......", 0},
        // The end of the helix: 0.5 for A-U and G-U
        {"G-C helix", "GGACAAAAGUCC", h, 0},
        {"A-U helix", "AGACAAAAGUCU", h, 50},
        {"G-U helix", "GGACAAAAGUCU", h, 50},
        // Dangles and terminal mismatches
        {"3' dangle", "GGACAAAAGUCCG", h + ".", -170},
        {"5' dangle", "AGGACAAAAGUCC", "." + h, -20},
        {"mismatch no better than the 3' dangle", "AGGACAAAAGUCCG", "." + h + ".", -170},
        {"terminal mismatch", "GCGACAAAAGUCGA", "." + h + ".", -130},
        {"terminal mismatch among more nucleotides", "AAGCGACAAAAGUCGAAA", "..." + h + "...", -130},
        {"helix of one pair", "GCAAAAGCA", ".(....)..", -110},
        // Two helices without a nucleotide between them stack coaxially
        {"flush stack", "GGACAAAAGUCCGGACAAAAGUCC", h + h, -240},
        {"flush stack of A-U ends", "AGACAAAAGUCUAGACAAAAGUCU", h + h, -30},
        {"flush stack, nucleotides around", "AGGACAAAAGUCCGGACAAAAGUCCA", "." + h + h + ".", -240},
        {"three helices: one flush stack", "GGACAAAAGUCCGGACAAAAGUCCGGACAAAAGUCC", h + h + h, -240},
        // ... or across one nucleotide, which needs a nucleotide on the other side of a helix
        {"stack across a nucleotide, mismatch after", "GGACAAAAGUCCAGGACAAAAGUCCA",
         h + "." + h + ".", -370},
        {"stack across a nucleotide, mismatch before", "AGGACAAAAGUCCAGGACAAAAGUCC",
         "." + h + "." + h, -370},
        {"one nucleotide between and none outside: 3' dangle", "GGACAAAAGUCCAGGACAAAAGUCC",
         h + "." + h, -170},
        {"two nucleotides between: both dangles", "GGACAAAAGUCCAAGGACAAAAGUCC", h + ".." + h, -190},
        {"three nucleotides between", "GGACAAAAGUCCAGAGGACAAAAGUCC", h + "..." + h, -190},
        {"three helices: stack across a nucleotide, then a dangle",
         "AGGACAAAAGUCCAGGACAAAAGUCCGGACAAAAGUCCA", "." + h + "." + h + h + ".", -540},
        {"four different helices", "CAGACAAAAGUCUGGCGGAAACGCCGACAAAAGUCGAAGGACAAAAGUCUU",
         ".((((....)))).(((....)))((((....))))..((((....)))).", -410},
        // N does not stack
        {"N on both sides", "NGGACAAAAGUCCN", "." + h + ".", 0},
        // A helix end RNAstructure does not allow: N next to the pair of a helix of one pair
        // (efn2 prints 1398.9), a G-A pair (1400.0)
        {"N inside a helix of one pair", "GCNAAAGCA", ".(....)..", INF},
        {"pair that cannot form", "GGACAAAAGUCAA", h + ".", INF},
    };
    for (const StructureCase& c : cases) {
        EXPECT_EQ(exterior(c.sequence, c.structure, rsp), c.expected)
            << c.what << ": " << c.sequence;
    }

    // Nothing stacks on a pseudoknot, so the helices on its two sides are scored apart.
    const std::string hairpin = "GGACAAAAGUCC";
    const std::string knot = "GGCAAGCGAAGCCAACGC";
    const std::string knot_structure = "(((..[[[..)))..]]]";
    EXPECT_EQ(exterior("A" + hairpin + "A" + knot + "A" + hairpin + "G",
                       "." + h + "." + knot_structure + "." + h + ".", rsp),
              exterior("A" + hairpin + "A", "." + h + ".", rsp) +
                  exterior("A" + hairpin + "G", "." + h + ".", rsp));
    EXPECT_EQ(exterior(hairpin + knot + hairpin, h + knot_structure + h, rsp), 0);
}

// "Multiloop" of efn2 -w. The loop is closed by the second pair of the structure (2-N in
// efn2's output), except where noted.
TEST(RNAstructureEnergies, MultibranchLoopsMatchEfn2) {
    const knotergy::rnastructure_param rsp = knotergy::RNAstructureParams::load();

    const std::string h = "((((....))))";
    const std::vector<LoopCase> cases = {
        // 9.3 - 0.6 per helix, plus a strain of 3.1 when an odd number of helices has at most
        // one unpaired nucleotide among them. The rest is the stacking of the helices.
        {"three helices, nothing unpaired", "GCGGACAAAAGUCCGGACAAAAGUCCGC", "((" + h + h + "))", 1,
         820},
        {"three helices, one nucleotide", "GCGGACAAAAGUCCAGGACAAAAGUCCGC",
         "((" + h + "." + h + "))", 1, 710},
        {"three helices, two nucleotides", "GCGGACAAAAGUCCAAGGACAAAAGUCCGC",
         "((" + h + ".." + h + "))", 1, 460},
        {"four helices, nothing unpaired", "GCGGACAAAAGUCCGGACAAAAGUCCGGACAAAAGUCCGC",
         "((" + h + h + h + "))", 1, 210},
        {"four helices, one nucleotide", "GCGGACAAAAGUCCGGACAAAAGUCCAGGACAAAAGUCCGC",
         "((" + h + h + "." + h + "))", 1, 260},
        {"five helices, nothing unpaired", "GCGGACAAAAGUCCGGACAAAAGUCCGGACAAAAGUCCGGACAAAAGUCCGC",
         "((" + h + h + h + h + "))", 1, 460},
        {"five helices, one nucleotide", "GCAGGACAAAAGUCCGGACAAAAGUCCGGACAAAAGUCCGGACAAAAGUCCGC",
         "((." + h + h + h + h + "))", 1, 330},
        // Asymmetry: 0.9 times the mean difference between neighbouring gaps, which is 0,
        // 2/3, 4/3 and 10/3 here; the mean counts up to 2.
        {"a nucleotide in every gap", "GCAGGACAAAAGUCCAGGACAAAAGUCCAGC",
         "((." + h + "." + h + ".))", 1, 210},
        {"gaps of 1, 0 and 1", "GCAGGACAAAAGUCCGGACAAAAGUCCAGC", "((." + h + h + ".))", 1, 400},
        {"gaps of 2, 1 and 0", "GCAGGGACAAAAGUCCAGGACAAAAGUCCGC", "((.." + h + "." + h + "))", 1,
         340},
        {"all the nucleotides in one gap", "GCAAAAAGGACAAAAGUCCGGACAAAAGUCCGC",
         "((....." + h + h + "))", 1, 520},
        // More than 8 unpaired nucleotides: 1.1 ln(n / 8), which is 0.1 for 9 and 1.0 for 20
        {"eight nucleotides", "GCAAAGGACAAAAGUCCAAGGACAAAAGUCCAAAGC",
         "((..." + h + ".." + h + "...))", 1, 300},
        {"nine nucleotides", "GCAAAGGACAAAAGUCCAAAGGACAAAAGUCCAAAGC",
         "((..." + h + "..." + h + "...))", 1, 250},
        {"twenty nucleotides", "GCAAAAAAAGGACAAAAGUCCAAAAAAGGACAAAAGUCCAAAAAAAGC",
         "((......." + h + "......" + h + ".......))", 1, 400},
        // Other pairs at the ends of the helices
        {"A-U closing pair", "GAGGACAAAAGUCCAGGACAAAAGUCCUC", "((" + h + "." + h + "))", 1, 790},
        {"U-G closing pair", "GUAGGACAAAAGUCCAGGACAAAAGUCCAGC", "((." + h + "." + h + ".))", 1,
         260},
        {"branches that end in A-U and G-U", "GCAGACAAAAGUCUAGGACAAAAGUCUAGC",
         "((" + h + "." + h + ".))", 1, 590},
        {"G and A around C-G branches", "GCGCGACAAAAGUCGAGCGACAAAAGUCGAGC",
         "((." + h + ".." + h + ".))", 1, 330},
        {"a branch of one pair", "GCACAAAAGAGGACAAAAGUCCGC", "((.(....)." + h + "))", 1, 460},
        {"five different helices", "GUCAGACAAAAGUCUGGCGGAAACGCCGACAAAAGUCGAAGGACAAAAGUCUUAC",
         "((.((((....)))).(((....)))((((....))))..((((....)))).))", 1, 90},
        {"a branch with a bulge and an internal loop", "GCAGCGAAGGACGUUCGCGAGCCUCGCGCGGAAACGCGGC",
         "((.((((.((.((....))..))))))(((....))).))", 1, 400},
        // A multibranch loop that is a branch of another (closures 2-42 and 4-28)
        {"loop around another loop", "GCAAGGACAAAAGUCCGCGGAAACGCAUAGGACAAAAGUCCGC",
         "((.(((((....))))(((....))).).((((....))))))", 1, 520},
        {"loop inside another loop", "GCAAGGACAAAAGUCCGCGGAAACGCAUAGGACAAAAGUCCGC",
         "((.(((((....))))(((....))).).((((....))))))", 3, 790},
        // N does not stack
        {"N in every gap", "GCNGGACAAAAGUCCNGGACAAAAGUCCNGC", "((." + h + "." + h + ".))", 1, 750},
        // A helix end RNAstructure does not allow: a G-A pair (efn2 prints 1408.8), N next to
        // the pair of a helix of one pair, inside a branch (1404.6) or outside the closing
        // pair (1407.1, closure 2-28)
        {"pair that cannot form", "GCGGACAAAAGUCAAGGACAAAAGUCCGC", "((" + h + "." + h + "))", 1,
         INF},
        {"N inside a branch of one pair", "GCACNAAAGAGGACAAAAGUCCGC", "((.(....)." + h + "))", 1,
         INF},
        {"N outside a closing helix of one pair", "NCGGACAAAAGUCCAGGACAAAAGUCCGA",
         ".(" + h + "." + h + ").", 1, INF},
    };
    for (const LoopCase& c : cases) {
        EXPECT_EQ(multibranch(c.sequence, c.structure, c.i, rsp), c.expected)
            << c.what << ": " << c.sequence;
    }

    // Only for multibranch loops: not a stack, a hairpin or the exterior loop
    EXPECT_THROW((void) multibranch("GGACAAAAGUCC", h, 0, rsp), knotergy::DetailedException);
    EXPECT_THROW((void) multibranch("GGACAAAAGUCC", h, 3, rsp), knotergy::DetailedException);
    const knotergy::ProcessedRNAEntry rna = process("GGACAAAAGUCC", h);
    knotergy::LoopFactory factory(rna);
    EXPECT_THROW(
        (void) RNAstructureFunctions::multibranch_energy(factory.get_root_node(), rna, rsp),
        knotergy::DetailedException);
}

// A pseudoknotted child counts as a helix of the loop, but it has no helix end and nothing
// stacks on it. efn2 has no rule for this case, so the expected values are put together from
// the parameters.
TEST(RNAstructureEnergies, PseudoknotsInMultibranchLoops) {
    const knotergy::rnastructure_param rsp = knotergy::RNAstructureParams::load();
    const auto& m = rsp.data.misc();

    const std::string knot = "GGCAAGCGAAGCCAACGC";
    const std::string knot_structure = "(((..[[[..)))..]]]";
    // Energy of the loop closed by the second pair, and the end of that pair seen from inside
    const auto loop_and_end = [&](const std::string& sequence, const std::string& structure) {
        const knotergy::ProcessedRNAEntry rna = process(sequence, structure);
        const std::size_t j = sequence.size() - 2;
        return std::pair<int, int>{multibranch(sequence, structure, 1, rsp),
                                   RNAstructureFunctions::helix_end(
                                       j, 1, rna.get_sequence(), rna.get_modified_sequence(), rsp)};
    };

    // Around one pseudoknot: two helices
    const auto [alone, alone_end] = loop_and_end("GC" + knot + "GC", "((" + knot_structure + "))");
    EXPECT_EQ(alone, 10 * (m.efn2a + 2 * m.efn2c) + alone_end);

    // Around two: three helices with nothing unpaired, which are strained
    const auto [two, two_end] =
        loop_and_end("GC" + knot + knot + "GC", "((" + knot_structure + knot_structure + "))");
    EXPECT_EQ(two, 10 * (m.efn2a + 3 * m.efn2c + m.strain) + two_end);

    // The nucleotides next to a pseudoknot can only stack on the closing pair.
    const std::string sequence = "GCA" + knot + "AGC";
    const std::string structure = "((." + knot_structure + ".))";
    const auto [flanked, flanked_end] = loop_and_end(sequence, structure);
    const knotergy::ProcessedRNAEntry rna = process(sequence, structure);
    const std::size_t j = sequence.size() - 2;
    const int stacking =
        std::min({0,
                  RNAstructureFunctions::dangle5_energy(1, j, rna.get_sequence(),
                                                        rna.get_modified_sequence(), rsp, true),
                  RNAstructureFunctions::dangle3_energy(1, j, rna.get_sequence(),
                                                        rna.get_modified_sequence(), rsp, true),
                  RNAstructureFunctions::mismatch_energy(
                      1, j, rna.get_sequence(), rna.get_modified_sequence(), rsp, false, true)});
    EXPECT_LT(stacking, 0);
    EXPECT_EQ(flanked, 10 * (m.efn2a + 2 * m.efn2c) + flanked_end + stacking);
}

// "Total energy" of efn2 -w: the exterior loop plus every loop inside the structure.
TEST(RNAstructureEnergies, TotalsMatchEfn2) {
    const knotergy::rnastructure_param rsp = knotergy::RNAstructureParams::load();

    const std::vector<StructureCase> cases = {
        {"one hairpin", "GGACAAAAGUCC", "((((....))))", -380},
        {"two hairpins", "AGGACAAAAGUCCAGCGGAAACGCG", ".((((....)))).(((....))).", -1020},
        {"a bulge and an internal loop", "CGCGAAGGACGUUCGCGAGCCUCGCAGGACAAAAGUCCAA",
         ".((((.((.((....))..)))))).((((....))))..", -1330},
        {"three branches", "GCGGAAACGCGCGAAGGACGUUCGCGAGCCUCGCUAGACAAAAGUCUC",
         "(((....)))((((.((.((....))..)))))).((((....)))).", -1350},
        {"one multibranch loop", "AGCGGACAAAAGUCCAGGACAAAAGUCCGCG",
         ".((((((....)))).((((....)))))).", -560},
        {"three helices with nothing unpaired", "GACGCGGAAACGCCGACAAAAGUCGGUC",
         "((((((....)))((((....)))))))", -290},
        {"a multibranch loop inside another",
         "GGCAAGGACAAAAGUCCGCGGAAACGCAUAGGACAAAAGUCCGCAGCGGAAACGC",
         ".((.(((((....))))(((....))).).((((....)))))).(((....)))", -690},
        {"two multibranch loops, a bulge and an internal loop",
         "GCGGACAAAAGUCCAAGCGGAAACGCAGCUCUGCGAAGGACGUUCGCGAGCCUCGCAAGACAAAAGUCUCAAAAGAGC",
         "((((((....))))..(((....))).)).((((((.((.((....))..)))))).((((....))))(....))).", -1180},
    };
    for (const StructureCase& c : cases) {
        EXPECT_EQ(total(c.sequence, c.structure, rsp), c.expected) << c.what << ": " << c.sequence;
    }
}

// efn2 -a m6A -w, efn2 -a PSU -w and efn2 -a 1 -w
TEST(RNAstructureEnergies, ExteriorLoopsWithModifiedBases) {
    const knotergy::rnastructure_param rsp =
        knotergy::RNAstructureParams::load("turner2004", {"m6A", "Ψ", "m1Ψ"});

    const std::string h = "((((....))))";
    const std::vector<StructureCase> cases = {
        {"m6A after a G-C helix", "GGACAAAAGUCC6", h + ".", -130},
        {"m6A on both sides of an A-U helix", "6AGACAAAAGUCU6", "." + h + ".", -80},
        {"m6A on both sides of a G-C helix", "6GGACAAAAGUCC6", "." + h + ".", -180},
        {"A on both sides of an m6A-U helix", "A6GACAAAAGUCUA", "." + h + ".", -130},
        {"flush stack on an m6A-U pair", "GGACAAAAGUCC6GACAAAAGUCU", h + h, -120},
        {"stack across m6A", "GGACAAAAGUCC6GGACAAAAGUCCA", h + "." + h + ".", -400},
        {"Ψ after a G-C helix", "GGACAAAAGUCCP", h + ".", -150},
        {"flush stack on an A-Ψ pair", "GGACAAAAGUCCAGACAAAAGUCP", h + h, -200},
        {"stack across Ψ", "GGACAAAAGUCCPGGACAAAAGUCCA", h + "." + h + ".", -400},
        {"m1Ψ after a G-C helix", "GGACAAAAGUCC1", h + ".", -120},
        {"flush stack on an A-m1Ψ pair", "GGACAAAAGUCCAGACAAAAGUC1", h + h, -260},
        {"stack across m1Ψ", "GGACAAAAGUCC1GGACAAAAGUCCA", h + "." + h + ".", -400},
    };
    for (const StructureCase& c : cases) {
        EXPECT_EQ(exterior(c.sequence, c.structure, rsp), c.expected)
            << c.what << ": " << c.sequence;
    }

    // Total energy = -8.1
    EXPECT_EQ(total("6GGACAAAAGUCC6G6GGA6ACUCA", ".((((....)))).(((....))).", rsp), -810);
}

// efn2 -a m6A -w, efn2 -a PSU -w and efn2 -a 1 -w
TEST(RNAstructureEnergies, MultibranchLoopsWithModifiedBases) {
    const knotergy::rnastructure_param rsp =
        knotergy::RNAstructureParams::load("turner2004", {"m6A", "Ψ", "m1Ψ"});

    const std::string h = "((((....))))";
    const std::string gaps = "((." + h + "." + h + ".))";
    const std::string between = "((" + h + "." + h + "))";
    const std::vector<LoopCase> cases = {
        {"m6A in every gap", "GC6GGACAAAAGUCC6GGACAAAAGUCC6GC", gaps, 1, 220},
        {"m6A next to the closing pair", "GC6GGACAAAAGUCCGGACAAAAGUCC6GC", "((." + h + h + ".))", 1,
         390},
        {"a branch that ends in m6A-U", "GCA6GACAAAAGUCUAGGACAAAAGUCCGC",
         "((." + h + "." + h + "))", 1, 440},
        {"m6A-U closing pair", "G6GGACAAAAGUCCAGGACAAAAGUCCUC", between, 1, 790},
        {"m6A between two helices", "GCGGACAAAAGUCC6GGACAAAAGUCCGC", between, 1, 750},
        {"Ψ in every gap", "GCPGGACAAAAGUCCPGGACAAAAGUCCPGC", gaps, 1, 170},
        {"a branch that ends in Ψ-A", "GCAPGACAAAAGUCAAGGACAAAAGUCCGC", "((." + h + "." + h + "))",
         1, 500},
        {"A-Ψ closing pair", "GAGGACAAAAGUCCAGGACAAAAGUCCPC", between, 1, 720},
        {"Ψ between two helices", "GCGGACAAAAGUCCPGGACAAAAGUCCGC", between, 1, 730},
        {"m1Ψ in every gap", "GC1GGACAAAAGUCC1GGACAAAAGUCC1GC", gaps, 1, 260},
        {"a branch that ends in m1Ψ-A", "GCA1GACAAAAGUCAAGGACAAAAGUCCGC",
         "((." + h + "." + h + "))", 1, 470},
        {"A-m1Ψ closing pair", "GAGGACAAAAGUCCAGGACAAAAGUCC1C", between, 1, 730},
        {"m1Ψ between two helices", "GCGGACAAAAGUCC1GGACAAAAGUCCGC", between, 1, 760},
    };
    for (const LoopCase& c : cases) {
        EXPECT_EQ(multibranch(c.sequence, c.structure, c.i, rsp), c.expected)
            << c.what << ": " << c.sequence;
    }

    // Total energy = -5.4, -6.9 and -6.5
    const std::string structure = ".((" + h + ".(((....))).)).";
    EXPECT_EQ(total("6GCGGACAAAAGUCC6G6GGA6ACUCAGCA", structure, rsp), -540);
    EXPECT_EQ(total("PGCGGACAAAAGUCCPGAGGAPACPCAGCA", structure, rsp), -690);
    EXPECT_EQ(total("1GCGGACAAAAGUCC1GAGGA1AC1CAGCA", structure, rsp), -650);
}

// Nucleotides that are already encoded go straight to ThermoRNA's functions, which take the
// codes of a loop instead of a sequence and positions, and from_rnastructure() converts the
// result. Here a sequence is encoded once and used for every loop.
TEST(RNAstructureEnergies, EncodedNucleotidesGiveTheSameEnergies) {
    using Element = rs::LoopElement;
    const knotergy::rnastructure_param rsp =
        knotergy::RNAstructureParams::load("turner2004", {"m6A", "Ψ", "m1Ψ"});
    const rs::Datatable& d = rsp.data;
    const auto convert = [](int energy) {
        return knotergy::RNAstructureUtils::from_rnastructure(energy);
    };

    //                             0         1         2
    //                             0123456789012345678901234567
    const std::string structure = "((((((....))))((((....))))))";
    // Some loops of the structure: the stack of 0-27 and 1-26, the hairpin closed by 5-10, the
    // stack of 2-13 and 3-12 as an internal loop without unpaired nucleotides, the multibranch
    // loop closed by 1-26 and the exterior loop.
    const auto energies = [&](const std::vector<rs::Code>& s) {
        return std::vector<int>{
            convert(d.stack(s[0], s[27], s[1], s[26])),
            convert(rs::hairpin_energy(4, &s[5], s[4], s[3], s[11], d)),
            convert(rs::internal_loop_energy(0, 0, s[2], s[13], s[3], s[12], s[3], s[12], s[2],
                                             s[13], s[1], s[14], s[4], s[11], 1, d)),
            convert(rs::multibranch_energy(s[1], s[26], s[0], s[27],
                                           {Element::helix(s[2], s[13], s[3], s[12]),
                                            Element::helix(s[14], s[25], s[15], s[24])},
                                           d)),
            convert(rs::exterior_energy({Element::helix(s[0], s[27], s[1], s[26])}, d))};
    };
    const auto from_the_sequence = [&](const std::string& sequence) {
        const knotergy::ProcessedRNAEntry rna = process(sequence, structure);
        return std::vector<int>{
            stack(rna, 0, 27, rsp), hairpin(rna, 5, 10, rsp), internal(sequence, structure, 2, rsp),
            multibranch(sequence, structure, 1, rsp), exterior(sequence, structure, rsp)};
    };

    // efn2 -w: Stack = -3.4, Hairpin = 4.1, Stack = -3.3, Multiloop = 8.2, Exterior loop = 0.0
    const std::string plain = "GCGGACAAAAGUCCGGACAAAAGUCCGC";
    EXPECT_EQ(energies(d.encode(plain)), (std::vector<int>{-340, 410, -330, 820, 0}));
    EXPECT_EQ(from_the_sequence(plain), energies(d.encode(plain)));

    // With modified bases, the symbols to encode are those of the modified sequence. Here an
    // m6A-U pair closes the multibranch loop and a pseudouridine follows the first hairpin.
    const std::string modified = "G6GGACAAAAGPCCGGACAAAAGUCCUC";
    const knotergy::ProcessedRNAEntry rna = process(modified, structure);
    ASSERT_TRUE(rna.has_modified_bases());
    EXPECT_EQ(from_the_sequence(modified), energies(d.encode(rna.get_modified_sequence())));
    EXPECT_NE(from_the_sequence(modified), from_the_sequence("GAGGACAAAAGUCCGGACAAAAGUCCUC"));
}

// A symbol RNAstructure's sets do not know is an error that names the symbol, with or without
// modified bases in the sequence.
// The two arguments of multibranch_energy and external_energy that leave efn2's default:
// the simple multibranch rule (efn2 --simple), and loops without coaxial stacking. efn2 has no
// option for the second. Its expected values are what RNAstructure's own folding recursion
// gives a structure when its coaxial stacking is disabled (as with Fold --disablecoax) and the
// structure is the only one it is allowed to find; that recursion uses the simple rule.
TEST(RNAstructureEnergies, SimpleRuleAndNoCoaxialStacking) {
    const knotergy::rnastructure_param rsp = knotergy::RNAstructureParams::load();
    const rs::MultibranchRule simple = rs::MultibranchRule::Simple;
    const rs::CoaxialStacking off = rs::CoaxialStacking::Off;

    // Three helices with nothing unpaired: efn2 -w has Multiloop = 8.2, efn2 -w -s 5.1. Without
    // coaxial stacking nothing is left that can stack: the offset and the penalty per helix.
    const std::string sequence = "GCGGACAAAAGUCCGGACAAAAGUCCGC";
    const std::string structure = "((((((....))))((((....))))))";
    EXPECT_EQ(multibranch(sequence, structure, 1, rsp), 820);
    EXPECT_EQ(multibranch(sequence, structure, 1, rsp, simple), 510);
    EXPECT_EQ(multibranch(sequence, structure, 1, rsp, simple, off), 750);
    // The rule and the stacking are independent: efn2's rule adds the strain (3.1) either way.
    EXPECT_EQ(multibranch(sequence, structure, 1, rsp, rs::MultibranchRule::Efn2, off), 1060);

    // Two helices one nucleotide apart: efn2 stacks them across it (Exterior loop = -3.7).
    // Without coaxial stacking each takes the A after it as a 3' dangle, twice -1.7.
    const std::string two = "GGACAAAAGUCCAGGACAAAAGUCCA";
    EXPECT_EQ(exterior(two, "((((....)))).((((....)))).", rsp), -370);
    EXPECT_EQ(exterior(two, "((((....)))).((((....)))).", rsp, off), -340);

    // The structures of TotalsMatchEfn2: their totals from efn2 -s, and from RNAstructure's
    // folding recursion without coaxial stacking.
    struct Case {
        const char* what;
        std::string sequence;
        std::string structure;
        int simple;
        int folding;
    };
    const std::vector<Case> cases = {
        {"one hairpin", "GGACAAAAGUCC", "((((....))))", -380, -380},
        {"two hairpins", "AGGACAAAAGUCCAGCGGAAACGCG", ".((((....)))).(((....))).", -1020, -990},
        {"a bulge and an internal loop", "CGCGAAGGACGUUCGCGAGCCUCGCAGGACAAAAGUCCAA",
         ".((((.((.((....))..)))))).((((....))))..", -1330, -1300},
        {"three branches", "GCGGAAACGCGCGAAGGACGUUCGCGAGCCUCGCUAGACAAAAGUCUC",
         "(((....)))((((.((.((....))..)))))).((((....)))).", -1350, -1190},
        {"one multibranch loop", "AGCGGACAAAAGUCCAGGACAAAAGUCCGCG",
         ".((((((....)))).((((....)))))).", -930, -690},
        {"three helices with nothing unpaired", "GACGCGGAAACGCCGACAAAAGUCGGUC",
         "((((((....)))((((....)))))))", -600, -270},
        {"a multibranch loop inside another",
         "GGCAAGGACAAAAGUCCGCGGAAACGCAUAGGACAAAAGUCCGCAGCGGAAACGC",
         ".((.(((((....))))(((....))).).((((....)))))).(((....)))", -1120, -620},
        {"two multibranch loops, a bulge and an internal loop",
         "GCGGACAAAAGUCCAAGCGGAAACGCAGCUCUGCGAAGGACGUUCGCGAGCCUCGCAAGACAAAAGUCUCAAAAGAGC",
         "((((((....))))..(((....))).)).((((((.((.((....))..)))))).((((....))))(....))).", -1350,
         -760},
    };
    for (const Case& c : cases) {
        EXPECT_EQ(total(c.sequence, c.structure, rsp, simple), c.simple)
            << c.what << ": " << c.sequence;
        EXPECT_EQ(total(c.sequence, c.structure, rsp, simple, off), c.folding)
            << c.what << ": " << c.sequence;
    }
}

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

TEST(EnergyFunctions, HairpinFollowsTheSelectedModel) {
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

    knotergy::AllParams vienna_params(mp, pkp, vp, rsp, knotergy::PKFreeEnergyModel::ViennaRNA2);
    knotergy::AllParams rs_params(mp, pkp, vp, rsp, knotergy::PKFreeEnergyModel::efn2);
    knotergy::EnergyFunctions vienna_energy(vienna_params);
    knotergy::EnergyFunctions rs_energy(rs_params);

    // ViennaRNA with unmodified input
    EXPECT_EQ(vienna_energy.hairpin(3, 8, sequence, {}, is_inf),
              knotergy::ViennaFunctions::hairpin_energy(3, 8, sequence, is_inf, vp));

    // ViennaRNA with the JSON modified-base corrections
    EXPECT_EQ(vienna_energy.hairpin(3, 8, sequence, mod_sequence, is_inf),
              knotergy::ModHairpin::find_mod_hairpin_energy(3, 8, sequence, mod_sequence, vp, mp,
                                                            is_inf));

    // RNAstructure selected explicitly
    EXPECT_EQ(rs_energy.hairpin(3, 8, sequence, mod_sequence, is_inf), 410);
    EXPECT_FALSE(is_inf);
}

TEST(EnergyFunctions, InternalFollowsTheSelectedModel) {
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

    knotergy::AllParams vienna_params(mp, pkp, vp, rsp, knotergy::PKFreeEnergyModel::ViennaRNA2);
    knotergy::AllParams rs_params(mp, pkp, vp, rsp, knotergy::PKFreeEnergyModel::efn2);
    knotergy::EnergyFunctions vienna_energy(vienna_params);
    knotergy::EnergyFunctions rs_energy(rs_params);

    // ViennaRNA with unmodified input
    EXPECT_EQ(vienna_energy.internal(3, 16, 5, 14, sequence, {}),
              knotergy::ViennaFunctions::internal_loop_energy(3, 16, 5, 14, sequence, vp));

    // ViennaRNA with the JSON modified-base corrections
    EXPECT_EQ(vienna_energy.internal(3, 16, 5, 14, sequence, mod_sequence),
              knotergy::ModInternal::find_mod_internal_energy(3, 16, 5, 14, mod_sequence, sequence,
                                                              vp, mp));

    // RNAstructure selected explicitly
    EXPECT_EQ(rs_energy.internal(3, 16, 5, 14, sequence, mod_sequence), 90);
}

TEST(EnergyFunctions, DanglesFollowTheSelectedModel) {
    knotergy::vrna_md_param vp = knotergy::ViennaParams::load_energy_parameters(turner_file, 2);
    const knotergy::pk_param pkp =
        knotergy::PseudoknotParams::load_pk_param(pkp_file, knotergy::RoundMethod::None);
    std::vector<knotergy::modified_base_param> json_params =
        knotergy::ModParams::load_modified_energy_parameters(mod_folder);
    const knotergy::all_mod_params mp{json_params};
    const knotergy::rnastructure_param rsp =
        knotergy::RNAstructureParams::load("turner2004", {"m6A"});

    // m6A on both sides of an A-U pair: outside it, at the end of a helix (pair 1-12), and
    // inside it, where it closes a multibranch loop (pair 1-32). The JSON file has both dangles.
    const std::string helix = ".((((....)))).";
    const std::string multiloop = "((..((((....))))..((((....))))..))";
    const knotergy::ProcessedRNAEntry rna = process("6AGACAAAAGUCU6", helix);
    const knotergy::ProcessedRNAEntry closed =
        process("GU6NGGACAAAAGUCCNNGGACAAAAGUCCN6AC", multiloop);
    // The same sequences without modified bases, as ViennaRNA alone sees them
    const knotergy::ProcessedRNAEntry plain_rna =
        knotergy::RNAProcessor::process_rna(knotergy::RNAEntry(rna.get_sequence(), helix));
    const knotergy::ProcessedRNAEntry plain_closed =
        knotergy::RNAProcessor::process_rna(knotergy::RNAEntry(closed.get_sequence(), multiloop));
    const knotergy::all_mod_params no_mods;

    knotergy::AllParams vienna_params(mp, pkp, vp, rsp, knotergy::PKFreeEnergyModel::ViennaRNA2);
    knotergy::AllParams rs_params(mp, pkp, vp, rsp, knotergy::PKFreeEnergyModel::efn2);
    knotergy::EnergyFunctions vienna_energy(vienna_params);
    knotergy::EnergyFunctions rs_energy(rs_params);
    const auto d5 = [&](knotergy::EnergyFunctions& energy, const knotergy::ProcessedRNAEntry& entry,
                        std::size_t i, std::size_t j, bool is_closing) {
        return energy.d5(i, j, entry.get_sequence(), entry.get_modified_sequence(), is_closing);
    };
    const auto d3 = [&](knotergy::EnergyFunctions& energy, const knotergy::ProcessedRNAEntry& entry,
                        std::size_t i, std::size_t j, bool is_closing) {
        return energy.d3(i, j, entry.get_sequence(), entry.get_modified_sequence(), is_closing);
    };

    // ViennaRNA with unmodified input: its tables, as the modified-base code reads them when
    // nothing is modified
    const int vienna_d5 = d5(vienna_energy, plain_rna, 1, 12, false);
    const int vienna_d3 = d3(vienna_energy, plain_rna, 1, 12, false);
    EXPECT_EQ(vienna_d5, knotergy::ModBaseUtils::get_dangle5_energy(1, 12, plain_rna, vp, no_mods));
    EXPECT_EQ(vienna_d3, knotergy::ModBaseUtils::get_dangle3_energy(1, 12, plain_rna, vp, no_mods));
    EXPECT_EQ(d5(vienna_energy, plain_closed, 1, 32, true),
              knotergy::ModBaseUtils::get_dangle5_energy(1, 32, plain_closed, vp, no_mods, true));
    EXPECT_EQ(d3(vienna_energy, plain_closed, 1, 32, true),
              knotergy::ModBaseUtils::get_dangle3_energy(1, 32, plain_closed, vp, no_mods, true));

    // ViennaRNA with the JSON modified-base parameters
    EXPECT_EQ(d5(vienna_energy, rna, 1, 12, false),
              knotergy::ModBaseUtils::get_dangle5_energy(1, 12, rna, vp, mp));
    EXPECT_EQ(d3(vienna_energy, rna, 1, 12, false),
              knotergy::ModBaseUtils::get_dangle3_energy(1, 12, rna, vp, mp));
    EXPECT_EQ(d5(vienna_energy, closed, 1, 32, true),
              knotergy::ModBaseUtils::get_dangle5_energy(1, 32, closed, vp, mp, true));
    EXPECT_EQ(d3(vienna_energy, closed, 1, 32, true),
              knotergy::ModBaseUtils::get_dangle3_energy(1, 32, closed, vp, mp, true));
    EXPECT_NE(d5(vienna_energy, rna, 1, 12, false), vienna_d5);
    EXPECT_NE(d3(vienna_energy, rna, 1, 12, false), vienna_d3);

    // RNAstructure selected explicitly (efn2 -a m6A -w, see DanglesWithModifiedBases)
    EXPECT_EQ(d5(rs_energy, rna, 1, 12, false), -90);
    EXPECT_EQ(d3(rs_energy, rna, 1, 12, false), -100);
    EXPECT_EQ(d5(rs_energy, closed, 1, 32, true), -90);
    EXPECT_EQ(d3(rs_energy, closed, 1, 32, true), -100);
}

TEST(EnergyFunctions, MismatchAndTerminalFollowTheSelectedModel) {
    knotergy::vrna_md_param vp = knotergy::ViennaParams::load_energy_parameters(turner_file, 2);
    const knotergy::pk_param pkp =
        knotergy::PseudoknotParams::load_pk_param(pkp_file, knotergy::RoundMethod::None);
    std::vector<knotergy::modified_base_param> json_params =
        knotergy::ModParams::load_modified_energy_parameters(mod_folder);
    const knotergy::all_mod_params mp{json_params};
    const knotergy::rnastructure_param rsp =
        knotergy::RNAstructureParams::load("turner2004", {"m6A"});

    // m6A on both sides of a G-C pair: at the end of a helix (pair 1-12) and where it closes a
    // multibranch loop (pair 1-32). And a helix that ends in an m6A-U pair (1-12).
    const std::string helix = ".((((....)))).";
    const std::string multiloop = "((..((((....))))..((((....))))..))";
    const knotergy::ProcessedRNAEntry rna = process("6GGACAAAAGUCC6", helix);
    const knotergy::ProcessedRNAEntry closed =
        process("GC6NGGACAAAAGUCCNNGGACAAAAGUCCN6GC", multiloop);
    const knotergy::ProcessedRNAEntry end = process("A6GACAAAAGUCUA", helix);
    // The same sequences without modified bases, as ViennaRNA alone sees them
    const knotergy::ProcessedRNAEntry plain_rna =
        knotergy::RNAProcessor::process_rna(knotergy::RNAEntry(rna.get_sequence(), helix));
    const knotergy::ProcessedRNAEntry plain_closed =
        knotergy::RNAProcessor::process_rna(knotergy::RNAEntry(closed.get_sequence(), multiloop));
    const knotergy::all_mod_params no_mods;

    knotergy::AllParams vienna_params(mp, pkp, vp, rsp, knotergy::PKFreeEnergyModel::ViennaRNA2);
    knotergy::AllParams rs_params(mp, pkp, vp, rsp, knotergy::PKFreeEnergyModel::efn2);
    knotergy::EnergyFunctions vienna_energy(vienna_params);
    knotergy::EnergyFunctions rs_energy(rs_params);
    const auto mismatch = [&](knotergy::EnergyFunctions& energy,
                              const knotergy::ProcessedRNAEntry& entry, std::size_t i,
                              std::size_t j, bool is_external, bool is_closing) {
        return energy.mismatch(i, j, entry.get_sequence(), entry.get_modified_sequence(),
                               is_external, is_closing);
    };
    const auto terminal = [&](knotergy::EnergyFunctions& energy, bool is_closing,
                              bool modified = true) {
        return energy.terminal(
            1, 12, end.get_sequence(),
            modified ? end.get_modified_sequence() : std::vector<std::string_view>{}, is_closing);
    };

    // ViennaRNA with unmodified input: its tables, as the modified-base code reads them when
    // nothing is modified
    const int vienna_mismatch = mismatch(vienna_energy, plain_rna, 1, 12, true, false);
    EXPECT_EQ(vienna_mismatch,
              knotergy::ModBaseUtils::get_mismatch_energy(1, 12, plain_rna, vp, no_mods, true));
    EXPECT_EQ(mismatch(vienna_energy, plain_rna, 1, 12, false, false),
              knotergy::ModBaseUtils::get_mismatch_energy(1, 12, plain_rna, vp, no_mods, false));
    EXPECT_EQ(
        mismatch(vienna_energy, plain_closed, 1, 32, false, true),
        knotergy::ModBaseUtils::get_mismatch_energy(1, 32, plain_closed, vp, no_mods, false, true));
    EXPECT_EQ(terminal(vienna_energy, false, false),
              vp.p->TerminalAU);  // the unmodified pair is A-U
    EXPECT_EQ(terminal(vienna_energy, true, false), vp.p->TerminalAU);

    // ViennaRNA with the JSON modified-base parameters
    EXPECT_EQ(mismatch(vienna_energy, rna, 1, 12, true, false),
              knotergy::ModBaseUtils::get_mismatch_energy(1, 12, rna, vp, mp, true));
    EXPECT_EQ(mismatch(vienna_energy, rna, 1, 12, false, false),
              knotergy::ModBaseUtils::get_mismatch_energy(1, 12, rna, vp, mp, false));
    EXPECT_EQ(mismatch(vienna_energy, closed, 1, 32, false, true),
              knotergy::ModBaseUtils::get_mismatch_energy(1, 32, closed, vp, mp, false, true));
    EXPECT_NE(mismatch(vienna_energy, rna, 1, 12, true, false), vienna_mismatch);
    EXPECT_EQ(terminal(vienna_energy, false),
              0);  // the JSON file has no penalty for the end of an m6A-U pair
    EXPECT_EQ(terminal(vienna_energy, true), 0);

    // RNAstructure selected explicitly. efn2 -a m6A -w: Exterior loop = -1.8 with
    // the m6As and 0.0 with N; Multiloop = 5.7 and 7.5; Exterior loop = 0.0 for N6GACAAAAGUCUN.
    EXPECT_EQ(mismatch(rs_energy, rna, 1, 12, true, false), -180);
    EXPECT_EQ(mismatch(rs_energy, closed, 1, 32, false, true), -180);
    EXPECT_EQ(terminal(rs_energy, false), 0);
    // Without m6A, the end of an A-U helix costs 0.5 (Exterior loop = 0.5 for AGACAAAAGUCU).
    const knotergy::ProcessedRNAEntry au = process("AGACAAAAGUCU", "((((....))))");
    EXPECT_EQ(rs_energy.terminal(0, 11, au.get_sequence(), au.get_modified_sequence()), 50);
}

TEST(EnergyFunctions, MultibranchFollowsTheSelectedModel) {
    knotergy::vrna_md_param vp = knotergy::ViennaParams::load_energy_parameters(turner_file, 2);
    const knotergy::pk_param pkp =
        knotergy::PseudoknotParams::load_pk_param(pkp_file, knotergy::RoundMethod::None);
    std::vector<knotergy::modified_base_param> json_params =
        knotergy::ModParams::load_modified_energy_parameters(mod_folder);
    const knotergy::all_mod_params mp{json_params};
    const knotergy::rnastructure_param rsp =
        knotergy::RNAstructureParams::load("turner2004", {"m6A"});

    // m6A in every gap of a loop of three helices (efn2 -a m6A -w: Multiloop = 2.2)
    const std::string structure = "((.((((....)))).((((....)))).))";
    const knotergy::ProcessedRNAEntry rna = process("GC6GGACAAAAGUCC6GGACAAAAGUCC6GC", structure);
    knotergy::LoopFactory factory(rna);
    const knotergy::LoopNode& loop = *factory.get_node_table()[1];
    // The same sequence without modified bases, as ViennaRNA alone sees it
    const knotergy::ProcessedRNAEntry plain_rna =
        knotergy::RNAProcessor::process_rna(knotergy::RNAEntry(rna.get_sequence(), structure));
    knotergy::LoopFactory plain_factory(plain_rna);
    const knotergy::LoopNode& plain_loop = *plain_factory.get_node_table()[1];

    knotergy::AllParams vienna_params(mp, pkp, vp, rsp, knotergy::PKFreeEnergyModel::ViennaRNA2);
    knotergy::AllParams rs_params(mp, pkp, vp, rsp, knotergy::PKFreeEnergyModel::efn2);
    knotergy::EnergyFunctions vienna_energy(vienna_params);
    knotergy::EnergyFunctions rs_energy(rs_params);

    // ViennaRNA with unmodified input
    const int vienna = vienna_energy.multibranch(plain_loop, plain_rna);
    EXPECT_EQ(vienna, knotergy::ViennaFunctions::multibranch_energy(plain_loop, plain_rna, vp));

    // ViennaRNA with the JSON modified-base corrections
    EXPECT_EQ(vienna_energy.multibranch(loop, rna),
              knotergy::ModMultiloop::find_mod_multiloop_energy(loop, rna, vp, mp));
    EXPECT_NE(vienna_energy.multibranch(loop, rna), vienna);

    // RNAstructure selected explicitly
    EXPECT_EQ(rs_energy.multibranch(loop, rna), 220);
}

TEST(EnergyFunctions, ExternalFollowsTheSelectedModel) {
    knotergy::vrna_md_param vp = knotergy::ViennaParams::load_energy_parameters(turner_file, 2);
    const knotergy::pk_param pkp =
        knotergy::PseudoknotParams::load_pk_param(pkp_file, knotergy::RoundMethod::None);
    std::vector<knotergy::modified_base_param> json_params =
        knotergy::ModParams::load_modified_energy_parameters(mod_folder);
    const knotergy::all_mod_params mp{json_params};
    const knotergy::rnastructure_param rsp =
        knotergy::RNAstructureParams::load("turner2004", {"m6A"});

    // m6A on both sides of a G-C helix (efn2 -a m6A -w: Exterior loop = -1.8)
    const std::string structure = ".((((....)))).";
    const knotergy::ProcessedRNAEntry rna = process("6GGACAAAAGUCC6", structure);
    knotergy::LoopFactory factory(rna);
    const std::vector<knotergy::LoopNode*>& children = factory.get_root_node().children;
    // The same sequence without modified bases, as ViennaRNA alone sees it
    const knotergy::ProcessedRNAEntry plain_rna =
        knotergy::RNAProcessor::process_rna(knotergy::RNAEntry(rna.get_sequence(), structure));
    knotergy::LoopFactory plain_factory(plain_rna);

    knotergy::AllParams vienna_params(mp, pkp, vp, rsp, knotergy::PKFreeEnergyModel::ViennaRNA2);
    knotergy::AllParams rs_params(mp, pkp, vp, rsp, knotergy::PKFreeEnergyModel::efn2);
    knotergy::EnergyFunctions vienna_energy(vienna_params);
    knotergy::EnergyFunctions rs_energy(rs_params);

    // ViennaRNA with unmodified input
    const int vienna = vienna_energy.external(plain_factory.get_root_node().children, plain_rna);
    EXPECT_EQ(vienna, knotergy::ViennaFunctions::external_energy(
                          plain_factory.get_root_node().children, plain_rna, vp));

    // ViennaRNA with the JSON modified-base corrections
    EXPECT_EQ(vienna_energy.external(children, rna),
              knotergy::ModExternal::find_mod_external_energy(children, rna, vp, mp));
    EXPECT_NE(vienna_energy.external(children, rna), vienna);

    // RNAstructure selected explicitly
    EXPECT_EQ(rs_energy.external(children, rna), -180);
}

// efn2 --simple: select the simple multibranch rule via the model enum.
// ViennaRNA dangle settings do not change RNAstructure energies.
TEST(EnergyFunctions, SimpleMultibranchRuleWithRNAstructureParameters) {
    const knotergy::pk_param pkp =
        knotergy::PseudoknotParams::load_pk_param(pkp_file, knotergy::RoundMethod::None);
    const knotergy::all_mod_params mp{};
    const knotergy::rnastructure_param rsp = knotergy::RNAstructureParams::load();

    // efn2 -w has Multiloop = 8.2 for this loop and efn2 -w -s 5.1; Exterior loop = -3.7 for
    // the two helices below with either.
    const knotergy::ProcessedRNAEntry junction =
        process("GCGGACAAAAGUCCGGACAAAAGUCCGC", "((((((....))))((((....))))))");
    knotergy::LoopFactory junction_factory(junction);
    const knotergy::LoopNode& loop = *junction_factory.get_node_table()[1];
    const knotergy::ProcessedRNAEntry helices =
        process("GGACAAAAGUCCAGGACAAAAGUCCA", "((((....)))).((((....)))).");
    knotergy::LoopFactory helices_factory(helices);
    const std::vector<knotergy::LoopNode*>& children = helices_factory.get_root_node().children;

    for (const int dangles : {0, 1, 2, 3}) {
        knotergy::vrna_md_param vp =
            knotergy::ViennaParams::load_energy_parameters(turner_file, dangles);
        knotergy::AllParams normal_params(mp, pkp, vp, rsp, knotergy::PKFreeEnergyModel::efn2);
        knotergy::AllParams simple_params(mp, pkp, vp, rsp,
                                          knotergy::PKFreeEnergyModel::enf2Simple);
        knotergy::EnergyFunctions normal(normal_params);
        knotergy::EnergyFunctions simple(simple_params);
        EXPECT_EQ(normal.multibranch(loop, junction), 820) << dangles;
        EXPECT_EQ(normal.external(children, helices), -370) << dangles;
        EXPECT_EQ(simple.multibranch(loop, junction), 510) << dangles;
        EXPECT_EQ(simple.external(children, helices), -370) << dangles;
    }
}
