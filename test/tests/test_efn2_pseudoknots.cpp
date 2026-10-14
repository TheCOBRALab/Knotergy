#include "energy/efn2/Efn2Pseudoknot.hpp"
#include "io/input/CliArgs.hpp"
#include "test_utils.hpp"

#include <algorithm>
#include <initializer_list>
#include <utility>

namespace {

const std::string h3 = "(((....[[[....)))....]]]";

knotergy::ProcessedRNAEntry process_structure(const std::string& structure) {
    const auto pairs = knotergy::RNAProcessor::compute_pair_table(structure);
    std::string sequence(structure.size(), 'A');
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        if (pairs[i] != knotergy::NULL_INDEX && pairs[i] > i) {
            sequence[i] = 'G';
            sequence[pairs[i]] = 'C';
        }
    }
    return knotergy::RNAProcessor::process_rna(knotergy::RNAEntry(sequence, structure));
}

std::vector<knotergy::Efn2PkFeatures> features_for(const std::string& structure) {
    const auto rna = process_structure(structure);
    const auto original_pairs = rna.get_pair_table();
    knotergy::LoopFactory factory(rna);
    std::vector<const knotergy::LoopNode*> stack{&factory.get_root_node()};
    std::vector<knotergy::Efn2PkFeatures> features;
    while (!stack.empty()) {
        const auto* node = stack.back();
        stack.pop_back();
        if (node->loop_type == knotergy::LoopType::Pseudoknot) {
            features.push_back(knotergy::get_efn2_features(*node));
        }
        for (auto it = node->children.rbegin(); it != node->children.rend(); ++it) {
            stack.push_back(*it);
        }
    }
    EXPECT_EQ(rna.get_pair_table(), original_pairs);
    return features;
}

void expect_features(const knotergy::Efn2PkFeatures& actual, int ss, int ne,
                     std::initializer_list<std::pair<std::size_t, int>> spans) {
    std::array<int, 16> expected{};
    for (const auto& [length, count] : spans) expected[length] = count;
    EXPECT_EQ(actual.ss, ss);
    EXPECT_EQ(actual.ne, ne);
    EXPECT_EQ(actual.il, expected);
}

// Geometry and penalty goldens were obtained from RNAstructure 6.6's
// DistanceCounter/pkEnergyCalc, with its ignored length-zero repeat bucket omitted.
TEST(Efn2Pseudoknots, HTypeGeometryAndPenalty) {
    const auto features = features_for(h3);
    ASSERT_EQ(features.size(), 1u);
    expect_features(features[0], 12, 0, {{3, 2}});
    EXPECT_NEAR(knotergy::efn2_pk_penalty(features[0]), 675.6566515239583, 1e-8);
}

TEST(Efn2Pseudoknots, SingleBaseBulgesOnEitherArmStayInTheHelix) {
    for (const auto& structure : {"((.(....[[[....)))....]]]", "(((....[[[....)).)....]]]"}) {
        SCOPED_TRACE(structure);
        const auto features = features_for(structure);
        ASSERT_EQ(features.size(), 1u);
        expect_features(features[0], 12, 0, {{3, 2}});
    }
}

TEST(Efn2Pseudoknots, LargerBulgesSplitTheHelixAndCountTheirBases) {
    for (const auto& structure : {"((..(....[[[....)))....]]]", "(((....[[[....))..)....]]]"}) {
        SCOPED_TRACE(structure);
        const auto features = features_for(structure);
        ASSERT_EQ(features.size(), 1u);
        expect_features(features[0], 14, 0, {{1, 1}, {2, 1}, {3, 1}});
        EXPECT_NEAR(knotergy::efn2_pk_penalty(features[0]), 681.0519253179138, 1e-8);
    }
}

TEST(Efn2Pseudoknots, NestedHairpinInteriorIsSkipped) {
    const auto features = features_for("(((.((...))..[[[....)))....]]]");
    ASSERT_EQ(features.size(), 1u);
    expect_features(features[0], 11, 1, {{3, 2}});
    EXPECT_NEAR(knotergy::efn2_pk_penalty(features[0]), 686.4302888569793, 1e-8);
}

TEST(Efn2Pseudoknots, ABranchWithinABandSplitsTheHelix) {
    const auto features = features_for("((.((...)).(....[[[....)))....]]]");
    ASSERT_EQ(features.size(), 1u);
    expect_features(features[0], 14, 1, {{1, 1}, {2, 1}, {3, 1}});
}

TEST(Efn2Pseudoknots, KissingHelixIsCountedOnceAcrossBothArms) {
    const auto features = features_for("(((...[[[...)))...(((...]]]...)))");
    ASSERT_EQ(features.size(), 1u);
    expect_features(features[0], 15, 0, {{3, 3}});
    EXPECT_NEAR(knotergy::efn2_pk_penalty(features[0]), 709.8219078473745, 1e-8);
}

TEST(Efn2Pseudoknots, LongHelicesUseTheLastBucket) {
    const auto features =
        features_for(std::string(16, '(') + "....[[[...." + std::string(16, ')') + "....]]]");
    ASSERT_EQ(features.size(), 1u);
    expect_features(features[0], 12, 0, {{3, 1}, {15, 1}});
}

TEST(Efn2Pseudoknots, SinglePairSegmentsHaveNoSpanningPenalty) {
    const auto features = features_for("(....[....)....]");
    ASSERT_EQ(features.size(), 1u);
    expect_features(features[0], 12, 0, {{1, 2}});
    EXPECT_NEAR(knotergy::efn2_pk_penalty(features[0]), 217.997885126005, 1e-8);
}

TEST(Efn2Pseudoknots, ZeroContourUsesTheFiniteFallback) {
    const auto features = features_for("((([[[)))]]]");
    ASSERT_EQ(features.size(), 1u);
    expect_features(features[0], 0, 0, {{3, 2}});
    EXPECT_NEAR(knotergy::efn2_pk_penalty(features[0]), 492.6587663979533, 1e-8);
}

TEST(Efn2Pseudoknots, NestedPseudoknotUsesRetainedAndRemovedChildPairs) {
    const auto features = features_for("(((." + h3 + "...[[[....)))....]]]");
    ASSERT_EQ(features.size(), 2u);
    expect_features(features[0], 19, 1, {{3, 2}});
    expect_features(features[1], 12, 0, {{3, 2}});
}

TEST(Efn2Pseudoknots, NestedPseudoknotWithinABandSplitsItsParentHelix) {
    const auto features = features_for("((." + h3 + "...(....[[[....)))....]]]");
    ASSERT_EQ(features.size(), 2u);
    expect_features(features[0], 23, 1, {{1, 1}, {2, 1}, {3, 1}});
    expect_features(features[1], 12, 0, {{3, 2}});
}

TEST(Efn2Pseudoknots, SeparateClosedRegionsHaveIndependentPenalties) {
    const auto features = features_for(h3 + "...." + h3);
    ASSERT_EQ(features.size(), 2u);
    EXPECT_NEAR(knotergy::efn2_pk_penalty(features[0]) + knotergy::efn2_pk_penalty(features[1]),
                1351.3133030479166, 1e-8);
}

TEST(Efn2Pseudoknots, LinksIncludeBothEndpointsAndInnermostPairs) {
    for (const auto& structure : {h3, std::string("(....[....)....]")}) {
        const auto rna = process_structure(structure);
        knotergy::LoopFactory factory(rna);
        const auto& node = *factory.get_root_node().children.front();
        for (const auto& band : node.bands) {
            const auto& pairs = band.base_pairs();
            ASSERT_FALSE(pairs.empty());
            EXPECT_EQ(pairs.front().parent, nullptr);
            EXPECT_EQ(pairs.back().next, nullptr);
            for (std::size_t i = 0; i + 1 < pairs.size(); ++i) {
                EXPECT_EQ(pairs[i].next, &pairs[i + 1]);
                EXPECT_EQ(pairs[i + 1].parent, &pairs[i]);
            }
        }
    }
}

TEST(Efn2Pseudoknots, EnergyModeUsesThePenaltyAndUnscaledTurnerStacks) {
    const auto rna = process_structure(h3);
    knotergy::LoopFactory factory(rna);
    auto vp = knotergy::ViennaParams::load_energy_parameters(turner_file, 0, 37.0, 1.021, true);
    const knotergy::pk_param dp_params;
    knotergy::ComputeEnergy energy(factory.get_root_node(), rna, vp, dp_params, {}, false, false,
                                   knotergy::VerbosityLevel::Quiet,
                                   knotergy::PseudoknotModel::Efn2Polymer);
    EXPECT_NEAR(energy.getEnergy(), -6.443433484760417, 1e-8);
    const auto& node = *factory.get_root_node().children.front();
    ASSERT_TRUE(node.efn2_features.has_value());
    expect_features(*node.efn2_features, 12, 0, {{3, 2}});
}

TEST(Efn2Pseudoknots, CommandLineSelectsTheModel) {
    char program[] = "Knotergy";
    char flag[] = "--pk-model";
    char model[] = "efn2";
    char* argv[] = {program, flag, model};
    knotergy::CliArgs args;
    EXPECT_EQ(knotergy::CliArgs::parse(3, argv, args), knotergy::ParseStatus::Ok);
    EXPECT_EQ(args.pk_model, knotergy::PseudoknotModel::Efn2Polymer);
}

}  // namespace
