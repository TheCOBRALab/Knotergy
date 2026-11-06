// Loading efn2's (RNAstructure's) parameter tables through ThermoRNA, from Knotergy's
// params/RNAstructure (one folder per parameter set).
//
// MatchesRNAstructureDataTables also compares them with an RNAstructure installation when
// $DATAPATH (or $THERMORNA_RNASTRUCTURE_DATA) names its data_tables folder.

#include "test_utils.hpp"

#include <RNAstructure/params/parameter_set.hpp>
#include <io/parameters/RNAstructureParams.hpp>

#include <cstdlib>
#include <random>
#include <string>
#include <vector>

namespace rs = thermorna::rnastructure;

namespace {

std::string rnastructure_data() {
    for (const char* name : {"THERMORNA_RNASTRUCTURE_DATA", "DATAPATH"}) {
        const char* value = std::getenv(name);
        if (value != nullptr && *value != '\0') return value;
    }
    return "";
}

bool ends_with(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() &&
           text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// Every lookup of the 4-index tables and a sample of the internal-loop tables agree.
void expect_same_tables(const rs::Datatable& a, const rs::Datatable& b) {
    ASSERT_EQ(a.size(), b.size());
    const auto n = static_cast<rs::Code>(a.size());
    for (rs::Code i = 0; i < n; ++i) {
        for (rs::Code j = 0; j < n; ++j) {
            for (rs::Code k = 0; k < n; ++k) {
                for (rs::Code l = 0; l < n; ++l) {
                    ASSERT_EQ(a.stack(i, j, k, l), b.stack(i, j, k, l));
                    ASSERT_EQ(a.helix_ends(i, j, k, l), b.helix_ends(i, j, k, l));
                    ASSERT_EQ(a.tstkh(i, j, k, l), b.tstkh(i, j, k, l));
                    ASSERT_EQ(a.tstki(i, j, k, l), b.tstki(i, j, k, l));
                    ASSERT_EQ(a.tstkm(i, j, k, l), b.tstkm(i, j, k, l));
                    ASSERT_EQ(a.tstack(i, j, k, l), b.tstack(i, j, k, l));
                    ASSERT_EQ(a.coax(i, j, k, l), b.coax(i, j, k, l));
                    ASSERT_EQ(a.coaxstack(i, j, k, l), b.coaxstack(i, j, k, l));
                    ASSERT_EQ(a.tstackcoax(i, j, k, l), b.tstackcoax(i, j, k, l));
                    ASSERT_EQ(a.dangle(i, j, k, 1), b.dangle(i, j, k, 1));
                    ASSERT_EQ(a.dangle(i, j, k, 2), b.dangle(i, j, k, 2));
                }
            }
        }
    }
    std::mt19937 rng(42);
    std::uniform_int_distribution<int> pick(0, static_cast<int>(n) - 1);
    for (int trial = 0; trial < 20000; ++trial) {
        rs::Code c[8];
        for (rs::Code& x : c) x = static_cast<rs::Code>(pick(rng));
        ASSERT_EQ(a.iloop11(c[0], c[1], c[2], c[3], c[4], c[5]),
                  b.iloop11(c[0], c[1], c[2], c[3], c[4], c[5]));
        ASSERT_EQ(a.iloop21(c[0], c[1], c[2], c[3], c[4], c[5], c[6]),
                  b.iloop21(c[0], c[1], c[2], c[3], c[4], c[5], c[6]));
        ASSERT_EQ(a.iloop22(c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]),
                  b.iloop22(c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]));
        ASSERT_EQ(a.tloop(c), b.tloop(c));
    }
    EXPECT_EQ(a.misc().efn2a, b.misc().efn2a);
    EXPECT_EQ(a.hairpin(), b.hairpin());
}

}  // namespace

TEST(RNAstructureParams, LoadsTurner2004WithSeveralModifications) {
    const knotergy::rnastructure_param p =
        knotergy::RNAstructureParams::load("turner2004", {"m6A", "Ψ", "m1Ψ"});

    const knotergy::ParamSourceInfo& info = p.get_source_info();
    EXPECT_EQ(info.label, "RNAstructure");
    EXPECT_EQ(info.status, knotergy::ParamStatus::LoadedUserFile);
    EXPECT_EQ(info.resolved_name, "turner2004 + m6A + Ψ + m1Ψ");
    EXPECT_TRUE(ends_with(info.resolved_path, "params/RNAstructure/turner2004"));
    EXPECT_FALSE(info.messages.empty());  // ambiguous 'M', PSU's short coaxstack rows

    const rs::Datatable& data = p.data;
    EXPECT_EQ(data.size(), 9u);
    EXPECT_EQ(data.misc().efn2a, 93);       // efn2 multibranch offset, 0.1 kcal/mol
    EXPECT_EQ(data.stack(1, 4, 1, 4), -9);  // 5'AA/3'UU
    EXPECT_DOUBLE_EQ(data.temperature(), 310.15);
    EXPECT_EQ(data.origin(data.code("P")), "Ψ");

    EXPECT_DOUBLE_EQ(p.pseudoknot.values[0], 6.5);  // pseudconst.dat: a
    EXPECT_DOUBLE_EQ(p.pseudoknot.values[15], 30000.0);
}

TEST(RNAstructureParams, DefaultsToTurner2004) {
    const knotergy::rnastructure_param p = knotergy::RNAstructureParams::load();
    EXPECT_EQ(p.get_source_info().status, knotergy::ParamStatus::Defaulted);
    EXPECT_EQ(p.get_source_info().resolved_name, "turner2004");
    EXPECT_EQ(p.data.size(), 6u);

    rs::LoadOptions options;
    options.temperature = 300.0;
    const knotergy::rnastructure_param q = knotergy::RNAstructureParams::load(options);
    EXPECT_EQ(q.get_source_info().resolved_name, "turner2004");
    EXPECT_DOUBLE_EQ(q.data.temperature(), 300.0);
}

// efn2's alphabet names ("rna", "PSU", "1") and other spellings find the same folders.
TEST(RNAstructureParams, RNAstructureNamesFindTheSameSets) {
    const knotergy::rnastructure_param named =
        knotergy::RNAstructureParams::load("turner2004", {"m6A", "Ψ", "m1Ψ"});
    const knotergy::rnastructure_param efn2 =
        knotergy::RNAstructureParams::load("rna", {"m6a", "PSU", "1"});
    EXPECT_EQ(efn2.get_source_info().resolved_name, "turner2004 + m6A + Ψ + m1Ψ");
    expect_same_tables(named.data, efn2.data);

    const std::string folder = knotergy::resolve_rnastructure_set("m1Ψ");
    EXPECT_TRUE(ends_with(folder, "params/RNAstructure/m1Ψ"));
    EXPECT_EQ(knotergy::resolve_rnastructure_set(folder), folder);
    EXPECT_EQ(knotergy::resolve_rnastructure_set("Full1"),
              knotergy::resolve_rnastructure_set("m1Ψ_full"));
}

TEST(RNAstructureParams, EveryShippedSetLoads) {
    const std::vector<rs::ParameterSetLocation> sets =
        rs::list_set_folders(knotergy::default_rnastructure_param_path());
    ASSERT_EQ(sets.size(), 8u);
    for (const rs::ParameterSetLocation& set : sets) {
        const knotergy::rnastructure_param p = knotergy::RNAstructureParams::load(set.name);
        EXPECT_EQ(p.get_source_info().resolved_name, set.name);
        EXPECT_GE(p.data.size(), 6u) << set.name;
    }
    // Every U is m1Ψ in m1Ψ_full.
    const knotergy::rnastructure_param full = knotergy::RNAstructureParams::load("m1Ψ_full");
    EXPECT_EQ(full.data.encode("GG1UC"), full.data.encode("GGUUC"));
}

// Knotergy's modified-base symbols ('6' = m6A, 'P' = pseudouridine, '1' = m1Ψ) are also
// symbols of RNAstructure's m6A, PSU and 1 sets, so the grapheme views Knotergy already
// computes can be encoded directly.
TEST(RNAstructureParams, EncodesKnotergyGraphemeViews) {
    const knotergy::rnastructure_param p =
        knotergy::RNAstructureParams::load("turner2004", {"m6A", "Ψ", "m1Ψ"});

    std::vector<knotergy::modified_base_param> json_params =
        knotergy::ModParams::load_modified_energy_parameters(mod_folder);
    const knotergy::all_mod_params mp{json_params};
    const knotergy::ProcessedRNAEntry entry =
        knotergy::RNAProcessor::process_rna(knotergy::RNAEntry("GG6AAPUCC1", "((......))"), mp);
    ASSERT_TRUE(entry.has_modified_bases());

    const std::vector<rs::Code> codes = p.data.encode(entry.get_modified_sequence());
    EXPECT_EQ(codes, p.data.encode(entry.get_raw_sequence()));
    EXPECT_EQ(p.data.origin(codes[2]), "m6A");
    EXPECT_EQ(p.data.origin(codes[5]), "Ψ");
    EXPECT_EQ(p.data.origin(codes[9]), "m1Ψ");
    // The unmodified fallback agrees with Knotergy's own unmodified sequence.
    for (std::size_t i = 0; i < codes.size(); ++i) {
        EXPECT_EQ(p.data.symbol(p.data.fallback(codes[i])),
                  std::string(1, entry.get_sequence()[i]));
    }
}

TEST(RNAstructureParams, MultiByteSymbols) {
    rs::LoadOptions options;
    options.modifications = {{knotergy::resolve_rnastructure_set("Ψ"), "", {}, {{"Ψ", "P"}}},
                             {knotergy::resolve_rnastructure_set("m6A")}};
    const knotergy::rnastructure_param p = knotergy::RNAstructureParams::load(options);

    const std::string raw = "GGΨA6ΨCC";
    const std::vector<std::string_view> views =
        knotergy::ProcessedRNAEntry::parse_modified_sequence(raw);
    EXPECT_EQ(p.data.encode(views), p.data.encode("GGPA6PCC"));
    EXPECT_EQ(p.data.encode(raw), p.data.encode("GGPA6PCC"));
}

TEST(RNAstructureParams, TemperatureInCelsius) {
    const knotergy::rnastructure_param p =
        knotergy::RNAstructureParams::load("turner2004", {}, 25.0);
    EXPECT_NEAR(p.data.temperature(), 298.15, 1e-9);

    // m6A has no enthalpy tables: efn2 refuses other temperatures; Knotergy may opt in.
    EXPECT_THROW((void) knotergy::RNAstructureParams::load("turner2004", {"m6A"}, 25.0),
                 knotergy::DetailedException);
    const knotergy::rnastructure_param kept =
        knotergy::RNAstructureParams::load("turner2004", {"m6A"}, 25.0, true);
    EXPECT_DOUBLE_EQ(kept.data.modification(0).temperature(), 310.15);

    // RNAstructure 6.6 has no dna.helix_ends.dh, so DNA cannot be rescaled either.
    EXPECT_THROW((void) knotergy::RNAstructureParams::load("dna", {}, 25.0),
                 knotergy::DetailedException);
}

TEST(RNAstructureParams, ErrorsUseKnotergyExceptions) {
    EXPECT_THROW((void) knotergy::RNAstructureParams::load("/nonexistent/data_tables"),
                 knotergy::DetailedException);
    EXPECT_THROW((void) knotergy::RNAstructureParams::load("turner2004", {"nope"}),
                 knotergy::DetailedException);
    try {
        (void) knotergy::RNAstructureParams::load("turner2005");
        FAIL();
    } catch (const knotergy::DetailedException& e) {
        EXPECT_NE(std::string(e.what()).find("'turner2004' (rna)"), std::string::npos) << e.what();
    }
}

// The shipped folders hold RNAstructure 6.6's data_tables files, so both layouts give the
// same tables (only the set names differ).
TEST(RNAstructureParams, MatchesRNAstructureDataTables) {
    if (rnastructure_data().empty()) {
        GTEST_SKIP() << "Set DATAPATH to RNAstructure's data_tables directory";
    }
    const knotergy::rnastructure_param flat =
        knotergy::RNAstructureParams::load(rnastructure_data(), {"m6A", "PSU", "1"});
    EXPECT_EQ(flat.get_source_info().resolved_name, "rna + m6A + PSU + 1");
    const knotergy::rnastructure_param shipped =
        knotergy::RNAstructureParams::load("turner2004", {"m6A", "Ψ", "m1Ψ"});
    expect_same_tables(flat.data, shipped.data);
}
