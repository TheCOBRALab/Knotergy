#include "energy/ComputeEnergy.hpp"
#include "io/input/CliArgs.hpp"
#include "io/input/RNAInputManager.hpp"
#include "io/output/OutputManager.hpp"
#include "io/output/colors.hpp"
#include "io/parameters/ModParams.hpp"
#include "io/parameters/PseudoknotParams.hpp"
#include "io/parameters/ViennaParams.hpp"
#include "loop_tree/LoopFactory.hpp"
#include "preprocessing/ProcessedRNAEntry.hpp"
#include "preprocessing/RNAEntry.hpp"
#include "utils/common.hpp"

#include <cstdio>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

namespace {

std::vector<knotergy::modified_base_param> load_modified_params(
    const std::vector<std::string>& mod_param_paths) {
    std::vector<knotergy::modified_base_param> modified_params;

    for (const std::string& mod_path : mod_param_paths) {
        std::vector<knotergy::modified_base_param> additional_mp =
            knotergy::ModParams::load_modified_energy_parameters(mod_path);

        modified_params.insert(modified_params.end(), additional_mp.begin(), additional_mp.end());
    }

    return modified_params;
}

int run_knotergy(const knotergy::CliArgs& args) {
    const bool use_color = knotergy::should_use_color();

    // ------------------------- Load Parameters -----------------------
    knotergy::vrna_md_param vp = knotergy::ViennaParams::load_energy_parameters(
        args.vienna_param_file, args.dangle, args.sequence, args.temperature, args.disable_cache);

    const knotergy::pk_param pkp =
        knotergy::PseudoknotParams::load_pk_param(args.pseudo_param_file, args.round_method());

    const knotergy::all_mod_params mp = load_modified_params(args.mod_param_paths);

    // ------------------------- Reading Inputs From File -----------------------------
    const std::vector<knotergy::RNAEntry> inputs =
        knotergy::RNAInputManager::get_all_inputs(args.input_file, args.sequence, args.structure);

    // ------------------------- Print Parameter Report -----------------------
    knotergy::OutputManager::print_parameter_report(vp, pkp, mp, args.efn2_correction);

    // ------------------------- Main Processing Loop ----------------------------
    for (const knotergy::RNAEntry& rna : inputs) {
        std::cout << "\n--------- Name: " << rna.name << " ---------" << '\n';
        if (args.show_input) {
            std::cout << "Sequence : " << rna.sequence << '\n';
            std::cout << "Structure: " << rna.structure << '\n';
        }

        // Preprocess the RNA entry to compute pair_table, closed regions, etc.
        const knotergy::ProcessedRNAEntry& processed_rna =
            knotergy::RNAProcessor::process_rna(rna, mp);

        // Build loop tree.
        knotergy::LoopFactory factory(processed_rna);

        // Compute the energy.
        knotergy::ComputeEnergy energy_calculator(factory.get_root_node(), processed_rna, vp, pkp,
                                                  mp, args.efn2_correction, args.verbose);

        // Output results.
        if (energy_calculator.getInfiniteEnergyFlag()) {
            printf("\nENERGY: Infinite (%.4f kcal/mol)\n", energy_calculator.getEnergy());
        } else if (use_color) {
            printf("\nENERGY:%s %.4f kcal/mol%s\n", ANSI_COLOR_GREEN, energy_calculator.getEnergy(),
                   ANSI_COLOR_RESET);
        } else {
            printf("\nENERGY: %.4f kcal/mol\n", energy_calculator.getEnergy());
        }
    }

    return EXIT_SUCCESS;
}

}  // namespace

int main(int argc, char** argv) {
    const bool verbose = knotergy::has_verbose_flag(argc, argv);

    try {
        knotergy::CliArgs args;

        switch (knotergy::CliArgs::parse(argc, argv, args)) {
            case knotergy::ParseStatus::ExitSuccess: return EXIT_SUCCESS;
            case knotergy::ParseStatus::ExitFailure: return EXIT_FAILURE;
            case knotergy::ParseStatus::Ok:          break;
        }

        args.prompt_for_missing_input();

        if (!args.validate()) return EXIT_FAILURE;

        return run_knotergy(args);

    } catch (const knotergy::DetailedException& e) {
        if (verbose) {
            std::cerr << e.detailed_message() << '\n';
        } else {
            std::cerr << ERROR << " " << e.what() << '\n';
        }

        return EXIT_FAILURE;

    } catch (const std::exception& e) {
        std::cerr << ERROR << " " << e.what() << '\n';
        return EXIT_FAILURE;

    } catch (...) {
        std::cerr << ERROR << " Unknown fatal error\n";
        return EXIT_FAILURE;
    }
}
