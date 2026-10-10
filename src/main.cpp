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

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

namespace {

// ------------- Parameter loading helper functions ----------------
knotergy::vrna_md_param load_vienna_params(const knotergy::CliArgs& args) {
    return knotergy::ViennaParams::load_energy_parameters(
        args.vienna_param_file, args.dangle, args.temperature, args.salt, args.disable_cache);
}

knotergy::pk_param load_pseudoknot_params(const knotergy::CliArgs& args) {
    return knotergy::PseudoknotParams::load_pk_param(args.pseudo_param_file, args.round_method());
}

knotergy::all_mod_params load_modified_params(const knotergy::CliArgs& args) {
    return knotergy::ModParams::load_modified_params(args.mod_param_paths);
}

struct EnergyResult {
    double energy = 0.0;
    bool infinite_energy = false;
    std::exception_ptr error;
};

EnergyResult evaluate_entry(const knotergy::RNAEntry& rna, const knotergy::CliArgs& args,
                            knotergy::vrna_md_param& vp, const knotergy::pk_param& pkp,
                            const knotergy::all_mod_params& mp) {
    const knotergy::ProcessedRNAEntry processed_rna = knotergy::RNAProcessor::process_rna(rna, mp);
    knotergy::LoopFactory factory(processed_rna);
    knotergy::ComputeEnergy calculator(factory.get_root_node(), processed_rna, vp, pkp, mp,
                                      args.pk_dangles, args.efn2_correction, args.verbosity);
    return {calculator.getEnergy(), calculator.getInfiniteEnergyFlag(), {}};
}

void print_per_entry_header(const knotergy::RNAEntry& rna, const knotergy::CliArgs& args) {
    std::cout << "\n--------- Name: " << rna.name << " ---------" << '\n';
    if (args.show_input) {
        std::cout << "Sequence : " << rna.sequence << '\n';
        std::cout << "Structure: " << rna.structure << '\n';
    }
}

void print_results(const EnergyResult& result, const bool use_color) {
    if (result.infinite_energy) {
        if (use_color) {
            printf("\nENERGY:%s Infinite (%.4f kcal/mol)%s\n", ANSI_COLOR_RED,
                   result.energy, ANSI_COLOR_RESET);
        } else {
            printf("\nENERGY: Infinite (%.4f kcal/mol)\n", result.energy);
        }
    } else {
        if (use_color) {
            printf("\nENERGY:%s %.4f kcal/mol%s\n", ANSI_COLOR_GREEN, result.energy,
                   ANSI_COLOR_RESET);
        } else {
            printf("\nENERGY: %.4f kcal/mol\n", result.energy);
        }
    }
}

int run_knotergy(const knotergy::CliArgs& args) {
    const bool use_color = knotergy::should_use_color();

    // ------------------------- Load Parameters -----------------------
    // Not const cuz ViennaRNA accidentally didn't mark their functions as const
    knotergy::vrna_md_param vp = load_vienna_params(args);
    const knotergy::pk_param pkp = load_pseudoknot_params(args);
    const knotergy::all_mod_params mp = load_modified_params(args);

    // The header of knotergy which includes the breakdown of the parameters and settings
    knotergy::OutputManager::print_parameter_report(vp, pkp, mp, args.efn2_correction);

    // ------------------------- Reading Inputs From File -----------------------------
    const std::vector<knotergy::RNAEntry> inputs =
        knotergy::RNAInputManager::get_all_inputs(args.input_file, args.sequence, args.structure);

#ifdef KNOTERGY_ENABLE_PARALLEL
    // Detailed breakdowns still use stdout, so evaluate verbose input serially.
    const int jobs = args.verbosity == knotergy::VerbosityLevel::Quiet
                         ? std::max(1, std::min(static_cast<int>(inputs.size()), args.jobs))
                         : 1;
    if (jobs > 1) {
        std::cout << "Running in parallel with " << jobs << " threads\n";
        std::vector<EnergyResult> results(inputs.size());

#pragma omp parallel for num_threads(jobs) schedule(dynamic, 32)
        for (std::size_t i = 0; i < inputs.size(); ++i) {
            try {
                results[i] = evaluate_entry(inputs[i], args, vp, pkp, mp);
            } catch (...) {
                results[i].error = std::current_exception();
            }
        }

        // Workers compute; the main thread prints complete records in input order.
        for (std::size_t i = 0; i < inputs.size(); ++i) {
            print_per_entry_header(inputs[i], args);
            if (results[i].error) std::rethrow_exception(results[i].error);
            print_results(results[i], use_color);
        }
        return EXIT_SUCCESS;
    }
#endif
    // ------------------------- Main Processing Loop ----------------------------
    for (const knotergy::RNAEntry& rna : inputs) {
        print_per_entry_header(rna, args);

        print_results(evaluate_entry(rna, args, vp, pkp, mp), use_color);
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
