#pragma once

#include "io/parameters/PseudoknotParams.hpp"
#include "utils/common.hpp"

#include <string>
#include <vector>

namespace knotergy {

// What the caller should do once the command line has been read.
enum class ParseStatus {
    Ok,           // Carry on; the CliArgs instance is populated.
    ExitSuccess,  // A terminal flag such as --help or --version was handled.
    ExitFailure,  // The command line was malformed; a message was already printed.
};

// Every setting Knotergy takes from the command line.
struct CliArgs {
    // Input.
    std::string sequence;
    std::string structure;
    std::string input_file;

    // Parameter files.
    std::string vienna_param_file;
    std::string pseudo_param_file = default_pk_param_path();
    std::vector<std::string> mod_param_paths;

    // Energy model.
    double temperature = 37.0;  // Celsius.
    int dangle = 2;
    int round_value = 0;  // 0 means no rounding; see RoundMethod.

    // Behaviour.
    bool verbose = false;
    bool efn2_correction = false;
    bool show_input = false;
    bool disable_cache = false;

    // Reads argv into `out`. Problems are reported on stderr, so the caller
    // only has to act on the returned status.
    static ParseStatus parse(int argc, char** argv, CliArgs& out);

    // Prompts on stdin for whichever of sequence/structure is still missing,
    // unless an input file was supplied.
    void prompt_for_missing_input();

    // Checks that referenced files exist and numeric options are in range.
    // Reports on stderr and returns false when something is wrong.
    bool validate() const;

    RoundMethod round_method() const { return static_cast<RoundMethod>(round_value); }
};

// Prints the usage block to stdout.
void print_help();

// True when -v/--verbose appears anywhere in argv. Checked before parsing so
// that a failure during parsing can still honour the flag.
bool has_verbose_flag(int argc, char** argv);

}  // namespace knotergy
