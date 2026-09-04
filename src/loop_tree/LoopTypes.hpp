#pragma once

namespace knotergy {
/**
 * @brief Enumeration of loop types in RNA secondary structures.
 */
enum class LoopType { Unknown, Stack, Hairpin, Internal, Multibranch, External, Pseudoknot };

/**
 * @brief Enumeration of pseudoknot nesting types.
 *
 * Within Band ((..(...).(.[[[.)..))]]]
 * This hairpin     ^   ^ is within a band
 * Nested      (((..(...)..[[[...)))]]]
 * This hairpin     ^   ^ is not within a band
 */
enum class PseudoNestedType { None, WithinBand, OutsideBandIntervals };

/**
 * @brief Get a human-readable name for a loop type.
 *
 * @param t The loop type.
 * @return String representation of the loop type.
 */
[[nodiscard]] static inline const char* loop_name(LoopType t) {
    switch (t) {
        case LoopType::Stack:       return "Stack        ";
        case LoopType::Hairpin:     return "Hairpin  loop";
        case LoopType::Internal:    return "Interior loop";
        case LoopType::Multibranch: return "Multi    loop";
        case LoopType::External:    return "External loop";
        case LoopType::Pseudoknot:  return "Pseudo   loop";
        case LoopType::Unknown:     return "Unknown  loop";
    }
    return "Unknown  loop";
}

[[nodiscard]] static inline const char* pk_nested_type(PseudoNestedType t) {
    switch (t) {
        case PseudoNestedType::None:                 return "None";
        case PseudoNestedType::WithinBand:           return "WithinBand";
        case PseudoNestedType::OutsideBandIntervals: return "OutsideBandIntervals";
    }
    return "Unknown pk nesting type";
};

}  // namespace knotergy