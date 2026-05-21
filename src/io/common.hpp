#pragma once

#include <dirent.h>
#include <sys/stat.h>

#include <cerrno>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace knotergy {
/**
 * @brief Custom exception class with detailed context information.
 *
 * Extends std::runtime_error to include file, line, and function information
 * for easier debugging.
 */
class DetailedException : public std::runtime_error {
   public:
    /**
     * @brief Construct a detailed exception.
     *
     * @param message Error message.
     * @param file Source file where error occurred.
     * @param line Line number where error occurred.
     * @param func Function name where error occurred.
     */
    DetailedException(const std::string& message, const char* file, int line, const char* func)
        : std::runtime_error(format_message(message, file, line, func)) {}

   private:
    /**
     * @brief Format the error message with context.
     *
     * @param message Error message.
     * @param file Source file.
     * @param line Line number.
     * @param func Function name.
     * @return Formatted error message string.
     */
    [[nodiscard]] static std::string format_message(const std::string& message, const char* file,
                                                    int line, const char* func) {
        std::ostringstream oss;

        oss << "\n"
            << "══════════════════════════════════════\n"
            << "❌ ERROR\n"
            << "══════════════════════════════════════\n"
            << message << "\n\n"
            << "📍 Location:\n"
            << "  Function : " << func << '\n'
            << "  File     : " << file << '\n'
            << "  Line     : " << line << '\n'
            << "══════════════════════════════════════";

        return oss.str();
    }
};

/// Macro to throw a DetailedException with current file, line, and function.
#define THROW_ERROR(msg) throw DetailedException((msg), __FILE__, __LINE__, __func__)

/// Maximum value of size_t, used as a sentinel for "no index" or "invalid index".
constexpr size_t NULL_INDEX = static_cast<size_t>(-1);

/**
 * @brief Trims leading and trailing whitespace from a string.
 *
 * This function modifies the input string in-place to remove any leading
 * and trailing whitespace characters, including spaces, tabs, newlines,
 * carriage returns, form feeds, and vertical tabs.
 *
 * @param s The string to be trimmed.
 */
inline void trim(std::string& s) {
    const char* ws = " \t\n\r\f\v\"";

    // Trim left
    size_t start = s.find_first_not_of(ws);
    if (start == std::string::npos) {
        s.clear();  // string is all whitespace
        return;
    }

    // Trim right
    size_t end = s.find_last_not_of(ws);

    s = s.substr(start, end - start + 1);
}

}  // namespace knotergy