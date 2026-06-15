#pragma once

#include <sys/stat.h>

#include <dirent.h>

#include <cerrno>
#include <cstddef>
#include <cstring>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace knotergy {

/**
 * @brief Custom exception class with both short and detailed error messages.
 *
 * what() returns only the short user-facing message.
 * detailed_message() returns the full debug message with file, line, and function.
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
        : std::runtime_error(message),
          message_(message),
          file_(file),
          line_(line),
          func_(func),
          detailed_message_(format_message(message, file, line, func)) {}

    /**
     * @brief Return the short user-facing error message.
     */
    [[nodiscard]] const std::string& message() const noexcept { return message_; }

    /**
     * @brief Return the full detailed debug message.
     */
    [[nodiscard]] const std::string& detailed_message() const noexcept { return detailed_message_; }

    /**
     * @brief Return source file where the exception was thrown.
     */
    [[nodiscard]] const std::string& file() const noexcept { return file_; }

    /**
     * @brief Return line number where the exception was thrown.
     */
    [[nodiscard]] int line() const noexcept { return line_; }

    /**
     * @brief Return function where the exception was thrown.
     */
    [[nodiscard]] const std::string& function() const noexcept { return func_; }

   private:
    std::string message_;
    std::string file_;
    int line_;
    std::string func_;
    std::string detailed_message_;

    /**
     * @brief Format the detailed error message with context.
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
#define THROW_ERROR(msg) throw ::knotergy::DetailedException((msg), __FILE__, __LINE__, __func__)

/// Maximum value of size_t, used as a sentinel for "no index" or "invalid index".
constexpr std::size_t NULL_INDEX = static_cast<std::size_t>(-1);
constexpr int NULL_ENERGY = std::numeric_limits<int>::max();

/**
 * @brief Trims leading and trailing whitespace from a string.
 *
 * This function modifies the input string in-place to remove any leading
 * and trailing whitespace characters, including spaces, tabs, newlines,
 * carriage returns, form feeds, vertical tabs, and quotes.
 *
 * @param s The string to be trimmed.
 */
inline void trim(std::string& s) {
    const char* ws = " \t\n\r\f\v\"";

    // Trim left
    std::size_t start = s.find_first_not_of(ws);
    if (start == std::string::npos) {
        s.clear();
        return;
    }

    // Trim right
    std::size_t end = s.find_last_not_of(ws);

    s = s.substr(start, end - start + 1);
}

}  // namespace knotergy