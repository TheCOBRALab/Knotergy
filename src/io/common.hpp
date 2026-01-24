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
class DetailedException : public std::runtime_error {
   public:
    DetailedException(const std::string& message, const char* file, int line, const char* func)
        : std::runtime_error(format_message(message, file, line, func)) {}

   private:
    static std::string format_message(const std::string& message, const char* file, int line,
                                      const char* func) {
        std::ostringstream oss;
        oss << "Error: " << message << "\n"
            << "Function: " << func << "\n"
            << "File: " << file << ":" << line;
        return oss.str();
    }
};

#define THROW_ERROR(msg) throw DetailedException((msg), __FILE__, __LINE__, __func__)

// Max size of size_t
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
    s.erase(0, s.find_first_not_of(" \t\n\r\f\v\""));
    s.erase(s.find_last_not_of(" \t\n\r\f\v\"") + 1);
}

}  // namespace knotergy