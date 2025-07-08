#pragma once

#include <sstream>
#include <stdexcept>
#include <string>

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

}  // namespace knotergy