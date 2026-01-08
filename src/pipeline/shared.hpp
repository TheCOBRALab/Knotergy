#pragma once

#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <dirent.h>
#include <cerrno>
#include <cstring>
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

// filesystem::exists not supported in older macOS Conda packages
inline bool file_exists(const std::string& name) {
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}

// Check if path is a file
inline bool is_file(const std::string& name) {
    struct stat buffer;
    if (stat(name.c_str(), &buffer) != 0) {
        return false;
    }
    return S_ISREG(buffer.st_mode);
}

// Check if path is a directory
inline bool is_directory(const std::string& name) {
    struct stat buffer;
    if (stat(name.c_str(), &buffer) != 0) {
        return false;
    }
    return S_ISDIR(buffer.st_mode);
}

// List files in a directory
inline std::vector<std::string> list_files_in_dir(const std::string& dir) {
    std::vector<std::string> out;

    DIR* d = ::opendir(dir.c_str());
    if (!d) {
        // up to you: throw, or return empty
        throw std::runtime_error("opendir failed: " + dir + " (" + std::strerror(errno) + ")");
    }

    while (dirent* e = ::readdir(d)) {
        // skip . and ..
        if (std::strcmp(e->d_name, ".") == 0 || std::strcmp(e->d_name, "..") == 0) continue;

        std::string full = dir;
        if (!full.empty() && full.back() != '/') full += '/';
        full += e->d_name;

        if (is_file(full)) out.push_back(full);
    }

    ::closedir(d);
    return out;
}

}  // namespace knotergy    