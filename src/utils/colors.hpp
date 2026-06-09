#pragma once

#include <string>

namespace knotergy {

#define ANSI_COLOR_BRIGHT "\x1b[1m"     // Bold / bright text
#define ANSI_COLOR_UNDERLINE "\x1b[4m"  // Underlined text

#define ANSI_COLOR_RED "\x1b[31m"      // Red:     #800000 or #FF0000
#define ANSI_COLOR_GREEN "\x1b[32m"    // Green:   #008000 or #00FF00
#define ANSI_COLOR_YELLOW "\x1b[33m"   // Yellow:  #808000 or #FFFF00
#define ANSI_COLOR_BLUE "\x1b[34m"     // Blue:    #000080 or #0000FF
#define ANSI_COLOR_MAGENTA "\x1b[35m"  // Magenta: #800080 or #FF00FF
#define ANSI_COLOR_CYAN "\x1b[36m"     // Cyan:    #008080 or #00FFFF

#define ANSI_COLOR_RED_B "\x1b[1;31m"      // Bright red:     #FF0000
#define ANSI_COLOR_GREEN_B "\x1b[1;32m"    // Bright green:   #00FF00
#define ANSI_COLOR_YELLOW_B "\x1b[1;33m"   // Bright yellow:  #FFFF00
#define ANSI_COLOR_BLUE_B "\x1b[1;34m"     // Bright blue:    #0000FF
#define ANSI_COLOR_MAGENTA_B "\x1b[1;35m"  // Bright magenta: #FF00FF
#define ANSI_COLOR_CYAN_B "\x1b[1;36m"     // Bright cyan:    #00FFFF

#define ANSI_COLOR_RESET "\x1b[0m"  // Reset / default terminal color

#define WARNING (ANSI_COLOR_YELLOW + std::string("[WARNING]") + ANSI_COLOR_RESET)
#define ERROR (ANSI_COLOR_RED_B + std::string("[ERROR]  ") + ANSI_COLOR_RESET)

}  // namespace knotergy