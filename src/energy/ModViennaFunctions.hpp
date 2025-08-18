#pragma once

namespace knotergy {
    class ModViennaFunctions {
       public:
        static bool is_modified(char c) {
            return !(c == 'A' || c == 'U' || c == 'G' || c == 'C' || c == 'T');
        }
       private:

    };
}