#pragma once
#include <cstddef>

#include "../helpers/common.hpp"
namespace knotergy {
struct Band {
    size_t left_border = NULL_INDEX;  //  i
    size_t left_inner;                //  i′
    size_t right_inner;               //  j′
    size_t right_border;              //  j
    bool contains(size_t idx) const {
        return (idx >= left_border && idx <= left_inner) ||
               (idx >= right_inner && idx <= right_border);
    }
};

inline std::ostream& operator<<(std::ostream& os, const Band& band) {
    os << "Band(" << band.left_border << ", " << band.left_inner << ", "
       << band.right_inner << ", " << band.right_border << ")";
    return os;
}

}  // namespace knotergy
