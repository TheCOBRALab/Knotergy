#pragma once
#include <cstddef>
#include <ostream>
#include <vector>

#include "../pipeline/shared.hpp"

namespace knotergy {

struct Pair {
    size_t i;
    size_t j;

    Pair() = default;
    Pair(size_t i, size_t j) : i(i), j(j) {}
    bool is_stack(Pair child) const {
        return i + 1 == child.i && j - 1 == child.j;
    }
};

// Visual representation of the band:
// (((((((((((((((............)))))))))))))))
// ^             ^            ^             ^
// left_border   left_inner   right_inner   right_border
class Band {
   public:
    Band(size_t lb, size_t li, size_t ri, size_t rb, const std::vector<size_t>& pairings)
        : left_border_{lb}, left_inner_{li}, right_inner_{ri}, right_border_{rb} {
        for (size_t idx = left_border_; idx <= left_inner_; ++idx) {
            size_t paired = pairings[idx];
            if (paired >= right_inner_ && paired <= right_border_) {
                base_pairs_.emplace_back(idx, paired);
            }
        }
    }

    bool contains(size_t idx) const {
        return (idx >= left_border_ && idx <= left_inner_) ||
               (idx >= right_inner_ && idx <= right_border_);
    }

    // === Read-only accessors ===
    size_t left_border() const { return left_border_; }
    size_t left_inner() const { return left_inner_; }
    size_t right_inner() const { return right_inner_; }
    size_t right_border() const { return right_border_; }

    const std::vector<Pair>& base_pairs() const { return base_pairs_; }

   private:
    size_t left_border_;   // i
    size_t left_inner_;    // i`
    size_t right_inner_;   // j`
    size_t right_border_;  // j

    std::vector<Pair> base_pairs_;
};

// === Operator overload for printing ===
inline std::ostream& operator<<(std::ostream& os, const Band& band) {
    os << "Band(" << band.left_border() << ", " << band.left_inner() << ", " << band.right_inner()
       << ", " << band.right_border() << ")";
    return os;
}

}  // namespace knotergy
