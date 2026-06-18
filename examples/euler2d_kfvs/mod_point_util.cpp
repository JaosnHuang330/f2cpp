#include "mod_point_util.h"

#include <cmath>

namespace mod_point_util {

// function equal(p1, p2) result(res)
// Returns true if p1 and p2 are within 1e-9 in each coordinate
bool equal(const point& p1, const point& p2) {
    bool res;
    if (std::abs(p1.x - p2.x) <= 1.0e-9 &&
        std::abs(p1.y - p2.y) <= 1.0e-9 &&
        std::abs(p1.z - p2.z) <= 1.0e-9) {
        res = true;
    } else {
        res = false;
    }
    return res;
}

// interface operator(.eq.) → operator==
bool operator==(const point& p1, const point& p2) {
    return equal(p1, p2);
}

// subroutine swap_point(p1, p2)
void swap_point(point& p1, point& p2) {
    point temp;
    temp = p1;
    p1   = p2;
    p2   = temp;
}

} // namespace mod_point_util
