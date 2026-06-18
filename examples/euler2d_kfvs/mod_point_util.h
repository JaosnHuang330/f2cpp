#pragma once

#include <cmath>
#include "mod_read_gmsh.h"

namespace mod_point_util {

using mod_read_gmsh::point;

// equal: compare two points within 1e-9 tolerance
bool equal(const point& p1, const point& p2);

// operator(.eq.) → operator==
bool operator==(const point& p1, const point& p2);

// swap_point: swap two points
void swap_point(point& p1, point& p2);

} // namespace mod_point_util
