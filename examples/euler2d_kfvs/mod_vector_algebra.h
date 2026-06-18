#pragma once

#include <cmath>

namespace mod_vector_algebra {

// Derived type: vector_t
struct vector_t {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

// Operator overloads for standard Fortran operators
// operator(+) → vector_sum
vector_t operator+(const vector_t& v1, const vector_t& v2);

// operator(-) → vector_subs
vector_t operator-(const vector_t& v1, const vector_t& v2);

// operator(*) → vector_dot_product
double operator*(const vector_t& v1, const vector_t& v2);

// Named functions for custom Fortran operators (.x., .abs., .times.)
// operator(.x.) → vector_cross_product
vector_t vector_cross_product(const vector_t& v1, const vector_t& v2);

// operator(.abs.) → vector_norm
double vector_norm(const vector_t& v);

// operator(.times.) → vector_times_scalar
vector_t vector_times_scalar(const vector_t& v1, double v2);

// Also expose the named versions of the standard operators for direct calls
vector_t vector_sum(const vector_t& v1, const vector_t& v2);
vector_t vector_subs(const vector_t& v1, const vector_t& v2);
double   vector_dot_product(const vector_t& v1, const vector_t& v2);

} // namespace mod_vector_algebra
