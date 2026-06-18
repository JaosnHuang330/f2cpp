#include "mod_vector_algebra.h"

#include <cmath>

namespace mod_vector_algebra {

// function vector_sum(v1, v2) result(v3)
vector_t vector_sum(const vector_t& v1, const vector_t& v2) {
    vector_t v3;
    v3.x = v1.x + v2.x;
    v3.y = v1.y + v2.y;
    v3.z = v1.z + v2.z;
    return v3;
}

// function vector_subs(v1, v2) result(v3)
vector_t vector_subs(const vector_t& v1, const vector_t& v2) {
    vector_t v3;
    v3.x = v1.x - v2.x;
    v3.y = v1.y - v2.y;
    v3.z = v1.z - v2.z;
    return v3;
}

// function vector_dot_product(v1, v2) result(d)
double vector_dot_product(const vector_t& v1, const vector_t& v2) {
    double d;
    d = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    return d;
}

// function vector_cross_product(v1, v2) result(v3)
vector_t vector_cross_product(const vector_t& v1, const vector_t& v2) {
    vector_t v3;
    v3.x = v1.y * v2.z - v1.z * v2.y;
    v3.y = v1.z * v2.x - v1.x * v2.z;
    v3.z = v1.x * v2.y - v1.y * v2.x;
    return v3;
}

// function vector_norm(v) result(l)
double vector_norm(const vector_t& v) {
    double l;
    l = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    return l;
}

// function vector_times_scalar(v1, v2) result(v3)
vector_t vector_times_scalar(const vector_t& v1, double v2) {
    vector_t v3;
    v3.x = v1.x * v2;
    v3.y = v1.y * v2;
    v3.z = v1.z * v2;
    return v3;
}

// Operator overloads delegating to named functions
vector_t operator+(const vector_t& v1, const vector_t& v2) {
    return vector_sum(v1, v2);
}

vector_t operator-(const vector_t& v1, const vector_t& v2) {
    return vector_subs(v1, v2);
}

double operator*(const vector_t& v1, const vector_t& v2) {
    return vector_dot_product(v1, v2);
}

} // namespace mod_vector_algebra
