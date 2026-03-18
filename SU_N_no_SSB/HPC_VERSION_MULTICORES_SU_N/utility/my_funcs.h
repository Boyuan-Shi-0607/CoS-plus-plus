#ifndef MY_FUNCS_H
#define MY_FUNCS_H

#include "../utility/types.h"
#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>

// Standard precision functions (for Real = double)
inline Real MY_COS(const Real x) {
    return std::cos(x);
}

inline Real MY_SIN(const Real x) {
    return std::sin(x);
}

inline Real MY_POW(const Real x, const Real y) {
    return std::pow(x, y);
}

inline Real MY_SQRT(const Real x) {
    return std::sqrt(x);
}

inline Real MY_ATAN2(const Real y, const Real x) {
    return std::atan2(y, x);
}

inline Real MY_EXP(const Real x) {
    return std::exp(x);
}

inline Real MY_ROUND(const Real x) {
    return std::round(x);
}

inline Real MY_LOG2(const Real x) {
    return std::log2(x);
}

inline Real MY_LOG(const Real x) {
    return std::log(x);
}

inline Real MY_FLOOR(const Real x) {
    return std::floor(x);
}

inline Real MY_ABS(const Real x) {
    return std::abs(x);
}

inline Real MY_CEIL(const Real x) {
    return std::ceil(x);
}

inline Real MY_ACOS(const Real x) {
    return std::acos(x);
}

inline Real MY_ASIN(const Real x) {
    return std::asin(x);
}

inline Real MY_ATAN(const Real x) {
    return std::atan(x);
}

inline Real MY_TAN(const Real x) {
    return std::tan(x);
}

// Helper function to check if a value is NaN
inline bool is_nan(const Real& x) {
    return std::isnan(x);
}

// Helper function to check if a value is infinite
inline bool is_inf(const Real& x) {
    return std::isinf(x);
}

// Helper function to check if a value is finite
inline bool is_finite(const Real& x) {
    return std::isfinite(x);
}

#endif //MY_FUNCS_H