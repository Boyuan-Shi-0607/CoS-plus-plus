#ifndef TYPES_H
#define TYPES_H

#include <map>
#include <set>
#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <iomanip>

using MP_Real    = boost::multiprecision::cpp_bin_float_100;
using MP_Complex = boost::multiprecision::cpp_complex_100;

// Get PI with full precision
inline MP_Real get_pi() {
    static MP_Real pi = boost::multiprecision::acos(MP_Real(-1));
    return pi;
}

inline MP_Real get_two_pi() {
    static MP_Real two_pi = MP_Real(2) * get_pi();
    return two_pi;
}

// DoubleDouble implementation for fast extended precision
struct DoubleDouble {
    double hi, lo;

    constexpr DoubleDouble() noexcept : hi(0), lo(0) {}

    // Careful conversion from int to avoid rounding errors
    constexpr explicit DoubleDouble(int x) noexcept : hi(static_cast<double>(x)), lo(0) {
        // For integers that fit exactly in double, lo is exactly 0
    }

    // Constructor for unsigned long long
    constexpr explicit DoubleDouble(unsigned long long x) noexcept : hi(static_cast<double>(x)), lo(0) {
        // For integers that fit exactly in double, lo is exactly 0
        // Note: for very large unsigned long long values, there may be precision loss
    }

    // Constructor for unsigned long long
    constexpr explicit DoubleDouble(long long x) noexcept : hi(static_cast<double>(x)), lo(0) {
        // For integers that fit exactly in double, lo is exactly 0
        // Note: for very large unsigned long long values, there may be precision loss
    }

    // Constructor for unsigned long long
    constexpr explicit DoubleDouble(unsigned long x) noexcept : hi(static_cast<double>(x)), lo(0) {
        // For integers that fit exactly in double, lo is exactly 0
        // Note: for very large unsigned long long values, there may be precision loss
    }

    constexpr explicit DoubleDouble(double x) noexcept : hi(x), lo(0) {}

    explicit DoubleDouble(const MP_Real& x) noexcept {
        hi = static_cast<double>(x);
        MP_Real residual = x - MP_Real(hi);
        lo = static_cast<double>(residual);
    }

    constexpr DoubleDouble(double h, double l) noexcept : hi(h), lo(l) {}

    // Arithmetic helpers
    static inline constexpr double quick_two_sum(double a, double b, double &err) noexcept {
        double s = a + b;
        err = b - (s - a);
        return s;
    }

    static inline constexpr double two_sum(double a, double b, double &err) noexcept {
        double s = a + b;
        double bb = s - a;
        err = (a - (s - bb)) + (b - bb);
        return s;
    }

    static inline constexpr double two_prod(double a, double b, double &err) noexcept {
        double p = a * b;
        err = std::fma(a, b, -p);
        return p;
    }

    // Renormalization
    constexpr void renormalize() noexcept {
        double s = hi + lo;
        lo = lo - (s - hi);
        hi = s;
    }

    // Arithmetic operations
    friend constexpr DoubleDouble operator+(const DoubleDouble &x, const DoubleDouble &y) noexcept {
        double s, e, t, f;
        s = two_sum(x.hi, y.hi, e);
        t = x.lo + y.lo;
        f = e + t;
        DoubleDouble r;
        r.hi = quick_two_sum(s, f, r.lo);
        return r;
    }

    friend constexpr DoubleDouble operator-(const DoubleDouble &x, const DoubleDouble &y) noexcept {
        double s, e, t, f;
        s = two_sum(x.hi, -y.hi, e);
        t = x.lo - y.lo;
        f = e + t;
        DoubleDouble r;
        r.hi = quick_two_sum(s, f, r.lo);
        return r;
    }

    friend constexpr DoubleDouble operator*(const DoubleDouble &x, const DoubleDouble &y) noexcept {
        double p, e;
        p = two_prod(x.hi, y.hi, e);
        e += x.hi * y.lo + x.lo * y.hi + x.lo * y.lo;
        DoubleDouble r;
        r.hi = quick_two_sum(p, e, r.lo);
        return r;
    }

    friend DoubleDouble operator/(const DoubleDouble &x, const DoubleDouble &y) {
        double q1 = x.hi / y.hi;
        DoubleDouble r = x - DoubleDouble(q1) * y;
        double q2 = (r.hi + r.lo) / y.hi;
        DoubleDouble q = DoubleDouble(q1) + DoubleDouble(q2);
        q.renormalize();
        return q;
    }

    // Division by int - careful to avoid rounding errors
    friend DoubleDouble operator/(const DoubleDouble &x, int n) {
        // Convert int to DoubleDouble exactly (no rounding error for reasonable ints)
        return x / DoubleDouble(n);
    }

    // Unary operators
    constexpr DoubleDouble operator-() const noexcept {
        return DoubleDouble(-hi, -lo);
    }

    // Compound assignments
    DoubleDouble& operator+=(const DoubleDouble &o) { return *this = *this + o; }
    DoubleDouble& operator-=(const DoubleDouble &o) { return *this = *this - o; }
    DoubleDouble& operator*=(const DoubleDouble &o) { return *this = *this * o; }
    DoubleDouble& operator/=(const DoubleDouble &o) { return *this = *this / o; }

    // Conversions
    constexpr double to_double() const noexcept { return hi + lo; }
    constexpr operator double() const noexcept { return to_double(); }

    MP_Real to_mp_real() const {
        return MP_Real(hi) + MP_Real(lo);
    }

    // Comparison operators
    friend bool operator<(const DoubleDouble& a, const DoubleDouble& b) {
        return (a.hi < b.hi) || (a.hi == b.hi && a.lo < b.lo);
    }

    friend bool operator>(const DoubleDouble& a, const DoubleDouble& b) {
        return b < a;
    }

    friend bool operator<=(const DoubleDouble& a, const DoubleDouble& b) {
        return !(a > b);
    }

    friend bool operator>=(const DoubleDouble& a, const DoubleDouble& b) {
        return !(a < b);
    }

    friend bool operator==(const DoubleDouble& a, const DoubleDouble& b) {
        return a.hi == b.hi && a.lo == b.lo;
    }

    friend bool operator!=(const DoubleDouble& a, const DoubleDouble& b) {
        return !(a == b);
    }
};

// For DoubleDouble operations that need PI
inline DoubleDouble DD_PI() {
    static DoubleDouble pi_dd(get_pi());
    return pi_dd;
}

inline DoubleDouble DD_TWO_PI() {
    static DoubleDouble two_pi_dd(get_two_pi());
    return two_pi_dd;
}

// Robust floor function for DoubleDouble
inline int dd_floor_to_int(const DoubleDouble& x) {
    // 0. Reject NaNs / infinities so the behaviour is defined.
    if (!std::isfinite(x.hi))
        throw std::domain_error("dd_floor_to_int: non-finite input");

    // 1. Floor the high component (quick).
    double hi_floor = std::floor(x.hi);

    // 2. If high word was already an integer,
    //    check whether the low word pushes the value below it.
    if (hi_floor == x.hi && x.lo < 0.0)      // low word negative
        hi_floor -= 1.0;                      // step down once

    // 3. Optional safety: ensure result fits in 'int'.
    constexpr auto INT_MIN_D = static_cast<double>(std::numeric_limits<int>::min());
    constexpr auto INT_MAX_D = static_cast<double>(std::numeric_limits<int>::max());
    if (hi_floor < INT_MIN_D || hi_floor > INT_MAX_D)
        throw std::overflow_error("dd_floor_to_int: result out of int range");

    return static_cast<int>(hi_floor);
}

// Alias for consistency
inline int DD_FLOOR_TO_INT(const DoubleDouble& x) {
    return dd_floor_to_int(x);
}

// Floor function that returns DoubleDouble
inline DoubleDouble dd_floor(const DoubleDouble& x) {
    int floor_val = dd_floor_to_int(x);
    return DoubleDouble(floor_val);
}

// Absolute value for DoubleDouble
inline DoubleDouble dd_abs(const DoubleDouble& x) {
    return (x.hi < 0.0) ? -x : x;
}

// For sqrt, exp, cos, sin etc - use cpp_dec_float
inline DoubleDouble dd_sqrt(const DoubleDouble& x) {
    MP_Real x_mp = x.to_mp_real();
    MP_Real result_mp = boost::multiprecision::sqrt(x_mp);
    return DoubleDouble(result_mp);
}

inline DoubleDouble dd_exp(const DoubleDouble& x) {
    MP_Real x_mp = x.to_mp_real();
    MP_Real result_mp = boost::multiprecision::exp(x_mp);
    return DoubleDouble(result_mp);
}

inline DoubleDouble dd_cos(const DoubleDouble& x) {
    MP_Real x_mp = x.to_mp_real();
    MP_Real result_mp = boost::multiprecision::cos(x_mp);
    return DoubleDouble(result_mp);
}

inline DoubleDouble dd_sin(const DoubleDouble& x) {
    MP_Real x_mp = x.to_mp_real();
    MP_Real result_mp = boost::multiprecision::sin(x_mp);
    return DoubleDouble(result_mp);
}

inline DoubleDouble dd_pow(const DoubleDouble& x, const DoubleDouble& y) {
    MP_Real x_mp = x.to_mp_real();
    MP_Real y_mp = y.to_mp_real();
    MP_Real result_mp = boost::multiprecision::pow(x_mp, y_mp);
    return DoubleDouble(result_mp);
}

inline DoubleDouble dd_atan2(const DoubleDouble& y, const DoubleDouble& x) {
    MP_Real y_mp = y.to_mp_real();
    MP_Real x_mp = x.to_mp_real();
    MP_Real result_mp = boost::multiprecision::atan2(y_mp, x_mp);
    return DoubleDouble(result_mp);
}

inline DoubleDouble dd_log(const DoubleDouble& x) {
    MP_Real x_mp = x.to_mp_real();
    MP_Real result_mp = boost::multiprecision::log(x_mp);
    return DoubleDouble(result_mp);
}

// Define Real as double for compatibility
using Real = double;
using Complex = std::complex<Real>;

// MP versions for data preparation
using MP_Real1DVector = std::vector<MP_Real>;
using MP_Real2DVector = std::vector<MP_Real1DVector>;
using MP_Real3DVector = std::vector<MP_Real2DVector>;
using MP_Real4DVector = std::vector<MP_Real3DVector>;
using MP_Real5DVector = std::vector<MP_Real4DVector>;

using MP_Complex1DVector = std::vector<MP_Complex>;
using MP_Complex2DVector = std::vector<MP_Complex1DVector>;
using MP_Complex3DVector = std::vector<MP_Complex2DVector>;
using MP_Complex4DVector = std::vector<MP_Complex3DVector>;
using MP_Complex5DVector = std::vector<MP_Complex4DVector>;

// Real number vectors
using Real1DVector = std::vector<Real>;
using Real2DVector = std::vector<Real1DVector>;
using Real3DVector = std::vector<Real2DVector>;
using Real4DVector = std::vector<Real3DVector>;
using Real5DVector = std::vector<Real4DVector>;

// DoubleDouble vectors
using DD1DVector = std::vector<DoubleDouble>;
using DD2DVector = std::vector<DD1DVector>;

// Complex number vectors
using Complex1DVector = std::vector<Complex>;
using Complex2DVector = std::vector<Complex1DVector>;
using Complex3DVector = std::vector<Complex2DVector>;
using Complex4DVector = std::vector<Complex3DVector>;
using Complex5DVector = std::vector<Complex4DVector>;

// Integer vectors
using Int1DVector = std::vector<int>;
using Int2DVector = std::vector<std::vector<int>>;
using Int3DVector = std::vector<std::vector<std::vector<int>>>;
using Int4DVector = std::vector<std::vector<std::vector<std::vector<int>>>>;
using Int5DVector = std::vector<std::vector<std::vector<std::vector<std::vector<int>>>>>;

// Output stream operator for DoubleDouble
inline std::ostream& operator<<(std::ostream& os, const DoubleDouble& x) {
    return os << std::fixed << std::setprecision(32) << x.to_double();
}

#endif // TYPES_H