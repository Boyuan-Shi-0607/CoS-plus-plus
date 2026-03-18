#ifndef BIPARTITE_LATTICE_PRINCIPAL_MINOR_CDET_FAST_PRINCIPLE_MINORS_H
#define BIPARTITE_LATTICE_PRINCIPAL_MINOR_CDET_FAST_PRINCIPLE_MINORS_H

#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <random>
#include <array>
#include <type_traits>
#include "parameters.h"
#include "../utility/types.h"

// Check if T is double or DoubleDouble
template<typename T>
static constexpr bool IS_DOUBLE = std::is_same_v<T, double>;

template<typename T>
static constexpr bool IS_DOUBLEDOUBLE = std::is_same_v<T, DoubleDouble>;

// Method for conditional compilation based on type
#define IF_T_IS_DOUBLE(fma_version, ordinary_version) \
    do { \
        if constexpr (IS_DOUBLE<T>) { \
            fma_version; \
        } else { \
            ordinary_version; \
        } \
    } while(0)

// Compute MSB using bit manipulation
static constexpr unsigned msbslow(unsigned x) {
    unsigned m = 0;
    if (x != 0) {
        m = 1;
        while (x != 0) {
            x >>= 1;
            m <<= 1;
        }
        m >>= 1;
    }
    return m;
}

static unsigned msb(unsigned x) {
    if (x == 0) return 0;
    unsigned result = 0;
    unsigned shiftCount = 0;
    while (x > 255) {
        x >>= 8;
        shiftCount += 8;
    }
    result = static_cast<unsigned>(MSB_table[x]) << shiftCount;
    return result;
}

// Flip bits in the lower 48 bits
static unsigned long long bitcmp48(unsigned long long val) {
    static const unsigned long long mask48 = ((1ULL << 48) - 1ULL);
    val &= mask48;
    return (~val) & mask48;
}

// Helper for absolute value
template<typename T>
inline T abs_helper(const T& x) {
    if constexpr (IS_DOUBLEDOUBLE<T>) {
        return DD_ABS(x);
    } else {
        return std::abs(x);
    }
}

// mat2pm implementation
template<typename T>
std::vector<T> mat2pm(const std::vector<T>& A, int n, T tol, T thresh = T(-1)) {
    if (n <= 0 || n > 48) {
        std::cerr << "mat2pm: n out of range [1..48]\n";
        return {};
    }
    if (static_cast<int>(A.size()) != n*n) {
        std::cerr << "mat2pm: matrix size mismatch\n";
        return {};
    }

    // Compute scale
    T sumAbs;
    if constexpr (IS_DOUBLEDOUBLE<T>) {
        sumAbs = DoubleDouble(0);
        for (const T& v : A) sumAbs += DD_ABS(v);
    } else {
        sumAbs = T(0);
        for (const T& v : A) sumAbs += std::abs(v);
    }

    T scale;
    if constexpr (IS_DOUBLEDOUBLE<T>) {
        scale = sumAbs / DoubleDouble(n * n);
        if (scale == DoubleDouble(0)) scale = DoubleDouble(1);
    } else {
        scale = sumAbs / (n * n);
        if (scale == T(0)) scale = T(1);
    }

    T ppivot = scale;
    if constexpr (IS_DOUBLEDOUBLE<T>) {
        if (thresh < DoubleDouble(0)) {
            thresh = tol * scale;
        }
    } else {
        if (thresh < T(0)) {
            thresh = tol * scale;
        }
    }

    // Define eps for the type
    T eps_val;
    if constexpr (IS_DOUBLEDOUBLE<T>) {
        eps_val = DoubleDouble(1e-28);
    } else {
        eps_val = 1e-15;
    }

    if (thresh < eps_val) {
        thresh = eps_val;
    }

    size_t pmSize = (size_t(1) << n) - 1;
    std::vector<T> pm(pmSize);
    if constexpr (IS_DOUBLEDOUBLE<T>) {
        for (auto& val : pm) val = DoubleDouble(0);
    } else {
        for (auto& val : pm) val = T(0);
    }

    // Buffer for BFS submatrices
    std::vector<T> q;
    q.reserve(A.size());
    q.insert(q.end(), A.begin(), A.end());

    T pivmin;
    if constexpr (IS_DOUBLEDOUBLE<T>) {
        pivmin = DoubleDouble(std::numeric_limits<double>::infinity());
    } else {
        pivmin = std::numeric_limits<T>::infinity();
    }

    std::vector<size_t> zeropivs;
    size_t ipm = 1;

    // BFS over principal submatrices
    for (int level = 0; level < n; ++level) {
        int n1 = n - level;
        int nq = 1 << level;
        std::vector<T> qq((size_t)nq * 2 * (n1 - 1) * (n1 - 1));
        if constexpr (IS_DOUBLEDOUBLE<T>) {
            for (auto& val : qq) val = DoubleDouble(0);
        }

        const T* qPtr = q.data();
        T* qqPtrA = qq.data();
        T* qqPtrB = qq.data() + (size_t)nq * (n1 - 1) * (n1 - 1);
        size_t ipm1_base = 1;

        for (int i = 0; i < nq; ++i) {
            T pivotVal = qPtr[0];
            pm[ipm - 1] = pivotVal;

            if (n1 > 1) {
                T abspiv = abs_helper(pivotVal);
                if (abspiv <= thresh) {
                    zeropivs.push_back(ipm);
                    pm[ipm - 1] += ppivot;
                    abspiv = abs_helper(pm[ipm - 1]);
                }

                if constexpr (IS_DOUBLEDOUBLE<T>) {
                    if (abspiv < pivmin) pivmin = abspiv;
                } else {
                    pivmin = std::min(pivmin, abspiv);
                }

                int dimOut = n1 - 1;
                T* bOut = qqPtrA + (size_t)i * dimOut * dimOut;
                for (int rr = 1; rr < n1; ++rr) {
                    for (int cc = 1; cc < n1; ++cc) {
                        bOut[(rr-1)*dimOut + (cc-1)] = qPtr[rr*n1 + cc];
                    }
                }

                T* cOut = qqPtrB + (size_t)i * dimOut * dimOut;
                for (int rr = 1; rr < n1; ++rr) {
                    T dVal;
                    if (abs_helper(pm[ipm - 1]) < eps_val) {
                        if constexpr (IS_DOUBLEDOUBLE<T>) {
                            dVal = DoubleDouble(0);
                        } else {
                            dVal = T(0);
                        }
                    } else {
                        dVal = qPtr[rr * n1] / pm[ipm - 1];
                    }
                    for (int cc = 1; cc < n1; ++cc) {
                        T valB = bOut[(rr-1)*dimOut + (cc-1)];
                        T valA = qPtr[cc];
                        if constexpr (IS_DOUBLEDOUBLE<T>) {
                            cOut[(rr-1)*dimOut + (cc-1)] = valB - dVal * valA;
                        } else {
                            IF_T_IS_DOUBLE(
                                cOut[(rr-1)*dimOut + (cc-1)] = std::fma(-dVal, valA, valB),
                                cOut[(rr-1)*dimOut + (cc-1)] = valB - dVal * valA
                            );
                        }
                    }
                }
            }

            if (i > 0) {
                pm[ipm - 1] *= pm[ipm1_base - 1];
                ++ipm1_base;
            }

            ++ipm;
            qPtr += n1 * n1;
        }

        q.swap(qq);
    }

    // Correct any pseudo-pivots
    for (int idx = (int)zeropivs.size() - 1; idx >= 0; --idx) {
        size_t mask = zeropivs[idx];
        unsigned delta = msb((unsigned)mask);
        if (delta == 0) continue;
        unsigned delta2 = 2 * delta;

        unsigned long long maskLL = (unsigned long long)mask;
        unsigned long long valCmp = bitcmp48(delta);
        unsigned long long ipm1LL = maskLL & valCmp;
        size_t ipm1val = (size_t)ipm1LL;

        T pmMask = pm[mask - 1];
        if (ipm1val == 0ULL) {
            pm[mask - 1] = pmMask - ppivot;
        } else {
            T pmIPM1 = pm[ipm1val - 1];
            if constexpr (IS_DOUBLEDOUBLE<T>) {
                pm[mask - 1] = pmMask - ppivot * pmIPM1;
            } else {
                IF_T_IS_DOUBLE(
                   pm[mask - 1] = std::fma(-ppivot, pmIPM1, pmMask),
                   pm[mask - 1] = pmMask - ppivot * pmIPM1
                );
            }
        }

        for (size_t j = mask + delta2; j <= pmSize; j += delta2) {
            if constexpr (IS_DOUBLEDOUBLE<T>) {
                pm[j - 1] = pm[j - 1] - ppivot * pm[j - delta2 - 1];
            } else {
                IF_T_IS_DOUBLE(
                    pm[j - 1] = std::fma(-ppivot, pm[j - delta2 - 1], pm[j-1]),
                    pm[j - 1] = pm[j-1] - ppivot * pm[j - delta2 - 1]
                );
            }
        }
    }

    return pm;
}

// getpm function
template<typename T>
T getpm(const std::vector<T>& pm, const std::vector<int>& indexSet) {
    unsigned long long idx = 0ULL;
    for (int v : indexSet) {
        idx |= (1ULL << (v - 1));
    }
    if (idx == 0ULL) {
        if constexpr (IS_DOUBLEDOUBLE<T>) {
            return DoubleDouble(0);
        } else {
            return T(0);
        }
    }

    size_t pmIndex = (size_t)idx - 1;
    if (pmIndex >= pm.size()) {
        std::cerr << "getpm: index out of range.\n";
        if constexpr (IS_DOUBLEDOUBLE<T>) {
            return DoubleDouble(0);
        } else {
            return T(0);
        }
    }
    return pm[pmIndex];
}

// Helper function for popcnt
inline unsigned popcnt(unsigned x) {
    return __builtin_popcount(x);
}

#endif