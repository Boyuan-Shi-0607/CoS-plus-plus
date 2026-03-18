#ifndef BIPARTITE_LATTICE_PRINCIPAL_MINOR_CDET_G_FUNCS_H
#define BIPARTITE_LATTICE_PRINCIPAL_MINOR_CDET_G_FUNCS_H
#define p_modulo(i, n) (((i) % (n) + (n)) % (n))

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <complex>
#include <numeric>
#include <type_traits>
#include "../utility/types.h"

// Check once at compile time if Real is double
static constexpr bool REAL_IS_DOUBLE = std::is_same_v<Real, double>;

// Method 2: For statements (like return statements)
#define IF_REAL_IS_DOUBLE(fma_version, ordinary_version) \
    do { \
        if constexpr (REAL_IS_DOUBLE) { \
            fma_version; \
        } else { \
            ordinary_version; \
        } \
    } while(0)

// Helper function to convert 5D indices to 1D index - optimized with precomputed strides
inline size_t index_5D_to_1D(const int k, const int a, const int b,
    const int R_b1_mod, const int R_b2_mod, const int N1, const int N2,
    int cheby_degree) {

    // Optimized: compute strides only once
    const size_t stride_R_b1 = N2;
    const size_t stride_b = N1 * N2;
    const size_t stride_a = 2 * stride_b;
    const size_t stride_k = 2 * stride_a;

    return k * stride_k + a * stride_a + b * stride_b + R_b1_mod * stride_R_b1 + R_b2_mod;
}

// DoubleDouble version of ChebyshevSum
inline DoubleDouble ChebyshevSum_dd(const int n, const DoubleDouble& x, const DD1DVector& g_vec_flattened,
                                   const int a, const int b, const int R_b1_mod, const int R_b2_mod,
                                   const int N1, const int N2, const int cheby_degree) {
    // Precompute strides for index calculation
    const size_t stride_R_b1 = N2;
    const size_t stride_b = N1 * N2;
    const size_t stride_a = 2 * stride_b;
    const size_t stride_k = 2 * stride_a;

    // Base index for this (a,b,R_b1_mod,R_b2_mod) combination
    const size_t base_idx = a * stride_a + b * stride_b + R_b1_mod * stride_R_b1 + R_b2_mod;

    // Special cases
    if (n == 0) {
        return g_vec_flattened[base_idx];
    }
    if (n == 1) {
        return g_vec_flattened[base_idx] + g_vec_flattened[base_idx + stride_k] * x;
    }

    // Clenshaw algorithm
    const DoubleDouble two_x = DoubleDouble(2) * x;
    DoubleDouble b_n_plus_1(0);
    DoubleDouble b_n_plus_2(0);
    DoubleDouble b_k;

    for (int k = n; k >= 1; --k) {
        b_k = g_vec_flattened[base_idx + k * stride_k] + two_x * b_n_plus_1 - b_n_plus_2;
        b_n_plus_2 = b_n_plus_1;
        b_n_plus_1 = b_k;
    }

    return g_vec_flattened[base_idx] + b_n_plus_1 * x - b_n_plus_2;
}

// DoubleDouble version of g_cheby_pos_flattened
inline DoubleDouble g_cheby_pos_flattened_dd(const int a, const int b, const DoubleDouble& tau, const DoubleDouble& eps,
                                            const int N1, const int N2, const int R_b1, const int R_b2,
                                            const DD1DVector& g_vec_flattened,
                                            const DoubleDouble& beta_s,
                                            const int degree,
                                            const int cheby_degree) {
    // Optimized t calculation
    const DoubleDouble inv_beta_range = DoubleDouble(1) / (beta_s - DoubleDouble(2) * eps);
    DoubleDouble t = DoubleDouble(2) * (tau - eps) * inv_beta_range - DoubleDouble(1);

    const int R_b1_mod = p_modulo(R_b1, N1);
    const int R_b2_mod = p_modulo(R_b2, N2);

    return ChebyshevSum_dd(degree, t, g_vec_flattened, a, b, R_b1_mod, R_b2_mod, N1, N2, cheby_degree);
}

// DoubleDouble version of g_func_flattened
inline DoubleDouble g_func_flattened_dd(const int a, const int b,
                                       const DoubleDouble& tau_1, const DoubleDouble& tau_2,
                                       const DoubleDouble& R1_b1, const DoubleDouble& R1_b2,
                                       const DoubleDouble& R2_b1, const DoubleDouble& R2_b2,
                                       const DD1DVector& g_vec_flattened,
                                       const DoubleDouble& beta_s,
                                       const int N1, const int N2, const DoubleDouble& eps,
                                       const int degree,
                                       const int cheby_degree) {

    const int R_b1 = DD_FLOOR_TO_INT(R1_b1) - DD_FLOOR_TO_INT(R2_b1);
    const int R_b2 = DD_FLOOR_TO_INT(R1_b2) - DD_FLOOR_TO_INT(R2_b2);
    const DoubleDouble tau_diff = tau_1 - tau_2;

    const DoubleDouble beta_s_minus_eps = beta_s - eps;
    const DoubleDouble minus_beta_s_plus_eps = -beta_s + eps;

    DoubleDouble val(0);

    if (tau_diff > eps && tau_diff < beta_s_minus_eps) {
        val = g_cheby_pos_flattened_dd(a, b, tau_diff, eps, N1, N2, R_b1, R_b2, g_vec_flattened, beta_s, degree, cheby_degree);
    }
    else if (tau_diff > minus_beta_s_plus_eps && tau_diff < -eps) {
        val = -g_cheby_pos_flattened_dd(a, b, tau_diff + beta_s, eps, N1, N2, R_b1, R_b2, g_vec_flattened, beta_s, degree, cheby_degree);
    }

    return val;
}

// DoubleDouble version of g_func_hartree_flattened
inline DoubleDouble g_func_hartree_flattened_dd(const int a, const int b,
                                               const DoubleDouble& R1_b1, const DoubleDouble& R1_b2,
                                               const DoubleDouble& R2_b1, const DoubleDouble& R2_b2,
                                               const DD1DVector& g_vec_flattened,
                                               const int N1, const int N2, const DoubleDouble& beta_s,
                                               const DoubleDouble& eps, const int degree,
                                               const int cheby_degree) {
    const int R_b1 = DD_FLOOR_TO_INT(R1_b1) - DD_FLOOR_TO_INT(R2_b1);
    const int R_b2 = DD_FLOOR_TO_INT(R1_b2) - DD_FLOOR_TO_INT(R2_b2);

    return -g_cheby_pos_flattened_dd(a, b, beta_s - eps, eps, N1, N2, R_b1, R_b2, g_vec_flattened, beta_s, degree, cheby_degree);
}

// Function to calculate factorial
inline unsigned long long factorial(int n) {
    if (n < 0) {
        std::cerr << "Factorial is not defined for negative numbers." << std::endl;
        return 0;
    }
    if (n > 11) {
        throw std::runtime_error("n > 11; not allowed.");
    }
    unsigned long long result = 1;
    for (int i = 1; i <= n; ++i) {
        result *= i;
    }
    return result;
}

inline Int1DVector findIndices(int subset_mask, int order_num) {
    Int1DVector indices;
    for (size_t i = 0; i < order_num; ++i) {
        if (subset_mask & (1 << i)) {
            indices.push_back(i);
        }
    }
    return indices;
}

inline Int1DVector findIndices_greens_func(int subset_mask, int order_num) {
    Int1DVector indices;
    for (size_t i = 0; i < order_num+1; ++i) {
        if (subset_mask & (1 << i)) {
            indices.push_back(i);
        }
    }
    return indices;
}

#endif