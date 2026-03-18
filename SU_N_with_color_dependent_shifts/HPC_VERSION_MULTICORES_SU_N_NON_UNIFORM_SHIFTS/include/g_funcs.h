#ifndef G_FUNCS_UNIFIED_H
#define G_FUNCS_UNIFIED_H

#include <iostream>
#include <vector>
#include <cmath>
#include "../utility/types.h" // Real, Real1DVector, Int1DVector, etc.

#define p_modulo(i, n) (((i) % (n) + (n)) % (n))

// ----------------- Chebyshev evaluator -----------------

inline Real ChebyshevSum(const int n, const Real& x, const Real1DVector& g_vec_flattened,
                         const int R_b1_mod, const int R_b2_mod,
                         const int N1, const int N2) {
    const size_t stride_R_b1 = N2;
    const size_t stride_k    = static_cast<size_t>(N1) * static_cast<size_t>(N2);
    const size_t base_idx    = static_cast<size_t>(R_b1_mod) * stride_R_b1 + static_cast<size_t>(R_b2_mod);

    if (n == 0) return g_vec_flattened[base_idx];
    if (n == 1) return g_vec_flattened[base_idx] + g_vec_flattened[base_idx + stride_k] * x;

    const Real two_x = 2.0 * x;
    Real b_np1 = 0.0, b_np2 = 0.0, b_k;
    for (int k = n; k >= 1; --k) {
        b_k  = g_vec_flattened[base_idx + static_cast<size_t>(k) * stride_k] + two_x * b_np1 - b_np2;
        b_np2 = b_np1;
        b_np1 = b_k;
    }
    return g_vec_flattened[base_idx] + b_np1 * x - b_np2;
}

inline Real g_cheby_pos_flattened(const Real& tau_pos, const Real& eps,
                                  const int N1, const int N2, const int R_b1, const int R_b2,
                                  const Real1DVector& g_vec_flattened,
                                  const Real& beta_s,
                                  const int degree) {
    // Map tau ∈ (eps, beta - eps) to t ∈ (-1, 1)
    const Real inv = 1.0 / (beta_s - 2.0 * eps);
    const Real t   = 2.0 * (tau_pos - eps) * inv - 1.0;
    const int R_b1_mod = p_modulo(R_b1, N1);
    const int R_b2_mod = p_modulo(R_b2, N2);
    return ChebyshevSum(degree, t, g_vec_flattened, R_b1_mod, R_b2_mod, N1, N2);
}

// ----------------- Unified Green's functions -----------------

/**
 * g(Δτ, ΔR_b1, ΔR_b2) for off-diagonal entries.
 * Anti-periodicity in τ:
 *   if Δτ ∈ (  eps,  β-eps):  +g(Δτ)
 *   if Δτ ∈ (-β+eps, -eps):   -g(Δτ+β)
 *   else: 0 (exclude ε-slivers)
 */
inline Real g_func_unified_tau_diff(const Real& tau_diff,
                                    const int R_b1_diff, const int R_b2_diff,
                                    const Real1DVector& g_vec_flattened,
                                    const Real& beta_s,
                                    const int N1, const int N2, const Real& eps,
                                    const int degree) {
    const Real beta_minus_eps   = beta_s - eps;
    const Real minus_beta_plus_eps = -beta_s + eps;

    if (tau_diff > eps && tau_diff < beta_minus_eps) {
        return g_cheby_pos_flattened(tau_diff, eps, N1, N2,
                                     R_b1_diff, R_b2_diff,
                                     g_vec_flattened, beta_s, degree);
    }
    if (tau_diff > minus_beta_plus_eps && tau_diff < -eps) {
        return -g_cheby_pos_flattened(tau_diff + beta_s, eps, N1, N2,
                                      R_b1_diff, R_b2_diff,
                                      g_vec_flattened, beta_s, degree);
    }
    return 0.0;
}

/** Hartree (diagonal) value: -g(β-ε, 0) */
inline Real g_func_hartree_unified(const Real1DVector& g_vec_flattened,
                                   const int N1, const int N2, const Real& beta_s,
                                   const Real& eps, const int degree) {
    return -g_cheby_pos_flattened(beta_s - eps, eps, N1, N2, 0, 0,
                                  g_vec_flattened, beta_s, degree);
}

// ----------------- Misc utilities -----------------

inline unsigned long long factorial(int n) {
    if (n < 0) {
        std::cerr << "Factorial is not defined for negative numbers.\n";
        return 0ULL;
    }
    if (n > 11) throw std::runtime_error("n > 11; not allowed.");
    unsigned long long r = 1ULL;
    for (int i = 1; i <= n; ++i) r *= static_cast<unsigned long long>(i);
    return r;
}

inline Int1DVector findIndices(int subset_mask, int order_num) {
    Int1DVector idx;
    for (int i = 0; i < order_num; ++i)
        if (subset_mask & (1 << i)) idx.push_back(i);
    return idx;
}

inline Int1DVector findIndices_greens_func(int subset_mask, int order_num) {
    Int1DVector idx;
    for (int i = 0; i < order_num + 1; ++i)
        if (subset_mask & (1 << i)) idx.push_back(i);
    return idx;
}

#endif // G_FUNCS_UNIFIED_H
