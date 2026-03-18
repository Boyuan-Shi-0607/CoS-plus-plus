#ifndef INTEGRAND_NEW_H
#define INTEGRAND_NEW_H

#include "g_funcs.h"
#include "dag_generator.h"
#include "../utility/types.h"
#include "../utility/my_funcs.h"
#include <cmath>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <cstdint>
#include <cmath>

constexpr Real CYCLE_FACTOR = 3.0;
thread_local Real1DVector tl_vec1;
thread_local Real1DVector tl_M_matrix;

/**
 * integrand_connected_dd:
 *  - Inputs (all in [0,1)):
 *      [ τ'_0 .. τ'_{N-1},   R'_0,b1 .. R'_{N-1},b1,   R'_0,b2 .. R'_{N-1},b2 ]
 *  - τ' mapping:
 *      τ'_0 =  ε + u_0 * (β - 2ε)                  (positive)
 *      τ'_k = (2u_k - 1) * (β - 2ε),  k >= 1       (signed)
 *  - Enforce polyhedron S:
 *      Let S_k = ∑_{ℓ=1}^k τ'_ℓ (S_0=0),
 *          L = τ'_0 - (β - ε),  U = τ'_0 - ε.
 *      Require L < S_k < U for all k=1..N-1.
 *      If violated, return 0 immediately.
 *  - Δτ_{ij} = S_j - S_i.
 *  - ΔR from integer-difference sums of R' as before.
 */

inline Real wrapToInterval(Real x, Real beta_s) {
    Real period = 2.0 * beta_s;
    // fmod can give negative results, so normalize it
    Real y = std::fmod(x + beta_s, period);
    if (y < 0) y += period;
    return y - beta_s;
}

inline Real integrand_connected(const Real beta_s,
                                   const Real U_s,
                                   const Real* RESTRICT args,
                                   const Real1DVector& g_vec_flattened,
                                   const int order_num,
                                   const int N1,
                                   const int N2,
                                   const Real eps,
                                   const int cheby_degree,
                                   const Real alpha_shift,
                                   const DAGData* dag) {

    // Ensure sampler provided [0,1) for all inputs
    const int nvars = 3 * order_num;
    for (int k = 0; k < nvars; ++k)
        if (args[k] < 0.0 || args[k] >= 1.0) return 0.0;

    // Special case N=1 (Hartree only)
    if (order_num == 1) {
        const Real g_val = g_func_hartree_unified(g_vec_flattened, N1, N2, beta_s, eps, cheby_degree);
        const Real g_sq  = g_val * g_val;
        const Real term1 = g_sq * (-U_s / 2.0) * (CYCLE_FACTOR * CYCLE_FACTOR - CYCLE_FACTOR) * beta_s;
        const Real term2 = -alpha_shift * CYCLE_FACTOR * g_val * beta_s;
        return term1 + term2;
    }

    if (dag == nullptr) throw std::runtime_error("DAG is required for order_num > 1.");

    const Real scale = std::sqrt(U_s * Real(N1) * Real(N2) * beta_s / 2.0);
    const Real shift = alpha_shift / (U_s * CYCLE_FACTOR);

    // Unpack unit inputs
    const Real* RESTRICT tau_p_unit = args;
    const Real* RESTRICT Rp_b1_unit = &args[order_num];
    const Real* RESTRICT Rp_b2_unit = &args[2 * order_num];

    // Map τ' from unit to physical
    Real1DVector tau_prime(order_num, 0.0);
    tau_prime[0] = eps + tau_p_unit[0] * (beta_s - 2.0 * eps);                 // τ'_0 ∈ (ε, β-ε)
    // fix tau_prime[0] at eps
    // tau_prime[0] = eps;
    for (int k = 1; k < order_num; ++k)
        tau_prime[k] = (2.0 * tau_p_unit[k] - 1.0) * (beta_s - 2.0 * eps);     // τ'_k ∈ (-(β-2ε), β-2ε)

    // Polyhedral constraint S:
    // S_k = sum_{ℓ=1}^k τ'_ℓ,  L = τ'_0 - (β - ε),  U = τ'_0 - ε
    Real1DVector S_prefix(order_num, 0.0); // S_0 = 0
    for (int k = 1; k < order_num; ++k) {
        S_prefix[k] = S_prefix[k - 1] + tau_prime[k];
    }

    /*
    const Real L = tau_prime[0] - (beta_s - eps);
    const Real U = tau_prime[0] - eps;

    // Bounds check: L < S_k < U, for all k = 1..N-1
    // If any violates, this τ' lies outside S ⇒ 0 contribution.
    for (int k = 1; k < order_num; ++k) {
        if (!(S_prefix[k] > L && S_prefix[k] < U)) {
            return 0.0;
        }
    }
    */

    // Build V (order_num x order_num)
    const int vec1_size = order_num * order_num;
    if (tl_vec1.size() < static_cast<size_t>(vec1_size)) tl_vec1.resize(vec1_size);

    const Real g_hartree = g_func_hartree_unified(g_vec_flattened, N1, N2, beta_s, eps, cheby_degree);

    Real* RESTRICT V = tl_vec1.data();
    for (int i = 0; i < order_num; ++i) {
        Real* RESTRICT row = V + i * order_num;
        for (int j = 0; j < order_num; ++j) {
            if (i == j) {
                row[j] = g_hartree;
            } else {
                // ΔR reconstruction from R' unit inputs (unchanged)
                int R_b1_diff = 0;
                int R_b2_diff = 0;

                if (j > i) { // R_i - R_j = sum_{k=i+1..j} R'_k
                    for (int k = i + 1; k <= j; ++k) {
                        R_b1_diff += static_cast<int>(std::floor(Rp_b1_unit[k] * N1));
                        R_b2_diff += static_cast<int>(std::floor(Rp_b2_unit[k] * N2));
                    }
                } else { // j < i, R_i - R_j = - sum_{k=j+1..i} R'_k
                    for (int k = j + 1; k <= i; ++k) {
                        R_b1_diff += static_cast<int>(std::floor(Rp_b1_unit[k] * N1));
                        R_b2_diff += static_cast<int>(std::floor(Rp_b2_unit[k] * N2));
                    }
                    R_b1_diff = -R_b1_diff;
                    R_b2_diff = -R_b2_diff;
                }

                // Δτ from τ': τ_i - τ_j = S_j - S_i
                const Real tau_diff = S_prefix[j] - S_prefix[i];

                row[j] = g_func_unified_tau_diff(
                    wrapToInterval(tau_diff, beta_s),
                    R_b1_diff, R_b2_diff,
                    g_vec_flattened, beta_s, N1, N2,
                    eps, cheby_degree
                );
            }
        }
    }

    // Scale V and build M
    for (int k = 0; k < vec1_size; ++k) V[k] *= scale;

    const int M_size  = 2 * order_num;
    const int M_total = M_size * M_size;
    if (tl_M_matrix.size() < static_cast<size_t>(M_total)) tl_M_matrix.resize(M_total);
    Real* RESTRICT M = tl_M_matrix.data();

    for (int i = 0; i < order_num; ++i) {
        const Real* RESTRICT src = V + i * order_num;
        Real* RESTRICT dst_top    = M + i * M_size;
        Real* RESTRICT dst_bottom = M + (i + order_num) * M_size;
        std::memcpy(dst_top,                 src, sizeof(Real) * order_num);
        std::memcpy(dst_top + order_num,     src, sizeof(Real) * order_num);
        std::memcpy(dst_bottom,              src, sizeof(Real) * order_num);
        std::memcpy(dst_bottom + order_num,  src, sizeof(Real) * order_num);
    }

    // Add shift on the diagonal
    const Real diag_add = shift * scale;
    for (int i = 0; i < M_size; ++i) M[i * M_size + i] += diag_add;

    // Evaluate DAG sink
    const Real res_raw = compute_sink_value_fast(*dag, M, M_size);

    // Normalization + sign
    Real res = res_raw / (Real(N1) * Real(N2) * Real(factorial(order_num)));
    if (order_num & 1) res = -res;
    return res;
}

#endif // INTEGRAND_NEW_H

