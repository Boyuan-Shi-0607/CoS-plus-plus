#ifndef INTEGRAND_TYPED_H
#define INTEGRAND_TYPED_H

#include "g_funcs.h"
#include "dag_generator.h"
#include "../utility/types.h"
#include "../utility/my_funcs.h"

#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>

// We keep thread-locals for per-flavour dense V matrices (flattened row-major)
// tl_vec1_multi[f] has size order_num*order_num and stores V_f
thread_local std::vector<Real1DVector> tl_vec1_multi;

inline Real wrapToInterval(Real x, Real beta_s) {
    Real period = 2.0 * beta_s;
    Real y = std::fmod(x + beta_s, period); // fmod can be negative
    if (y < 0) y += period;
    return y - beta_s;
}

/**
 * integrand_connected_typed:
 *  - Multi-flavour version that mirrors the typed DAG.
 *  - Inputs (all in [0,1), interleaved per vertex):
 *      [ τ'_0, R'_{0,b1}, R'_{0,b2},  τ'_1, R'_{1,b1}, R'_{1,b2},  ... , τ'_{N-1}, R'_{N-1,b1}, R'_{N-1,b2} ]
 *  - g_vecs_flattened: size N_f, each flavour is a flattened (kx,ky,τ) Green's table (opaque to this function).
 *  - alpha_shift_vec:  size N_f, per-flavour shift parameters.
 *  - dag: typed DAG built with build_dag_typed(order_num, N_f).
 */

inline Real integrand_connected_typed(const Real beta_s,
                                     const Real U_s,
                                     const Real* RESTRICT args,
                                     const std::vector<Real1DVector>& g_vecs_flattened, // N_f
                                     const int order_num,
                                     const int N1,
                                     const int N2,
                                     const Real eps,
                                     const int cheby_degree,
                                     const std::vector<Real>& alpha_shift_vec,          // N_f
                                     const typed_dag::DAGData* dag) {

    const int N_f = static_cast<int>(g_vecs_flattened.size());
    if (N_f <= 0) throw std::runtime_error("N_f must be > 0");
    if (static_cast<int>(alpha_shift_vec.size()) != N_f)
        throw std::runtime_error("alpha_shift_vec size must equal number of flavours (g_vecs_flattened.size())");

    // Ensure sampler provided [0,1) for all inputs
    const int nvars = 3 * order_num;
    for (int k = 0; k < nvars; ++k)
        if (args[k] < 0.0 || args[k] >= 1.0) return 0.0;

    // Helpers to access interleaved unit inputs
    auto tau_u = [&](int i)->Real { return args[3*i + 0]; };
    auto r1_u  = [&](int i)->Real { return args[3*i + 1]; };
    auto r2_u  = [&](int i)->Real { return args[3*i + 2]; };

    // Map τ' from unit to physical
    Real1DVector tau_prime(order_num, 0.0);
    tau_prime[0] = eps; // fix τ'_0 at ε (as in previous integrand)
    for (int k = 1; k < order_num; ++k)
        tau_prime[k] = (2.0 * tau_u(k) - 1.0) * (beta_s - 2.0 * eps);    // τ'_k ∈ (-(β-2ε), β-2ε)

    // Prefix sums S_k (S_0 = 0) so that Δτ_{ij} = S_j - S_i
    Real1DVector S_prefix(order_num, 0.0);
    for (int k = 1; k < order_num; ++k)
        S_prefix[k] = S_prefix[k - 1] + tau_prime[k];

    // Global scale (same for all flavours)
    const Real scale = std::sqrt(U_s * Real(N1) * Real(N2) * beta_s / 2.0);

    // ========== Special case N = 1 ==========
    if (order_num == 1) {
        // Per-flavour g_val (Hartree)
        std::vector<Real> g_vals(N_f);
        for (int f = 0; f < N_f; ++f)
            g_vals[f] = g_func_hartree_unified(g_vecs_flattened[f], N1, N2, beta_s, eps, cheby_degree);

        // term 1: sum over i != j (ordered pairs) of g_i * g_j
        Real S1 = 0.0, S2 = 0.0;
        for (int f = 0; f < N_f; ++f) { S1 += g_vals[f]; S2 += g_vals[f]*g_vals[f]; }
        const Real sum_i_neq_j = S1*S1 - S2; // equals N_f(N_f-1) terms of g_i g_j
        const Real term1 = sum_i_neq_j * (-U_s / 2.0) * beta_s;

        // term 2: sum over flavours with matching alpha index (one-to-one): − alpha_f * CYCLE_FACTOR * g_f * β
        Real term2 = 0.0;
        for (int f = 0; f < N_f; ++f) term2 += (-alpha_shift_vec[f]) * g_vals[f] * beta_s;

        return term1 + term2;
    }

    if (dag == nullptr) throw std::runtime_error("Typed DAG is required for order_num > 1.");

    // Ensure thread-local buffers
    if (static_cast<int>(tl_vec1_multi.size()) < N_f)
        tl_vec1_multi.resize(N_f);

    const int n = order_num;                 // matrix side per flavour for typed DAG
    const int mat_elems = n * n;             // flattened size per flavour

    // Build per-flavour dense V matrices (row-major), scaled by `scale`
    for (int f = 0; f < N_f; ++f) {
        auto &Vf = tl_vec1_multi[f];
        if (static_cast<int>(Vf.size()) < mat_elems) Vf.resize(mat_elems);
        const Real g_hartree_f = g_func_hartree_unified(g_vecs_flattened[f], N1, N2, beta_s, eps, cheby_degree);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Real val;
                if (i == j) {
                    val = g_hartree_f;
                } else {
                    int R_b1_diff = 0, R_b2_diff = 0;
                    if (j > i) {
                        for (int k = i + 1; k <= j; ++k) {
                            R_b1_diff += static_cast<int>(std::floor(r1_u(k) * N1));
                            R_b2_diff += static_cast<int>(std::floor(r2_u(k) * N2));
                        }
                    } else { // j < i
                        for (int k = j + 1; k <= i; ++k) {
                            R_b1_diff += static_cast<int>(std::floor(r1_u(k) * N1));
                            R_b2_diff += static_cast<int>(std::floor(r2_u(k) * N2));
                        }
                        R_b1_diff = -R_b1_diff;
                        R_b2_diff = -R_b2_diff;
                    }
                    const Real tau_diff = S_prefix[j] - S_prefix[i];
                    val = g_func_unified_tau_diff(
                        wrapToInterval(tau_diff, beta_s),
                        R_b1_diff, R_b2_diff,
                        g_vecs_flattened[f], beta_s, N1, N2,
                        eps, cheby_degree
                    );
                }
                Vf[i * n + j] = val * scale; // apply scale here
            }
        }
    }

    // Build per-flavour shift close values, scaled consistently
    std::vector<Real> shift_close(N_f);
    for (int f = 0; f < N_f; ++f) {
        shift_close[f] = scale * (alpha_shift_vec[f] / Real(2.0)) / (U_s / Real(2.0));
    }

    // Typed fast sink evaluation
    const Real res_raw = typed_dag::compute_sink_value_fast(*dag, tl_vec1_multi, n, shift_close);

    // Normalization + sign, kept as before
    Real res = res_raw / (Real(N1) * Real(N2) * Real(factorial(order_num)));
    if (order_num & 1) res = -res;
    return res;
}

#endif // INTEGRAND_TYPED_H
