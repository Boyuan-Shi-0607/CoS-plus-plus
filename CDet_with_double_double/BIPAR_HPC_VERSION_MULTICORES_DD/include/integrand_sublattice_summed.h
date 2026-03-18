#ifndef INTEGRAND_SUBLATTICE_SUMMED_H
#define INTEGRAND_SUBLATTICE_SUMMED_H

#include "fpm_sublattices_summed.h"
#include "fast_principle_minors.h"  // For popcnt and other utilities
#include "g_funcs.h"
#include "../utility/my_funcs.h"  // For is_finite
#include <vector>
#include <cmath>
#include <stdexcept>
#include <cstdlib>  // For aligned_alloc
#include <sstream>  // For stringstream

#include <memory>

// Global workspace manager for FPM computations
struct FPMWorkspaceManager {
    static constexpr size_t WORKSPACE_SIZE = FPM_WORKSPACE_SIZE;
    static constexpr size_t OUT_RAW_SIZE = FPM_OUT_RAW_SIZE;
    static constexpr size_t MINORS_SIZE = (1 << FPM_MAXN) - 1; // 2^12 - 1

    // Thread-local storage for workspace buffers
    struct ThreadLocalBuffers {
        std::vector<DoubleDouble> dd_workspace;
        std::vector<DoubleDouble> dd_out_raw_up;
        std::vector<DoubleDouble> dd_out_raw_down;
        std::vector<DoubleDouble> dd_minors_result;

        std::vector<double> d_workspace;
        std::vector<double> d_out_raw_up;
        std::vector<double> d_out_raw_down;
        std::vector<double> d_minors_result;

        ThreadLocalBuffers() {
            // Pre-allocate all buffers
            dd_workspace.resize(2 * WORKSPACE_SIZE);
            dd_out_raw_up.resize(OUT_RAW_SIZE);
            dd_out_raw_down.resize(OUT_RAW_SIZE);
            dd_minors_result.resize(MINORS_SIZE);

            d_workspace.resize(2 * WORKSPACE_SIZE);
            d_out_raw_up.resize(OUT_RAW_SIZE);
            d_out_raw_down.resize(OUT_RAW_SIZE);
            d_minors_result.resize(MINORS_SIZE);
        }
    };

    static thread_local ThreadLocalBuffers buffers;

    // Get buffers for DoubleDouble computation
    static void get_dd_buffers(DoubleDouble*& workspace, DoubleDouble*& out_raw_up,
                              DoubleDouble*& out_raw_down, DoubleDouble*& minors_result) {
        workspace = buffers.dd_workspace.data();
        out_raw_up = buffers.dd_out_raw_up.data();
        out_raw_down = buffers.dd_out_raw_down.data();
        minors_result = buffers.dd_minors_result.data();
    }

    // Get buffers for double computation
    static void get_double_buffers(double*& workspace, double*& out_raw_up,
                                  double*& out_raw_down, double*& minors_result) {
        workspace = buffers.d_workspace.data();
        out_raw_up = buffers.d_out_raw_up.data();
        out_raw_down = buffers.d_out_raw_down.data();
        minors_result = buffers.d_minors_result.data();
    }
};

// Define thread-local storage
thread_local FPMWorkspaceManager::ThreadLocalBuffers FPMWorkspaceManager::buffers;

// Modified DP Z-transform for pre-multiplied minors
template<typename T>
T compute_dp_ztransform_premultiplied(const T* minors, const int order_num) {
    const size_t N = size_t{1} << order_num;

    // Aligned alloc for a[] and dp[]
    size_t bytes = ((N * sizeof(T)) + 31) & ~size_t(31);
    auto *a = (T*)std::aligned_alloc(32, bytes);
    auto *dp = (T*)std::aligned_alloc(32, bytes);

    // Build a array - minors are already products of up and down
    if constexpr (std::is_same_v<T, DoubleDouble>) {
        a[0] = DoubleDouble(0);
    } else {
        a[0] = T(0);
    }

    // Copy pre-multiplied minors to a[]
    for (size_t m = 1; m < N; ++m) {
        a[m] = minors[m - 1];  // minors array is 0-indexed
    }

    // Initialize dp array
    for (size_t i = 0; i < N; ++i) {
        if constexpr (std::is_same_v<T, DoubleDouble>) {
            dp[i] = DoubleDouble(0);
        } else {
            dp[i] = T(0);
        }
    }

    // Base cases
    for (int i = 0; i < order_num; ++i) {
        dp[size_t{1} << i] = a[size_t{1} << i];
    }

    // Bucket by popcount
    std::vector<std::vector<size_t>> by_cnt(order_num + 1);
    for (size_t m = 0; m < N; ++m) {
        by_cnt[__builtin_popcount(static_cast<uint32_t>(m))].push_back(m);
    }

    // DP computation
    for (int bits = 2; bits <= order_num; ++bits) {
        for (auto mask : by_cnt[bits]) {
            T acc = a[mask];
            if (mask & 1) {
                for (size_t sub = (mask - 1) & mask; sub; sub = (sub - 1) & mask) {
                    if (!(sub & 1)) continue;
                    acc = acc - dp[sub] * a[mask ^ sub];
                }
            }
            dp[mask] = acc;
        }
    }

    T res = dp[N - 1];

    std::free(a);
    std::free(dp);

    return res;
}

// DIVISION OPTIMIZATION: Pre-computed constants structure
template<typename T>
struct IntegrandConstants {
    T inv_U_s;
    T normalization_factor;

    IntegrandConstants(const T& U_s, int order_num, int N1, int N2) {
        inv_U_s = T(1) / U_s;

        T denominator = T(factorial(order_num)) * T(2) * T(N1) * T(N2);
        normalization_factor = T(1) / denominator;
    }
};

// Sublattice-summed integrand function using DoubleDouble
inline DoubleDouble integrand_sublattice_summed_dd(const DoubleDouble beta_s,
                                                  const DoubleDouble U_s,
                                                  const DoubleDouble* args,
                                                  const DD1DVector& g_vec_1_flattened,
                                                  const DD1DVector& g_vec_2_flattened,
                                                  const int order_num,
                                                  const int N1,
                                                  const int N2,
                                                  const DoubleDouble eps_val,
                                                  const int cheby_deg,
                                                  const DoubleDouble alpha_up,
                                                  const DoubleDouble alpha_down,
                                                  const DoubleDouble delta_up,
                                                  const DoubleDouble delta_down,
                                                  const DoubleDouble tol) {

    // Check bounds for 3*order_num arguments
    for (int i = 0; i < 3 * order_num; ++i) {
        if (args[i] < DoubleDouble(0) || args[i] > DoubleDouble(1)) {
            return DoubleDouble(0);
        }
    }

    // Pre-compute expensive reciprocals once
    IntegrandConstants<DoubleDouble> constants(U_s, order_num, N1, N2);

    // Prepare 2n × 2n matrices
    const int dim = 2 * order_num;
    DD1DVector vec1(dim * dim, DoubleDouble(0));
    DD1DVector vec2(dim * dim, DoubleDouble(0));

    DoubleDouble scale = DD_SQRT(U_s * DoubleDouble(N1) * DoubleDouble(N2) * beta_s);

    // Fill the matrices with 2×2 blocks
    for (int i = 0; i < order_num; ++i) {
        for (int j = 0; j < order_num; ++j) {
            int base_row = 2 * i;
            int base_col = 2 * j;

            if (i == j) {
                // Hartree term - diagonal blocks
                for (int m = 0; m < 2; ++m) {
                    for (int n = 0; n < 2; ++n) {
                        int row = base_row + m;
                        int col = base_col + n;
                        int idx = row * dim + col;

                        vec1[idx] = g_func_hartree_flattened_dd(m, n,
                            DoubleDouble(0), DoubleDouble(0), DoubleDouble(0), DoubleDouble(0),
                            g_vec_1_flattened, N1, N2, beta_s, eps_val, cheby_deg, cheby_deg);
                        vec2[idx] = g_func_hartree_flattened_dd(m, n,
                            DoubleDouble(0), DoubleDouble(0), DoubleDouble(0), DoubleDouble(0),
                            g_vec_2_flattened, N1, N2, beta_s, eps_val, cheby_deg, cheby_deg);

                        if (m == n) {
                            DoubleDouble shift_down = alpha_down + ((m == 0) ? delta_down : -delta_down);
                            DoubleDouble shift_up = alpha_up + ((m == 0) ? delta_up : -delta_up);

                            vec1[idx] += shift_down * constants.inv_U_s;
                            vec2[idx] += shift_up * constants.inv_U_s;
                        }

                        vec1[idx] *= scale;
                        vec2[idx] *= scale;
                    }
                }
            } else {
                // Off-diagonal blocks
                DoubleDouble tau1 = beta_s * args[3 * i];
                DoubleDouble tau2 = beta_s * args[3 * j];
                DoubleDouble R1b1 = DoubleDouble(N1) * args[3 * i + 1];
                DoubleDouble R1b2 = DoubleDouble(N2) * args[3 * i + 2];
                DoubleDouble R2b1 = DoubleDouble(N1) * args[3 * j + 1];
                DoubleDouble R2b2 = DoubleDouble(N2) * args[3 * j + 2];

                for (int m = 0; m < 2; ++m) {
                    for (int n = 0; n < 2; ++n) {
                        int row = base_row + m;
                        int col = base_col + n;
                        int idx = row * dim + col;

                        vec1[idx] = g_func_flattened_dd(m, n,
                                                       tau1, tau2, R1b1, R1b2, R2b1, R2b2,
                                                       g_vec_1_flattened, beta_s, N1, N2,
                                                       eps_val, cheby_deg, cheby_deg);
                        vec2[idx] = g_func_flattened_dd(m, n,
                                                       tau1, tau2, R1b1, R1b2, R2b1, R2b2,
                                                       g_vec_2_flattened, beta_s, N1, N2,
                                                       eps_val, cheby_deg, cheby_deg);

                        vec1[idx] *= scale;
                        vec2[idx] *= scale;
                    }
                }
            }
        }
    }

    // Check for NaN
    check_vector_nan(vec1, "vec1");
    check_vector_nan(vec2, "vec2");

    // Get pre-allocated buffers
    DoubleDouble *workspace, *out_raw_up, *out_raw_down, *minors_premultiplied;
    FPMWorkspaceManager::get_dd_buffers(workspace, out_raw_up, out_raw_down, minors_premultiplied);

    // Use the pre-allocated buffers for FPM computation
    compute_postprocessed_minors(vec1, vec2, order_num, workspace, out_raw_up, out_raw_down, minors_premultiplied);

    // Check for NaN in result
    for (int i = 0; i < pow2_table[order_num] - 1; ++i) {
        if (!is_finite(minors_premultiplied[i])) {
            throw std::runtime_error("NaN detected in minors_premultiplied");
        }
    }

    // Use modified DP transform for pre-multiplied minors
    auto c_result = compute_dp_ztransform_premultiplied<DoubleDouble>(minors_premultiplied, order_num);

    // Use pre-computed constants
    DoubleDouble pow_1;
    if (order_num % 2 == 0) {
      pow_1 = DoubleDouble(1);
    } else {
      pow_1 = DoubleDouble(-1);
    }
    const DoubleDouble val = pow_1 * c_result * constants.normalization_factor;

    if (!is_finite(val)) {
        throw std::runtime_error("Error (integrand): Value is not finite");
    }
    return val;
}

// Template version for other types
template<typename T>
T integrand_sublattice_summed(const T beta_s,
                             const T U_s,
                             const T* args,
                             const std::vector<T>& g_vec_1_flattened,
                             const std::vector<T>& g_vec_2_flattened,
                             const int order_num,
                             const int N1,
                             const int N2,
                             const T eps_val,
                             const int cheby_deg,
                             const T alpha_up,
                             const T alpha_down,
                             const T delta_up,
                             const T delta_down,
                             const T tol) {

    // Check bounds
    for (int i = 0; i < 3 * order_num; ++i) {
        if (args[i] < T(0) || args[i] > T(1)) {
            return T(0);
        }
    }

    // Pre-compute expensive reciprocals
    IntegrandConstants<T> constants(U_s, order_num, N1, N2);

    // Prepare 2n × 2n matrices
    const int dim = 2 * order_num;
    std::vector<T> vec1(dim * dim, T(0));
    std::vector<T> vec2(dim * dim, T(0));

    T scale = std::sqrt(U_s * T(N1) * T(N2) * beta_s);

    // Fill matrices
    for (int i = 0; i < order_num; ++i) {
        for (int j = 0; j < order_num; ++j) {
            int base_row = 2 * i;
            int base_col = 2 * j;

            if (i == j) {
                // Hartree term
                for (int m = 0; m < 2; ++m) {
                    for (int n = 0; n < 2; ++n) {
                        int idx = (base_row + m) * dim + (base_col + n);

                        vec1[idx] = g_func_hartree_flattened(m, n,
                            T(0), T(0), T(0), T(0),
                            g_vec_1_flattened, N1, N2, beta_s, eps_val, cheby_deg, cheby_deg);
                        vec2[idx] = g_func_hartree_flattened(m, n,
                            T(0), T(0), T(0), T(0),
                            g_vec_2_flattened, N1, N2, beta_s, eps_val, cheby_deg, cheby_deg);

                        if (m == n) {
                            T shift_down = alpha_down + ((m == 0) ? delta_down : -delta_down);
                            T shift_up = alpha_up + ((m == 0) ? delta_up : -delta_up);

                            vec1[idx] += shift_down * constants.inv_U_s;
                            vec2[idx] += shift_up * constants.inv_U_s;
                        }

                        vec1[idx] *= scale;
                        vec2[idx] *= scale;
                    }
                }
            } else {
                // Off-diagonal blocks
                T tau1 = beta_s * args[3 * i];
                T tau2 = beta_s * args[3 * j];
                T R1b1 = T(N1) * args[3 * i + 1];
                T R1b2 = T(N2) * args[3 * i + 2];
                T R2b1 = T(N1) * args[3 * j + 1];
                T R2b2 = T(N2) * args[3 * j + 2];

                for (int m = 0; m < 2; ++m) {
                    for (int n = 0; n < 2; ++n) {
                        int idx = (base_row + m) * dim + (base_col + n);

                        vec1[idx] = g_func_flattened(m, n,
                                                    tau1, tau2, R1b1, R1b2, R2b1, R2b2,
                                                    g_vec_1_flattened, beta_s, N1, N2,
                                                    eps_val, cheby_deg, cheby_deg) * scale;
                        vec2[idx] = g_func_flattened(m, n,
                                                    tau1, tau2, R1b1, R1b2, R2b1, R2b2,
                                                    g_vec_2_flattened, beta_s, N1, N2,
                                                    eps_val, cheby_deg, cheby_deg) * scale;
                    }
                }
            }
        }
    }

    // Check for NaN
    check_vector_nan(vec1, "vec1");
    check_vector_nan(vec2, "vec2");

    // Get pre-allocated buffers (for double type)
    double *workspace, *out_raw_up, *out_raw_down, *minors_premultiplied;
    FPMWorkspaceManager::get_double_buffers(workspace, out_raw_up, out_raw_down, minors_premultiplied);

    // Use pre-allocated buffers
    compute_postprocessed_minors(vec1, vec2, order_num, workspace, out_raw_up, out_raw_down, minors_premultiplied);

    // Check for NaN
    for (int i = 0; i < pow2_table[order_num] - 1; ++i) {
        if (!std::isfinite(minors_premultiplied[i])) {
            throw std::runtime_error("NaN detected in minors_premultiplied");
        }
    }

    // Use modified DP transform
    auto c_result = compute_dp_ztransform_premultiplied<T>(minors_premultiplied, order_num);

    // Use pre-computed normalization factor
    const T val = std::pow(-1.0, order_num) * c_result * constants.normalization_factor;

    if (!std::isfinite(val)) {
        throw std::runtime_error("Error (integrand): Value is not finite");
    }
    return val;
}

#endif