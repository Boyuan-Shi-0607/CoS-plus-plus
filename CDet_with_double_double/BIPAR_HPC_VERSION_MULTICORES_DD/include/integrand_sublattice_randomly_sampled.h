#ifndef BIPARTITE_LATTICE_PRINCIPAL_MINOR_CDET_INTEGRAND_SUBLATTICE_RANDOMLY_SAMPLED_H
#define BIPARTITE_LATTICE_PRINCIPAL_MINOR_CDET_INTEGRAND_SUBLATTICE_RANDOMLY_SAMPLED_H

#include "fast_principle_minors.h"
#include <functional>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <exception>
#include <stdexcept>
#include <sstream>
#include "g_funcs.h"
#include <cstring>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <vector>
#include <random>
#include <chrono>
#include <cstdlib>

// Helper function to check vector for NaN values
template<typename T>
inline void check_vector_nan(const std::vector<T>& vec, const std::string& vector_name) {
    for (size_t i = 0; i < vec.size(); ++i) {
        if (is_nan(vec[i])) {
            std::stringstream ss;
            ss << "NaN detected in " << vector_name << "[" << i << "]";
            throw std::runtime_error(ss.str());
        }
    }
}

// Get epsilon threshold for different types
template<typename T>
constexpr T get_epsilon_threshold_1() {
    if constexpr (std::is_same_v<T, float>) {
        return static_cast<T>(1e-6);
    } else if constexpr (std::is_same_v<T, long double>) {
        return static_cast<T>(4e-18L);
    } else if constexpr (std::is_same_v<T, DoubleDouble>) {
        return DoubleDouble(1e-28);
    } else {  // double
        return static_cast<T>(5e-15);
    }
}

// Define a threshold for numerical stability
template<typename T>
constexpr T EPSILON_THRESHOLD = get_epsilon_threshold_1<T>();

/* ================= Z-Transform Implementation ================= */
template<typename T>
T compute_dp_ztransform(const std::vector<T>& up,
                        const std::vector<T>& down,
                        const int order_num) {
    const size_t N = size_t{1} << order_num;

    // aligned alloc for a[] and dp[]
    size_t bytes = ((N * sizeof(T)) + 31) & ~size_t(31);

    auto *a  = (T*)std::aligned_alloc(32, bytes);
    auto *dp = (T*)std::aligned_alloc(32, bytes);

    // build a
    if constexpr (IS_DOUBLEDOUBLE<T>) {
        a[0] = DoubleDouble(0);
    } else {
        a[0] = T(0);
    }
    for (size_t m = 1; m < N; ++m) {
        a[m] = up[m-1] * down[m-1];
    }

    // Initialize dp array properly
    for (size_t i = 0; i < N; ++i) {
        if constexpr (IS_DOUBLEDOUBLE<T>) {
            dp[i] = DoubleDouble(0);
        } else {
            dp[i] = T(0);
        }
    }

    // base
    for (int i = 0; i < order_num; ++i)
        dp[size_t{1} << i] = a[size_t{1} << i];

    // bucket by popcount
    std::vector<std::vector<size_t>> by_cnt(order_num + 1);
    for (size_t m = 0; m < N; ++m)
        by_cnt[popcnt(static_cast<uint32_t>(m))].push_back(m);

    // DP
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

// DoubleDouble version of integrand_temp
inline DoubleDouble integrand_temp_dd(const DoubleDouble beta_s,
                                     const DoubleDouble U_s,
                                     const DoubleDouble* args,
                                     const DD1DVector& g_vec_1_flattened,
                                     const DD1DVector& g_vec_2_flattened,
                                     const int order_num,
                                     const int N1,
                                     const int N2,
                                     const DoubleDouble eps,
                                     const int cheby_degree,
                                     const DoubleDouble alpha_shift_up,
                                     const DoubleDouble alpha_shift_down,
                                     const DoubleDouble Delta_up,
                                     const DoubleDouble Delta_down,
                                     const DoubleDouble tol) {

    // Check bounds
    for (int i = 0, L = 4 * order_num; i < L; ++i) {
        if (args[i] < DoubleDouble(0) || args[i] > DoubleDouble(1))
            return DoubleDouble(0);
    }

    // Build vectors
    DD1DVector vec1(order_num * order_num), vec2(order_num * order_num);
    std::vector<int> sub(order_num);

    for (int i = 0; i < order_num; ++i) {
        DoubleDouble val = DoubleDouble(2) * args[3 * order_num + i];
        sub[i] = DD_FLOOR_TO_INT(val);
    }

    DoubleDouble scale = DD_SQRT(U_s * DoubleDouble(N1) * DoubleDouble(N2) * beta_s);

    for (int i = 0; i < order_num; ++i) {
        for (int j = 0; j < order_num; ++j) {
            int idx = i * order_num + j;
            if (i == j) {
                vec1[idx] = g_func_hartree_flattened_dd(sub[i], sub[j],
                    DoubleDouble(0), DoubleDouble(0), DoubleDouble(0), DoubleDouble(0),
                    g_vec_1_flattened, N1, N2, beta_s, eps, cheby_degree, cheby_degree);
                vec2[idx] = g_func_hartree_flattened_dd(sub[i], sub[j],
                    DoubleDouble(0), DoubleDouble(0), DoubleDouble(0), DoubleDouble(0),
                    g_vec_2_flattened, N1, N2, beta_s, eps, cheby_degree, cheby_degree);

                DoubleDouble shift_down = alpha_shift_down + ((sub[i] == 0) ? Delta_down : -Delta_down);
                DoubleDouble shift_up   = alpha_shift_up   + ((sub[i] == 0) ? Delta_up   : -Delta_up);
                vec1[idx] += shift_down / U_s;
                vec2[idx] += shift_up   / U_s;
            } else {
                DoubleDouble tau1 = beta_s * args[3 * i];
                DoubleDouble tau2 = beta_s * args[3 * j];
                DoubleDouble R1b1 = DoubleDouble(N1) * args[3 * i + 1];
                DoubleDouble R1b2 = DoubleDouble(N2) * args[3 * i + 2];
                DoubleDouble R2b1 = DoubleDouble(N1) * args[3 * j + 1];
                DoubleDouble R2b2 = DoubleDouble(N2) * args[3 * j + 2];

                vec1[idx] = g_func_flattened_dd(sub[i], sub[j],
                                               tau1, tau2, R1b1, R1b2, R2b1, R2b2,
                                               g_vec_1_flattened, beta_s, N1, N2,
                                               eps, cheby_degree, cheby_degree);
                vec2[idx] = g_func_flattened_dd(sub[i], sub[j],
                                               tau1, tau2, R1b1, R1b2, R2b1, R2b2,
                                               g_vec_2_flattened, beta_s, N1, N2,
                                               eps, cheby_degree, cheby_degree);
            }
            vec1[idx] *= scale;
            vec2[idx] *= scale;
        }
    }

    // Check for NaN
    check_vector_nan(vec1, "vec1");
    check_vector_nan(vec2, "vec2");

    // Fast exit for order_num == 1
    if (order_num == 1) {
        DoubleDouble prod  = vec1[0] * vec2[0];
        DoubleDouble extra = U_s * DoubleDouble(N1) * DoubleDouble(N2) * beta_s
                   * ((sub[0] == 0 ? (alpha_shift_down + Delta_down) : (alpha_shift_down - Delta_down)) / U_s)
                   * ((sub[0] == 0 ? (alpha_shift_up   + Delta_up  ) : (alpha_shift_up   - Delta_up  )) / U_s);
        DoubleDouble num  = DD_POW(DoubleDouble(-2), DoubleDouble(1)) * (prod - extra);

        if (!is_finite(num)) {
            throw std::runtime_error("Error (integrand): Value is not finite");
        }
        return num / (DoubleDouble(factorial(1)) * DoubleDouble(2) * DoubleDouble(N1) * DoubleDouble(N2));
    }

    DD1DVector up_mat_minors = mat2pm<DoubleDouble>(vec1, order_num, tol);
    DD1DVector down_mat_minors = mat2pm<DoubleDouble>(vec2, order_num, tol);

    check_vector_nan(up_mat_minors, "up_mat_minors");
    check_vector_nan(down_mat_minors, "down_mat_minors");

    auto c_result = compute_dp_ztransform<DoubleDouble>(up_mat_minors, down_mat_minors, order_num);

    const DoubleDouble val = DD_POW(DoubleDouble(-2), DoubleDouble(order_num)) * c_result /
                            (DoubleDouble(factorial(order_num)) *
                             DoubleDouble(2) * DoubleDouble(N1) * DoubleDouble(N2));

    if (!is_finite(val)) {
        throw std::runtime_error("Error (integrand): Value is not finite");
    }
    return val;
}

/*
// DoubleDouble version of integrand_temp_density
inline DoubleDouble integrand_temp_density_dd(const DoubleDouble beta_s,
                                            const DoubleDouble U_s,
                                            const DoubleDouble* args,
                                            const DD1DVector& g_vec_1_flattened,
                                            const DD1DVector& g_vec_2_flattened,
                                            const int order_num,
                                            const int N1,
                                            const int N2,
                                            const DoubleDouble eps,
                                            const int cheby_degree,
                                            const DoubleDouble alpha_shift_up,
                                            const DoubleDouble alpha_shift_down,
                                            const DoubleDouble Delta_up,
                                            const DoubleDouble Delta_down,
                                            const DoubleDouble tol,
                                            const bool is_spin_up) {

    // Get the base result
    DoubleDouble base_result = integrand_temp_dd(beta_s, U_s, args,
                                               g_vec_1_flattened, g_vec_2_flattened,
                                               order_num, N1, N2, eps, cheby_degree,
                                               alpha_shift_up, alpha_shift_down,
                                               Delta_up, Delta_down, tol);

    // Extract the extra parameter for density calculation
    DoubleDouble s_param = args[4 * order_num];  // The (4*order_num + 1)-th parameter

    // Apply density-specific modification
    DoubleDouble density_factor = is_spin_up ? s_param : (DoubleDouble(1) - s_param);

    return base_result * density_factor;
}

// DoubleDouble version of integrand_temp_double_occupancy
inline DoubleDouble integrand_temp_double_occupancy_dd(const DoubleDouble beta_s,
                                                     const DoubleDouble U_s,
                                                     const DoubleDouble* args,
                                                     const DD1DVector& g_vec_1_flattened,
                                                     const DD1DVector& g_vec_2_flattened,
                                                     const int order_num,
                                                     const int N1,
                                                     const int N2,
                                                     const DoubleDouble eps,
                                                     const int cheby_degree,
                                                     const DoubleDouble alpha_shift_up,
                                                     const DoubleDouble alpha_shift_down,
                                                     const DoubleDouble Delta_up,
                                                     const DoubleDouble Delta_down,
                                                     const DoubleDouble tol) {

    // Get the base result
    DoubleDouble base_result = integrand_temp_dd(beta_s, U_s, args,
                                               g_vec_1_flattened, g_vec_2_flattened,
                                               order_num, N1, N2, eps, cheby_degree,
                                               alpha_shift_up, alpha_shift_down,
                                               Delta_up, Delta_down, tol);

    // Extract the extra parameter for double occupancy calculation
    DoubleDouble s_param = args[4 * order_num];  // The (4*order_num + 1)-th parameter

    // Apply double occupancy-specific modification
    DoubleDouble double_occ_factor = s_param * (DoubleDouble(1) - s_param);

    return base_result * double_occ_factor;
}
*/
#endif