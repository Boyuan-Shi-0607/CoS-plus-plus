#ifndef LOAD_DATA_H
#define LOAD_DATA_H

#include "../utility/resize_and_initialize.h"
#include "parameters.h"
#include "../utility/fft.h"
#include "cheby_decomp.h"
#include "G_0.h"
#include <random>
#include <chrono>
#include <future>
#include <vector>
#include <algorithm>

// Helper function to pre-allocate all map entries to ensure thread safety
inline void preallocate_map_entries() {
    // Pre-allocate entries for lnZ/double_occupancy modes
    if (lnZ_mode) {
        for (int ind = 0; ind < N1_vec.size(); ind++) {
            for (int ind_1 = 0; ind_1 < beta_vec.size(); ind_1++) {
                std::string key = std::to_string(ind) + "_" + std::to_string(ind_1);

                // Pre-allocate and resize g_vec entries - now 3D for scalar
                resizeMultiDimVector(g_vec[key], {cheby_degree+1, N1_vec[ind], N2_vec[ind]});

                // Pre-allocate flattened entries (empty vectors)
                g_vec_flattened[key] = Real1DVector();
            }
        }
    }

    // Pre-allocate entries for density mode
    if (density_mode) {
        for (int ind = 0; ind < N1_vec.size(); ind++) {
            for (int ind_1 = 0; ind_1 < beta_vec.size(); ind_1++) {
                std::string key = std::to_string(ind) + "_" + std::to_string(ind_1);

                // Pre-allocate and resize g_vec_mu_h entries - now 3D for scalar
                resizeMultiDimVector(g_vec_mu_h[key], {cheby_degree+1, N1_vec[ind], N2_vec[ind]});

                // Pre-allocate and resize g_vec_mu_l entries - now 3D for scalar
                resizeMultiDimVector(g_vec_mu_l[key], {cheby_degree+1, N1_vec[ind], N2_vec[ind]});

                // Pre-allocate flattened entries
                g_vec_mu_h_flattened[key] = Real1DVector();
                g_vec_mu_l_flattened[key] = Real1DVector();
            }
        }
    }
}

// Helper function to process a single (ind, ind_1) pair for lnZ/double_occupancy modes
inline void process_lnZ_pair(int ind, int ind_1) {
    Complex3DVector g_k_space;
    resizeMultiDimVector(g_k_space, {N1_vec[ind], N2_vec[ind], cheby_degree+1});

    // Calculate Green's functions for each k-point
    for (int n_k1 = 0; n_k1 < N1_vec[ind]; n_k1++) {
        for (int n_k2 = 0; n_k2 < N2_vec[ind]; n_k2++) {
            // Both spin up and down use the same chemical potential
            auto partialG0 = [=](Real tau) -> Real {
                return G_0(beta_vec[ind_1], mu-alpha_shift,
                          N1_vec[ind], N2_vec[ind], n_k1, n_k2, greens_func_exp_cut, tau);
            };
            Real1DVector g_coeffs = chebyshev_decomposition_scalar(partialG0,
                                                                   eps, beta_vec[ind_1]-eps, cheby_degree);

            // Store coefficients as complex numbers (with zero imaginary part)
            for (int k = 0; k <= cheby_degree; k++) {
                g_k_space[n_k1][n_k2][k] = Complex(g_coeffs[k], 0.0);
            }
        }
    }

    // Transform from k-space to real space
    std::string key = std::to_string(ind) + "_" + std::to_string(ind_1);
    transform_vec(g_k_space, g_vec[key]);

    // Flatten vectors
    flatten_3D_to_1D(g_vec[key], g_vec_flattened[key]);
}

// Helper functions to process mu_h and mu_l separately for density mode
inline void process_density_mu_h(int ind, int ind_1) {
    Complex3DVector g_mu_h_k_space;
    resizeMultiDimVector(g_mu_h_k_space, {N1_vec[ind], N2_vec[ind], cheby_degree+1});

    for (int n_k1 = 0; n_k1 < N1_vec[ind]; n_k1++) {
        for (int n_k2 = 0; n_k2 < N2_vec[ind]; n_k2++) {
            auto partialG0_mu_h = [=](Real tau) -> Real {
                return G_0(beta_vec[ind_1], mu_h-alpha_shift,
                          N1_vec[ind], N2_vec[ind], n_k1, n_k2, greens_func_exp_cut, tau);
            };
            Real1DVector g_mu_h_coeffs = chebyshev_decomposition_scalar(partialG0_mu_h,
                                                                        eps, beta_vec[ind_1]-eps, cheby_degree);

            for (int k = 0; k <= cheby_degree; k++) {
                g_mu_h_k_space[n_k1][n_k2][k] = Complex(g_mu_h_coeffs[k], 0.0);
            }
        }
    }

    std::string key = std::to_string(ind) + "_" + std::to_string(ind_1);
    transform_vec(g_mu_h_k_space, g_vec_mu_h[key]);
    flatten_3D_to_1D(g_vec_mu_h[key], g_vec_mu_h_flattened[key]);
}

inline void process_density_mu_l(int ind, int ind_1) {
    Complex3DVector g_mu_l_k_space;
    resizeMultiDimVector(g_mu_l_k_space, {N1_vec[ind], N2_vec[ind], cheby_degree+1});

    for (int n_k1 = 0; n_k1 < N1_vec[ind]; n_k1++) {
        for (int n_k2 = 0; n_k2 < N2_vec[ind]; n_k2++) {
            auto partialG0_mu_l = [=](Real tau) -> Real {
                return G_0(beta_vec[ind_1], mu_l-alpha_shift,
                          N1_vec[ind], N2_vec[ind], n_k1, n_k2, greens_func_exp_cut, tau);
            };
            Real1DVector g_mu_l_coeffs = chebyshev_decomposition_scalar(partialG0_mu_l,
                                                                        eps, beta_vec[ind_1]-eps, cheby_degree);

            for (int k = 0; k <= cheby_degree; k++) {
                g_mu_l_k_space[n_k1][n_k2][k] = Complex(g_mu_l_coeffs[k], 0.0);
            }
        }
    }

    std::string key = std::to_string(ind) + "_" + std::to_string(ind_1);
    transform_vec(g_mu_l_k_space, g_vec_mu_l[key]);
    flatten_3D_to_1D(g_vec_mu_l[key], g_vec_mu_l_flattened[key]);
}

// Generic parallel executor for nested loops
template<typename ProcessFunc>
void parallel_nested_loop(size_t n1_size, size_t beta_size, ProcessFunc process_func) {
    size_t total_tasks = n1_size * beta_size;
    size_t num_threads = std::min(static_cast<size_t>(cores_cpp), total_tasks);

    std::vector<std::future<void>> futures;
    futures.reserve(total_tasks);

    // Launch all tasks
    for (size_t ind = 0; ind < n1_size; ++ind) {
        for (size_t ind_1 = 0; ind_1 < beta_size; ++ind_1) {
            futures.push_back(std::async(std::launch::async, process_func, ind, ind_1));
        }
    }

    // Wait for all tasks to complete
    for (auto& future : futures) {
        future.wait();
    }
}

// Enhanced parallel executor for nested loops with h/l variants
template<typename ProcessFuncH, typename ProcessFuncL>
void parallel_nested_loop_with_variants(size_t n1_size, size_t beta_size,
                                      ProcessFuncH process_func_h, ProcessFuncL process_func_l) {
    size_t total_tasks = n1_size * beta_size * 2; // *2 for h and l variants
    size_t num_threads = std::min(static_cast<size_t>(cores_cpp), total_tasks);

    std::vector<std::future<void>> futures;
    futures.reserve(total_tasks);

    // Launch all tasks - both h and l variants
    for (size_t ind = 0; ind < n1_size; ++ind) {
        for (size_t ind_1 = 0; ind_1 < beta_size; ++ind_1) {
            // Launch h variant
            futures.push_back(std::async(std::launch::async, process_func_h, ind, ind_1));
            // Launch l variant
            futures.push_back(std::async(std::launch::async, process_func_l, ind, ind_1));
        }
    }

    // Wait for all tasks to complete
    for (auto& future : futures) {
        future.wait();
    }
}

inline void load_data() {
    // CRITICAL: Pre-allocate all map entries before parallel processing
    preallocate_map_entries();

    // Load data for lnZ and double occupancy modes
    if (lnZ_mode) {
        // Parallel calculation of Green's functions
        parallel_nested_loop(N1_vec.size(), beta_vec.size(), process_lnZ_pair);
    }

    // Load data for density mode with parallelized h/l variants
    if (density_mode) {
        // Parallel calculation of Green's functions for different chemical potentials
        parallel_nested_loop_with_variants(N1_vec.size(), beta_vec.size(),
                                         process_density_mu_h, process_density_mu_l);
    }
}

#endif