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

// Helper function to flatten a 5D MP_Real vector to a 1D MP_Real vector
inline void flatten_5D_to_1D_mp(const MP_Real5DVector& vec5D, MP_Real1DVector& vec1D) {
    vec1D.clear();
    for (size_t i = 0; i < vec5D.size(); ++i) {
        for (size_t j = 0; j < vec5D[i].size(); ++j) {
            for (size_t k = 0; k < vec5D[i][j].size(); ++k) {
                for (size_t l = 0; l < vec5D[i][j][k].size(); ++l) {
                    for (size_t m = 0; m < vec5D[i][j][k][l].size(); ++m) {
                        vec1D.push_back(vec5D[i][j][k][l][m]);
                    }
                }
            }
        }
    }
}

// Helper function to determine if we need data for a specific (ind, ind_1) pair
inline bool need_data_for_indices(int ind, int ind_1, size_t n1_size, size_t beta_size) {
    // We need data for:
    // 1. All grid sizes at beta index 0: (i, 0) for i in [0, n1_size-1]
    // 2. Final grid size at all beta indices: (n1_size-1, j) for j in [0, beta_size-1]
    // This forms an L-shaped pattern
    return (ind_1 == 0) || (ind == static_cast<int>(n1_size - 1));
}

// Helper function to pre-allocate only necessary map entries
inline void preallocate_map_entries() {
    size_t n1_size = N1_vec.size();
    size_t beta_size = beta_vec.size();

    // Pre-allocate entries for lnZ/double_occupancy modes
    if (lnZ_mode || double_occupancy_mode || density_direct_mode || double_occupancy_direct_mode) {
        for (int ind = 0; ind < n1_size; ind++) {
            for (int ind_1 = 0; ind_1 < beta_size; ind_1++) {
                if (!need_data_for_indices(ind, ind_1, n1_size, beta_size)) continue;

                std::string key_up = "up_" + std::to_string(ind) + "_" + std::to_string(ind_1);
                std::string key_down = "down_" + std::to_string(ind) + "_" + std::to_string(ind_1);

                // Pre-allocate and resize g_vec entries
                resizeMultiDimVector(g_vec[key_up], {cheby_degree+1, 2, 2, N1_vec[ind], N2_vec[ind]});
                resizeMultiDimVector(g_vec[key_down], {cheby_degree+1, 2, 2, N1_vec[ind], N2_vec[ind]});

                // Pre-allocate flattened entries (empty vectors)
                g_vec_flattened_mp[key_up] = MP_Real1DVector();
                g_vec_flattened_mp[key_down] = MP_Real1DVector();
            }
        }
    }

    // Pre-allocate entries for density mode
    if (density_mode) {
        for (int ind = 0; ind < n1_size; ind++) {
            for (int ind_1 = 0; ind_1 < beta_size; ind_1++) {
                if (!need_data_for_indices(ind, ind_1, n1_size, beta_size)) continue;

                std::string key_up = "up_" + std::to_string(ind) + "_" + std::to_string(ind_1);
                std::string key_down = "down_" + std::to_string(ind) + "_" + std::to_string(ind_1);

                // Pre-allocate and resize g_vec_mu_h entries
                resizeMultiDimVector(g_vec_mu_h[key_up], {cheby_degree+1, 2, 2, N1_vec[ind], N2_vec[ind]});
                resizeMultiDimVector(g_vec_mu_h[key_down], {cheby_degree+1, 2, 2, N1_vec[ind], N2_vec[ind]});

                // Pre-allocate and resize g_vec_mu_l entries
                resizeMultiDimVector(g_vec_mu_l[key_up], {cheby_degree+1, 2, 2, N1_vec[ind], N2_vec[ind]});
                resizeMultiDimVector(g_vec_mu_l[key_down], {cheby_degree+1, 2, 2, N1_vec[ind], N2_vec[ind]});

                // Pre-allocate flattened entries
                g_vec_mu_h_flattened_mp[key_up] = MP_Real1DVector();
                g_vec_mu_h_flattened_mp[key_down] = MP_Real1DVector();
                g_vec_mu_l_flattened_mp[key_up] = MP_Real1DVector();
                g_vec_mu_l_flattened_mp[key_down] = MP_Real1DVector();
            }
        }
    }

    // Pre-allocate entries for energy mode
    if (energy_mode) {
        for (int ind = 0; ind < n1_size; ind++) {
            for (int ind_1 = 0; ind_1 < beta_size; ind_1++) {
                if (!need_data_for_indices(ind, ind_1, n1_size, beta_size)) continue;

                std::string key_up = "up_" + std::to_string(ind) + "_" + std::to_string(ind_1);
                std::string key_down = "down_" + std::to_string(ind) + "_" + std::to_string(ind_1);

                // Pre-allocate and resize g_vec_beta_h entries
                resizeMultiDimVector(g_vec_beta_h[key_up], {cheby_degree+1, 2, 2, N1_vec[ind], N2_vec[ind]});
                resizeMultiDimVector(g_vec_beta_h[key_down], {cheby_degree+1, 2, 2, N1_vec[ind], N2_vec[ind]});

                // Pre-allocate and resize g_vec_beta_l entries
                resizeMultiDimVector(g_vec_beta_l[key_up], {cheby_degree+1, 2, 2, N1_vec[ind], N2_vec[ind]});
                resizeMultiDimVector(g_vec_beta_l[key_down], {cheby_degree+1, 2, 2, N1_vec[ind], N2_vec[ind]});

                // Pre-allocate flattened entries
                g_vec_beta_h_flattened_mp[key_up] = MP_Real1DVector();
                g_vec_beta_h_flattened_mp[key_down] = MP_Real1DVector();
                g_vec_beta_l_flattened_mp[key_up] = MP_Real1DVector();
                g_vec_beta_l_flattened_mp[key_down] = MP_Real1DVector();
            }
        }
    }

    // Pre-allocate entries for staggered magnetization mode
    if (staggered_magnetization_mode) {
        for (int ind = 0; ind < n1_size; ind++) {
            for (int ind_1 = 0; ind_1 < beta_size; ind_1++) {
                if (!need_data_for_indices(ind, ind_1, n1_size, beta_size)) continue;

                std::string key_up = "up_" + std::to_string(ind) + "_" + std::to_string(ind_1);
                std::string key_down = "down_" + std::to_string(ind) + "_" + std::to_string(ind_1);

                // Pre-allocate and resize g_vec_H_h entries
                resizeMultiDimVector(g_vec_H_h[key_up], {cheby_degree+1, 2, 2, N1_vec[ind], N2_vec[ind]});
                resizeMultiDimVector(g_vec_H_h[key_down], {cheby_degree+1, 2, 2, N1_vec[ind], N2_vec[ind]});

                // Pre-allocate and resize g_vec_H_l entries
                resizeMultiDimVector(g_vec_H_l[key_up], {cheby_degree+1, 2, 2, N1_vec[ind], N2_vec[ind]});
                resizeMultiDimVector(g_vec_H_l[key_down], {cheby_degree+1, 2, 2, N1_vec[ind], N2_vec[ind]});

                // Pre-allocate flattened entries
                g_vec_H_h_flattened_mp[key_up] = MP_Real1DVector();
                g_vec_H_h_flattened_mp[key_down] = MP_Real1DVector();
                g_vec_H_l_flattened_mp[key_up] = MP_Real1DVector();
                g_vec_H_l_flattened_mp[key_down] = MP_Real1DVector();
            }
        }
    }
}

// Helper function to process a single (ind, ind_1) pair for lnZ/double_occupancy modes
inline void process_lnZ_pair(int ind, int ind_1) {
    MP_Complex5DVector g_up_k_space;
    MP_Complex5DVector g_down_k_space;
    resizeMultiDimVector(g_up_k_space, {N1_vec[ind], N2_vec[ind], cheby_degree+1, 2, 2});
    resizeMultiDimVector(g_down_k_space, {N1_vec[ind], N2_vec[ind], cheby_degree+1, 2, 2});

    // Calculate Green's functions for each k-point
    for (int n_k1 = 0; n_k1 < N1_vec[ind]; n_k1++) {
        for (int n_k2 = 0; n_k2 < N2_vec[ind]; n_k2++) {
            // Up spin
            auto partialG0_up = [=](MP_Real tau) -> MP_Complex2DVector {
                return G_0_mp(beta_vec[ind_1], MP_Real(0), Delta_up, mu-alpha_shift_up,
                            N1_vec[ind], N2_vec[ind], n_k1, n_k2, greens_func_exp_cut, tau);
            };
            g_up_k_space[n_k1][n_k2] = chebyshev_decomposition_matrix_mp(partialG0_up,
                                                                        eps, beta_vec[ind_1]-eps, cheby_degree);

            // Down spin
            auto partialG0_down = [=](MP_Real tau) -> MP_Complex2DVector {
                return G_0_mp(beta_vec[ind_1], MP_Real(0), Delta_down, mu-alpha_shift_down,
                            N1_vec[ind], N2_vec[ind], n_k1, n_k2, greens_func_exp_cut, tau);
            };
            g_down_k_space[n_k1][n_k2] = chebyshev_decomposition_matrix_mp(partialG0_down,
                                                                          eps, beta_vec[ind_1]-eps, cheby_degree);
        }
    }

    // Transform from k-space to real space
    // Since we pre-allocated the map entries, these operations are now thread-safe
    transform_vec_mp(g_up_k_space, g_vec["up_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);
    transform_vec_mp(g_down_k_space, g_vec["down_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);

    // Flatten vectors to MP_Real version
    flatten_5D_to_1D_mp(g_vec["up_" + std::to_string(ind) + "_" + std::to_string(ind_1)],
                       g_vec_flattened_mp["up_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);
    flatten_5D_to_1D_mp(g_vec["down_" + std::to_string(ind) + "_" + std::to_string(ind_1)],
                       g_vec_flattened_mp["down_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);
}

// Helper functions to process mu_h and mu_l separately for density mode
inline void process_density_mu_h(int ind, int ind_1) {
    MP_Complex5DVector g_mu_h_up_k_space;
    MP_Complex5DVector g_mu_h_down_k_space;
    resizeMultiDimVector(g_mu_h_up_k_space, {N1_vec[ind], N2_vec[ind], cheby_degree+1, 2, 2});
    resizeMultiDimVector(g_mu_h_down_k_space, {N1_vec[ind], N2_vec[ind], cheby_degree+1, 2, 2});

    for (int n_k1 = 0; n_k1 < N1_vec[ind]; n_k1++) {
        for (int n_k2 = 0; n_k2 < N2_vec[ind]; n_k2++) {
            auto partialG0_mu_h_up = [=](MP_Real tau) -> MP_Complex2DVector {
                return G_0_mp(beta_vec[ind_1], MP_Real(0), Delta_up, mu_h-alpha_shift_up,
                            N1_vec[ind], N2_vec[ind], n_k1, n_k2, greens_func_exp_cut, tau);
            };
            g_mu_h_up_k_space[n_k1][n_k2] = chebyshev_decomposition_matrix_mp(partialG0_mu_h_up,
                                                                             eps, beta_vec[ind_1]-eps, cheby_degree);

            auto partialG0_mu_h_down = [=](MP_Real tau) -> MP_Complex2DVector {
                return G_0_mp(beta_vec[ind_1], MP_Real(0), Delta_down, mu_h-alpha_shift_down,
                            N1_vec[ind], N2_vec[ind], n_k1, n_k2, greens_func_exp_cut, tau);
            };
            g_mu_h_down_k_space[n_k1][n_k2] = chebyshev_decomposition_matrix_mp(partialG0_mu_h_down,
                                                                               eps, beta_vec[ind_1]-eps, cheby_degree);
        }
    }

    transform_vec_mp(g_mu_h_up_k_space, g_vec_mu_h["up_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);
    transform_vec_mp(g_mu_h_down_k_space, g_vec_mu_h["down_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);

    flatten_5D_to_1D_mp(g_vec_mu_h["up_" + std::to_string(ind) + "_" + std::to_string(ind_1)],
                       g_vec_mu_h_flattened_mp["up_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);
    flatten_5D_to_1D_mp(g_vec_mu_h["down_" + std::to_string(ind) + "_" + std::to_string(ind_1)],
                       g_vec_mu_h_flattened_mp["down_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);
}

inline void process_density_mu_l(int ind, int ind_1) {
    MP_Complex5DVector g_mu_l_up_k_space;
    MP_Complex5DVector g_mu_l_down_k_space;
    resizeMultiDimVector(g_mu_l_up_k_space, {N1_vec[ind], N2_vec[ind], cheby_degree+1, 2, 2});
    resizeMultiDimVector(g_mu_l_down_k_space, {N1_vec[ind], N2_vec[ind], cheby_degree+1, 2, 2});

    for (int n_k1 = 0; n_k1 < N1_vec[ind]; n_k1++) {
        for (int n_k2 = 0; n_k2 < N2_vec[ind]; n_k2++) {
            auto partialG0_mu_l_up = [=](MP_Real tau) -> MP_Complex2DVector {
                return G_0_mp(beta_vec[ind_1], MP_Real(0), Delta_up, mu_l-alpha_shift_up,
                            N1_vec[ind], N2_vec[ind], n_k1, n_k2, greens_func_exp_cut, tau);
            };
            g_mu_l_up_k_space[n_k1][n_k2] = chebyshev_decomposition_matrix_mp(partialG0_mu_l_up,
                                                                             eps, beta_vec[ind_1]-eps, cheby_degree);

            auto partialG0_mu_l_down = [=](MP_Real tau) -> MP_Complex2DVector {
                return G_0_mp(beta_vec[ind_1], MP_Real(0), Delta_down, mu_l-alpha_shift_down,
                            N1_vec[ind], N2_vec[ind], n_k1, n_k2, greens_func_exp_cut, tau);
            };
            g_mu_l_down_k_space[n_k1][n_k2] = chebyshev_decomposition_matrix_mp(partialG0_mu_l_down,
                                                                               eps, beta_vec[ind_1]-eps, cheby_degree);
        }
    }

    transform_vec_mp(g_mu_l_up_k_space, g_vec_mu_l["up_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);
    transform_vec_mp(g_mu_l_down_k_space, g_vec_mu_l["down_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);

    flatten_5D_to_1D_mp(g_vec_mu_l["up_" + std::to_string(ind) + "_" + std::to_string(ind_1)],
                       g_vec_mu_l_flattened_mp["up_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);
    flatten_5D_to_1D_mp(g_vec_mu_l["down_" + std::to_string(ind) + "_" + std::to_string(ind_1)],
                       g_vec_mu_l_flattened_mp["down_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);
}

// Helper functions to process beta_h and beta_l separately for energy mode
inline void process_energy_beta_h(int ind, int ind_1) {
    MP_Complex5DVector g_beta_h_up_k_space;
    MP_Complex5DVector g_beta_h_down_k_space;
    resizeMultiDimVector(g_beta_h_up_k_space, {N1_vec[ind], N2_vec[ind], cheby_degree+1, 2, 2});
    resizeMultiDimVector(g_beta_h_down_k_space, {N1_vec[ind], N2_vec[ind], cheby_degree+1, 2, 2});

    for (int n_k1 = 0; n_k1 < N1_vec[ind]; n_k1++) {
        for (int n_k2 = 0; n_k2 < N2_vec[ind]; n_k2++) {
            auto partialG0_beta_h_up = [=](MP_Real tau) -> MP_Complex2DVector {
                return G_0_mp(beta_h_vec[ind_1], MP_Real(0), Delta_up, mu-alpha_shift_up,
                            N1_vec[ind], N2_vec[ind], n_k1, n_k2, greens_func_exp_cut, tau);
            };
            g_beta_h_up_k_space[n_k1][n_k2] = chebyshev_decomposition_matrix_mp(partialG0_beta_h_up,
                                                                               eps, beta_h_vec[ind_1]-eps, cheby_degree);

            auto partialG0_beta_h_down = [=](MP_Real tau) -> MP_Complex2DVector {
                return G_0_mp(beta_h_vec[ind_1], MP_Real(0), Delta_down, mu-alpha_shift_down,
                            N1_vec[ind], N2_vec[ind], n_k1, n_k2, greens_func_exp_cut, tau);
            };
            g_beta_h_down_k_space[n_k1][n_k2] = chebyshev_decomposition_matrix_mp(partialG0_beta_h_down,
                                                                                 eps, beta_h_vec[ind_1]-eps, cheby_degree);
        }
    }

    transform_vec_mp(g_beta_h_up_k_space, g_vec_beta_h["up_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);
    transform_vec_mp(g_beta_h_down_k_space, g_vec_beta_h["down_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);

    flatten_5D_to_1D_mp(g_vec_beta_h["up_" + std::to_string(ind) + "_" + std::to_string(ind_1)],
                       g_vec_beta_h_flattened_mp["up_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);
    flatten_5D_to_1D_mp(g_vec_beta_h["down_" + std::to_string(ind) + "_" + std::to_string(ind_1)],
                       g_vec_beta_h_flattened_mp["down_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);
}

inline void process_energy_beta_l(int ind, int ind_1) {
    MP_Complex5DVector g_beta_l_up_k_space;
    MP_Complex5DVector g_beta_l_down_k_space;
    resizeMultiDimVector(g_beta_l_up_k_space, {N1_vec[ind], N2_vec[ind], cheby_degree+1, 2, 2});
    resizeMultiDimVector(g_beta_l_down_k_space, {N1_vec[ind], N2_vec[ind], cheby_degree+1, 2, 2});

    for (int n_k1 = 0; n_k1 < N1_vec[ind]; n_k1++) {
        for (int n_k2 = 0; n_k2 < N2_vec[ind]; n_k2++) {
            auto partialG0_beta_l_up = [=](MP_Real tau) -> MP_Complex2DVector {
                return G_0_mp(beta_l_vec[ind_1], MP_Real(0), Delta_up, mu-alpha_shift_up,
                            N1_vec[ind], N2_vec[ind], n_k1, n_k2, greens_func_exp_cut, tau);
            };
            g_beta_l_up_k_space[n_k1][n_k2] = chebyshev_decomposition_matrix_mp(partialG0_beta_l_up,
                                                                               eps, beta_l_vec[ind_1]-eps, cheby_degree);

            auto partialG0_beta_l_down = [=](MP_Real tau) -> MP_Complex2DVector {
                return G_0_mp(beta_l_vec[ind_1], MP_Real(0), Delta_down, mu-alpha_shift_down,
                            N1_vec[ind], N2_vec[ind], n_k1, n_k2, greens_func_exp_cut, tau);
            };
            g_beta_l_down_k_space[n_k1][n_k2] = chebyshev_decomposition_matrix_mp(partialG0_beta_l_down,
                                                                                 eps, beta_l_vec[ind_1]-eps, cheby_degree);
        }
    }

    transform_vec_mp(g_beta_l_up_k_space, g_vec_beta_l["up_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);
    transform_vec_mp(g_beta_l_down_k_space, g_vec_beta_l["down_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);

    flatten_5D_to_1D_mp(g_vec_beta_l["up_" + std::to_string(ind) + "_" + std::to_string(ind_1)],
                       g_vec_beta_l_flattened_mp["up_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);
    flatten_5D_to_1D_mp(g_vec_beta_l["down_" + std::to_string(ind) + "_" + std::to_string(ind_1)],
                       g_vec_beta_l_flattened_mp["down_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);
}

// Helper functions to process H_h and H_l separately for staggered magnetization mode
inline void process_staggered_magnetization_H_h(int ind, int ind_1) {
    MP_Complex5DVector g_H_h_up_k_space;
    MP_Complex5DVector g_H_h_down_k_space;
    resizeMultiDimVector(g_H_h_up_k_space, {N1_vec[ind], N2_vec[ind], cheby_degree+1, 2, 2});
    resizeMultiDimVector(g_H_h_down_k_space, {N1_vec[ind], N2_vec[ind], cheby_degree+1, 2, 2});

    for (int n_k1 = 0; n_k1 < N1_vec[ind]; n_k1++) {
        for (int n_k2 = 0; n_k2 < N2_vec[ind]; n_k2++) {
            auto partialG0_H_h_up = [=](MP_Real tau) -> MP_Complex2DVector {
                return G_0_mp(beta_vec[ind_1], H_h, Delta_up, mu-alpha_shift_up,
                            N1_vec[ind], N2_vec[ind], n_k1, n_k2, greens_func_exp_cut, tau);
            };
            g_H_h_up_k_space[n_k1][n_k2] = chebyshev_decomposition_matrix_mp(partialG0_H_h_up,
                                                                            eps, beta_vec[ind_1]-eps, cheby_degree);

            auto partialG0_H_h_down = [=](MP_Real tau) -> MP_Complex2DVector {
                return G_0_mp(beta_vec[ind_1], -H_h, Delta_down, mu-alpha_shift_down,
                            N1_vec[ind], N2_vec[ind], n_k1, n_k2, greens_func_exp_cut, tau);
            };
            g_H_h_down_k_space[n_k1][n_k2] = chebyshev_decomposition_matrix_mp(partialG0_H_h_down,
                                                                              eps, beta_vec[ind_1]-eps, cheby_degree);
        }
    }

    transform_vec_mp(g_H_h_up_k_space, g_vec_H_h["up_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);
    transform_vec_mp(g_H_h_down_k_space, g_vec_H_h["down_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);

    flatten_5D_to_1D_mp(g_vec_H_h["up_" + std::to_string(ind) + "_" + std::to_string(ind_1)],
                       g_vec_H_h_flattened_mp["up_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);
    flatten_5D_to_1D_mp(g_vec_H_h["down_" + std::to_string(ind) + "_" + std::to_string(ind_1)],
                       g_vec_H_h_flattened_mp["down_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);
}

inline void process_staggered_magnetization_H_l(int ind, int ind_1) {
    MP_Complex5DVector g_H_l_up_k_space;
    MP_Complex5DVector g_H_l_down_k_space;
    resizeMultiDimVector(g_H_l_up_k_space, {N1_vec[ind], N2_vec[ind], cheby_degree+1, 2, 2});
    resizeMultiDimVector(g_H_l_down_k_space, {N1_vec[ind], N2_vec[ind], cheby_degree+1, 2, 2});

    for (int n_k1 = 0; n_k1 < N1_vec[ind]; n_k1++) {
        for (int n_k2 = 0; n_k2 < N2_vec[ind]; n_k2++) {
            auto partialG0_H_l_up = [=](MP_Real tau) -> MP_Complex2DVector {
                return G_0_mp(beta_vec[ind_1], H_l, Delta_up, mu-alpha_shift_up,
                            N1_vec[ind], N2_vec[ind], n_k1, n_k2, greens_func_exp_cut, tau);
            };
            g_H_l_up_k_space[n_k1][n_k2] = chebyshev_decomposition_matrix_mp(partialG0_H_l_up,
                                                                            eps, beta_vec[ind_1]-eps, cheby_degree);

            auto partialG0_H_l_down = [=](MP_Real tau) -> MP_Complex2DVector {
                return G_0_mp(beta_vec[ind_1], -H_l, Delta_down, mu-alpha_shift_down,
                            N1_vec[ind], N2_vec[ind], n_k1, n_k2, greens_func_exp_cut, tau);
            };
            g_H_l_down_k_space[n_k1][n_k2] = chebyshev_decomposition_matrix_mp(partialG0_H_l_down,
                                                                              eps, beta_vec[ind_1]-eps, cheby_degree);
        }
    }

    transform_vec_mp(g_H_l_up_k_space, g_vec_H_l["up_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);
    transform_vec_mp(g_H_l_down_k_space, g_vec_H_l["down_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);

    flatten_5D_to_1D_mp(g_vec_H_l["up_" + std::to_string(ind) + "_" + std::to_string(ind_1)],
                       g_vec_H_l_flattened_mp["up_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);
    flatten_5D_to_1D_mp(g_vec_H_l["down_" + std::to_string(ind) + "_" + std::to_string(ind_1)],
                       g_vec_H_l_flattened_mp["down_" + std::to_string(ind) + "_" + std::to_string(ind_1)]);
}

// Optimized parallel executor that only processes needed indices
template<typename ProcessFunc>
void parallel_selective_loop(size_t n1_size, size_t beta_size, ProcessFunc process_func) {
    std::vector<std::future<void>> futures;

    // Only process L-shaped pattern:
    // 1. Row 0: all grid sizes at beta index 0
    for (size_t ind = 0; ind < n1_size; ++ind) {
        futures.push_back(std::async(std::launch::async, process_func, ind, 0));
    }

    // 2. Last column: final grid size at all beta indices (skip (n1_size-1, 0) to avoid duplicate)
    for (size_t ind_1 = 1; ind_1 < beta_size; ++ind_1) {
        futures.push_back(std::async(std::launch::async, process_func, n1_size - 1, ind_1));
    }

    // Wait for all tasks to complete
    for (auto& future : futures) {
        future.wait();
    }
}

// Optimized parallel executor for h/l variants
template<typename ProcessFuncH, typename ProcessFuncL>
void parallel_selective_loop_with_variants(size_t n1_size, size_t beta_size,
                                         ProcessFuncH process_func_h, ProcessFuncL process_func_l) {
    std::vector<std::future<void>> futures;

    // Only process L-shaped pattern:
    // 1. Row 0: all grid sizes at beta index 0
    for (size_t ind = 0; ind < n1_size; ++ind) {
        futures.push_back(std::async(std::launch::async, process_func_h, ind, 0));
        futures.push_back(std::async(std::launch::async, process_func_l, ind, 0));
    }

    // 2. Last column: final grid size at all beta indices (skip (n1_size-1, 0) to avoid duplicate)
    for (size_t ind_1 = 1; ind_1 < beta_size; ++ind_1) {
        futures.push_back(std::async(std::launch::async, process_func_h, n1_size - 1, ind_1));
        futures.push_back(std::async(std::launch::async, process_func_l, n1_size - 1, ind_1));
    }

    // Wait for all tasks to complete
    for (auto& future : futures) {
        future.wait();
    }
}

inline void load_data() {
    // CRITICAL: Pre-allocate all map entries before parallel processing
    preallocate_map_entries();

    size_t n1_size = N1_vec.size();
    size_t beta_size = beta_vec.size();

    std::cout << "Loading data in L-shaped pattern:\n";
    std::cout << "  - Row 0: " << n1_size << " entries (grid indices 0-" << n1_size-1 << " at beta index 0)\n";
    std::cout << "  - Column " << n1_size-1 << ": " << beta_size-1 << " entries (grid index " << n1_size-1 << " at beta indices 1-" << beta_size-1 << ")\n";
    std::cout << "  - Total entries: " << (n1_size + beta_size - 1) << " (vs " << n1_size * beta_size << " for full grid)\n";

    // Load data for lnZ and double occupancy modes
    if (lnZ_mode || double_occupancy_mode || density_direct_mode || double_occupancy_direct_mode) {
        // Parallel calculation of Green's functions - only for needed indices
        parallel_selective_loop(n1_size, beta_size, process_lnZ_pair);
    }

    // Load data for density mode with parallelized h/l variants
    if (density_mode) {
        // Parallel calculation of Green's functions for different chemical potentials
        parallel_selective_loop_with_variants(n1_size, beta_size,
                                            process_density_mu_h, process_density_mu_l);
    }

    // Load data for energy mode with parallelized h/l variants
    if (energy_mode) {
        // Parallel calculation of Green's functions for different beta values
        parallel_selective_loop_with_variants(n1_size, beta_size,
                                            process_energy_beta_h, process_energy_beta_l);
    }

    // Load data for staggered magnetization mode with parallelized h/l variants
    if (staggered_magnetization_mode) {
        // Parallel calculation of Green's functions for different H values
        parallel_selective_loop_with_variants(n1_size, beta_size,
                                            process_staggered_magnetization_H_h, process_staggered_magnetization_H_l);
    }

    // Compressibility mode uses the same data as density mode
    if (compressibility_mode) {
        // Uses g_vec_mu_h and g_vec_mu_l from density mode
        if (!density_mode && !density_direct_mode) {
            std::cerr << "Warning: Compressibility mode requires density mode data\n";
        }
    }
}

#endif