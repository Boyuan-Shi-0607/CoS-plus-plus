#ifndef LOAD_DATA_H
#define LOAD_DATA_H

#include "../utility/resize_and_initialize.h"
#include "parameters.h"
#include "../utility/fft.h"
#include "cheby_decomp.h"
#include "G_0.h"

#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>

// -----------------------------------------------------------------------------
// Pre-allocation for all (ind, ind_1) keys and all flavours
// -----------------------------------------------------------------------------
inline void preallocate_multi_flavour_maps() {
    if (N_f <= 0) throw std::runtime_error("N_f must be > 0");
    if (cheby_degree < 0) throw std::runtime_error("cheby_degree must be >= 0");
    if (N1_vec.size() != N2_vec.size())
        throw std::runtime_error("N1_vec and N2_vec must have the same size");
    if (alpha_shift_vec.size() != static_cast<size_t>(N_f))
        throw std::runtime_error("alpha_shift_vec size must equal N_f");

    const int n_sizes = static_cast<int>(N1_vec.size());
    const int n_betas = static_cast<int>(beta_vec.size());

    for (int ind = 0; ind < n_sizes; ++ind) {
        const int N1i = N1_vec[ind];
        const int N2i = N2_vec[ind];
        if (N1i <= 0 || N2i <= 0) throw std::runtime_error("N1_vec/N2_vec entries must be > 0");

        const int flattened_len = (cheby_degree + 1) * N1i * N2i;

        for (int ind_1 = 0; ind_1 < n_betas; ++ind_1) {
            const std::string key = std::to_string(ind) + "_" + std::to_string(ind_1);

            // lnZ / double-occupancy container
            if (lnZ_mode) {
                // [flavour][cheby+1][N1][N2]
                resizeMultiDimVector(g_vec[key], { N_f, (cheby_degree + 1), N1i, N2i });
                // [flavour][(cheby+1)*N1*N2]
                resizeMultiDimVector(g_vec_flattened[key], { N_f, flattened_len });
            }

            if (density_mode) {
                resizeMultiDimVector(g_vec_mu_h[key], { N_f, (cheby_degree + 1), N1i, N2i });
                resizeMultiDimVector(g_vec_mu_l[key], { N_f, (cheby_degree + 1), N1i, N2i });

                resizeMultiDimVector(g_vec_mu_h_flattened[key], { N_f, flattened_len });
                resizeMultiDimVector(g_vec_mu_l_flattened[key], { N_f, flattened_len });
            }
        }
    }
}

// -----------------------------------------------------------------------------
// Per-(ind, ind_1, flavour) workers (single-core / serial)
// -----------------------------------------------------------------------------
inline void process_lnZ_flavour_for_pair(int ind, int ind_1, int f) {
    const int N1i = N1_vec[ind];
    const int N2i = N2_vec[ind];
    const Real betai = beta_vec[ind_1];

    // Work in k-space Chebyshev coeffs then FFT to r-space
    Complex3DVector g_k_space;
    resizeMultiDimVector(g_k_space, { N1i, N2i, (cheby_degree + 1) });

    for (int n_k1 = 0; n_k1 < N1i; ++n_k1) {
        for (int n_k2 = 0; n_k2 < N2i; ++n_k2) {
            auto partialG0 = [=](Real tau) -> Real {
                return G_0(betai, mu_vec[f] - alpha_shift_vec[f],
                           N1i, N2i, n_k1, n_k2, greens_func_exp_cut, tau);
            };
            Real1DVector g_coeffs =
                chebyshev_decomposition_scalar(partialG0, eps, betai - eps, cheby_degree);

            for (int k = 0; k <= cheby_degree; ++k)
                g_k_space[n_k1][n_k2][k] = Complex(g_coeffs[k], 0.0);
        }
    }

    const std::string key = std::to_string(ind) + "_" + std::to_string(ind_1);
    // Write directly into the flavour slice [f]
    transform_vec(g_k_space, g_vec[key][f]);
    flatten_3D_to_1D(g_vec[key][f], g_vec_flattened[key][f]);
}

inline void process_density_mu_h_flavour_for_pair(int ind, int ind_1, int f) {
    const int N1i = N1_vec[ind];
    const int N2i = N2_vec[ind];
    const Real betai = beta_vec[ind_1];

    Complex3DVector g_k_space;
    resizeMultiDimVector(g_k_space, { N1i, N2i, (cheby_degree + 1) });

    for (int n_k1 = 0; n_k1 < N1i; ++n_k1) {
        for (int n_k2 = 0; n_k2 < N2i; ++n_k2) {
            auto partialG0 = [=](Real tau) -> Real {
                return G_0(betai, mu_h_vec[f] - alpha_shift_vec[f],
                           N1i, N2i, n_k1, n_k2, greens_func_exp_cut, tau);
            };
            Real1DVector g_coeffs =
                chebyshev_decomposition_scalar(partialG0, eps, betai - eps, cheby_degree);

            for (int k = 0; k <= cheby_degree; ++k)
                g_k_space[n_k1][n_k2][k] = Complex(g_coeffs[k], 0.0);
        }
    }

    const std::string key = std::to_string(ind) + "_" + std::to_string(ind_1);
    transform_vec(g_k_space, g_vec_mu_h[key][f]);
    flatten_3D_to_1D(g_vec_mu_h[key][f], g_vec_mu_h_flattened[key][f]);
}

inline void process_density_mu_l_flavour_for_pair(int ind, int ind_1, int f) {
    const int N1i = N1_vec[ind];
    const int N2i = N2_vec[ind];
    const Real betai = beta_vec[ind_1];

    Complex3DVector g_k_space;
    resizeMultiDimVector(g_k_space, { N1i, N2i, (cheby_degree + 1) });

    for (int n_k1 = 0; n_k1 < N1i; ++n_k1) {
        for (int n_k2 = 0; n_k2 < N2i; ++n_k2) {
            auto partialG0 = [=](Real tau) -> Real {
                return G_0(betai, mu_l_vec[f] - alpha_shift_vec[f],
                           N1i, N2i, n_k1, n_k2, greens_func_exp_cut, tau);
            };
            Real1DVector g_coeffs =
                chebyshev_decomposition_scalar(partialG0, eps, betai - eps, cheby_degree);

            for (int k = 0; k <= cheby_degree; ++k)
                g_k_space[n_k1][n_k2][k] = Complex(g_coeffs[k], 0.0);
        }
    }

    const std::string key = std::to_string(ind) + "_" + std::to_string(ind_1);
    transform_vec(g_k_space, g_vec_mu_l[key][f]);
    flatten_3D_to_1D(g_vec_mu_l[key][f], g_vec_mu_l_flattened[key][f]);
}

// -----------------------------------------------------------------------------
// Entry point (single core): loops over (ind, ind_1) then flavours
// -----------------------------------------------------------------------------
inline void load_data() {
    preallocate_multi_flavour_maps();

    const int n_sizes = static_cast<int>(N1_vec.size());
    const int n_betas = static_cast<int>(beta_vec.size());

    for (int ind = 0; ind < n_sizes; ++ind) {
        for (int ind_1 = 0; ind_1 < n_betas; ++ind_1) {

            if (lnZ_mode) {
                for (int f = 0; f < N_f; ++f) {
                    process_lnZ_flavour_for_pair(ind, ind_1, f);
                }
            }

            if (density_mode) {
                for (int f = 0; f < N_f; ++f) {
                    process_density_mu_h_flavour_for_pair(ind, ind_1, f);
                    process_density_mu_l_flavour_for_pair(ind, ind_1, f);
                }
            }
        }
    }
}

#endif // LOAD_DATA_H
