#ifndef G_0_H
#define G_0_H

#include "../utility/types.h"
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include "../utility/my_funcs.h"

// Real version of h_matrix - now returns a scalar energy
inline Real h_matrix(const int n_k1, const int n_k2,
                     const Real& mu_s,
                     const int N1, const int N2) {
    // Energy dispersion with Real
    const Real r_nk1 = Real(n_k1);
    const Real r_nk2 = Real(n_k2);
    const Real r_N1 = Real(N1);
    const Real r_N2 = Real(N2);

    const Real two_pi = 2.0 * M_PI;
    const Real k1 = two_pi * r_nk1 / r_N1;
    const Real k2 = two_pi * r_nk2 / r_N2;

    // Simple energy dispersion for square lattice
    Real energy = -2.0 * (std::cos(k1) + std::cos(k2)) - mu_s;

    return energy;
}

// Real version of G_0 - now returns a scalar
inline Real G_0(const Real& beta_s, const Real& mu_s,
                const int N1, const int N2, const int n_k1, const int n_k2,
                const Real& cut, const Real& tau) {

    const Real h_val = h_matrix(n_k1, n_k2, mu_s, N1, N2);

    if (tau > 0) {
        if (beta_s * h_val < cut) {
            return std::exp((beta_s - tau) * h_val) / (1.0 + std::exp(beta_s * h_val));
        } else {
            return std::exp(-tau * h_val);
        }
    } else {
        if (beta_s * h_val < cut) {
            return -std::exp(-tau * h_val) / (1.0 + std::exp(beta_s * h_val));
        } else {
            return -std::exp(-(tau + beta_s) * h_val);
        }
    }
}

#endif