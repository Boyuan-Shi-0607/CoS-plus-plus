/*
===============================================================================
Free multi-flavour fermions on a 2D lattice — grand-canonical ln Z
===============================================================================

Setup
-----
- Flavours: f = 1, …, N_f
- Momenta: k ≡ (n_{k1}, n_{k2}) with n_{k1}=0..N_1-1, n_{k2}=0..N_2-1
- Number operators: n_{k f} ∈ {0,1}
- One-particle energies (flavour-dependent shift α_f):
    ε_{k,f} = -2[ cos(2π n_{k1}/N_1) + cos(2π n_{k2}/N_2) ] - μ + α_f

Hamiltonian
-----------
    H = Σ_k Σ_{f=1}^{N_f} ε_{k,f} n_{k f}

Because H is diagonal in (k,f), the grand-canonical partition function factorizes:
For fixed k and each flavour f,
    Z_{k,f} = Σ_{n=0}^1 e^{-β ε_{k,f} n} = 1 + e^{-β ε_{k,f}}.

Thus per momentum k,
    Z_k = Π_{f=1}^{N_f} (1 + e^{-β ε_{k,f}}).

Total partition function:
    Z = Π_k Z_k
      = Π_k Π_{f=1}^{N_f} (1 + e^{-β ε_{k,f}}).

Taking logs:
    ln Z = Σ_k Σ_{f=1}^{N_f} ln(1 + e^{-β ε_{k,f}}).

In code we therefore compute:
    ln Z = sum over n_{k1}, n_{k2}, and flavours f of
           log1pexp(x_{k,f}) with x_{k,f} := -β ε_{k,f},

using a numerically stable "log1p(exp(x))" evaluation.

We also commonly report ln Z per momentum (or per site) as:
    (ln Z) / (N_1 N_2).

===============================================================================
*/

#ifndef LNZ_FREE_PART
#define LNZ_FREE_PART

#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/math/constants/constants.hpp>

using MP_Real    = boost::multiprecision::cpp_bin_float_100;
using MP_Complex = boost::multiprecision::cpp_complex_100;

// π with full precision (cached)
inline const MP_Real& get_pi() {
    static const MP_Real pi = boost::math::constants::pi<MP_Real>();
    return pi;
}

// Numerically stable log(1 + exp(x)) for high precision types.
inline MP_Real log1p_exp(const MP_Real& x) {
    using boost::multiprecision::exp;
    using boost::multiprecision::log1p;
    using boost::multiprecision::log;

    // Handle common extremes explicitly
    if (x <= MP_Real(-60)) {
        // exp(x) ~ 0 in 100-bit precision
        return MP_Real(0);
    } else if (x >= MP_Real(60)) {
        // log(1 + e^x) ~ x for large x
        return x;
    }

    if (x > 0) {
        // log(1 + e^x) = x + log(1 + e^{-x})
        return x + log1p(exp(-x));
    } else {
        // log(1 + e^x) = log(1 + e^x)
        return log1p(exp(x));
    }
}

// 2D lattice cosine dispersion part (flavour-independent)
inline MP_Real lattice_part_2d(int n_k1, int n_k2, int N1, int N2) {
    using boost::multiprecision::cos;
    MP_Real theta1 = MP_Real(2) * get_pi() * MP_Real(n_k1) / MP_Real(N1);
    MP_Real theta2 = MP_Real(2) * get_pi() * MP_Real(n_k2) / MP_Real(N2);
    return MP_Real(-2) * (cos(theta1) + cos(theta2));
}

// Full single-particle energy ε_{k,f}
inline MP_Real epsilon_kf(int n_k1, int n_k2, int N1, int N2,
                          const MP_Real& mu, const MP_Real& alpha_f) {
    // ε_{k,f} = (-2[cos(..)+cos(..)]) - μ + α_f
    return lattice_part_2d(n_k1, n_k2, N1, N2) - mu + alpha_f;
}

// ln Z contribution from a single momentum k, summing over flavours
inline MP_Real lnZ_k(int n_k1, int n_k2, int N1, int N2,
                     const MP_Real& beta, const MP_Real& mu,
                     const std::vector<MP_Real>& alpha_f) {
    MP_Real sum_logs = 0;
    for (const auto& a_f : alpha_f) {
        MP_Real eps = epsilon_kf(n_k1, n_k2, N1, N2, mu, a_f);
        MP_Real x   = -beta * eps;                 // x = -β ε_{k,f}
        sum_logs   += log1p_exp(x);                // ln(1 + e^{x})
    }
    return sum_logs; // = ln Z_k
}

// Total ln Z (grand-canonical), summed over all k, f
inline MP_Real calculate_lnZ_total(int N1, int N2,
                                   const MP_Real& beta, const MP_Real& mu,
                                   const std::vector<MP_Real>& alpha_f) {
    if (N1 <= 0 || N2 <= 0) {
        throw std::invalid_argument("N1 and N2 must be positive.");
    }
    if (alpha_f.empty()) {
        throw std::invalid_argument("alpha_f must contain at least one flavour.");
    }

    MP_Real lnZ = 0;
    for (int n_k1 = 0; n_k1 < N1; ++n_k1) {
        for (int n_k2 = 0; n_k2 < N2; ++n_k2) {
            lnZ += lnZ_k(n_k1, n_k2, N1, N2, beta, mu, alpha_f);
        }
    }
    return lnZ;
}

// Convenience: ln Z per momentum/site
inline MP_Real calculate_lnZ_per_site(int N1, int N2,
                                      const MP_Real& beta, const MP_Real& mu,
                                      const std::vector<MP_Real>& alpha_f) {
    MP_Real lnZ = calculate_lnZ_total(N1, N2, beta, mu, alpha_f);
    return lnZ / (MP_Real(N1) * MP_Real(N2));
}

#endif

