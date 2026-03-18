#ifndef LNZ_CALCULATION_H
#define LNZ_CALCULATION_H

#include "../utility/types.h"
#include <cmath>
#include <vector>
#include <numeric>

#include <limits>
#include <stdexcept>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <iomanip>

using MP_Real    = boost::multiprecision::cpp_bin_float_100;
using MP_Complex = boost::multiprecision::cpp_complex_100;

// Get PI with full precision
inline MP_Real get_pi() {
    static MP_Real pi = boost::multiprecision::acos(MP_Real(-1));
    return pi;
}

// Calculate binomial coefficient C(n, k)
MP_Real binomial_coefficient(int n, int k) {
    if (k > n || k < 0) return 0;
    if (k == 0 || k == n) return 1;
    
    MP_Real result = 1;
    for (int i = 0; i < k; ++i) {
        result *= (n - i);
        result /= (i + 1);
    }
    return result;
}

// Dispersion function
inline MP_Real dispersion(int n_k1, int n_k2, int N1, int N2, MP_Real mu_s, MP_Real alpha) {
    
    MP_Real cos1 = boost::multiprecision::cos(MP_Real(2.0) * get_pi() * n_k1 / N1);
    MP_Real cos2 = boost::multiprecision::cos(MP_Real(2.0) * get_pi() * n_k2 / N2);
    
    return - MP_Real(2.0) * (cos1 + cos2) - mu_s + alpha;
}

// Calculate lnZ for a single (n_k1, n_k2) pair
inline MP_Real lnZ_analytical(int n_k1, int n_k2, int N1, int N2, int N_f,
                              MP_Real beta_s, MP_Real mu_s, MP_Real alpha) {

    // Since we're doing brute force sum of exp due to high precision
    MP_Real sum = 1.0; // First term: 1.0 * exp(0.0) = 1.0
    
    // Add remaining terms
    for (int i = 1; i <= N_f; ++i) {
        MP_Real coeff = binomial_coefficient(N_f, i);
        MP_Real exponent = -i * beta_s * dispersion(n_k1, n_k2, N1, N2, mu_s, alpha);
        sum += coeff * boost::multiprecision::exp(exponent);
    }
    
    return boost::multiprecision::log(sum);
}

// Main function to calculate total lnZ
inline MP_Real calculate_lnZ(int N1, int N2, int N_f, MP_Real beta_s, MP_Real mu_s, MP_Real alpha) {
    MP_Real lnZ = 0.0;
    
    // Sum over all (n_k1, n_k2) combinations
    for (int n_k1 = 0; n_k1 < N1; ++n_k1) {
        for (int n_k2 = 0; n_k2 < N2; ++n_k2) {
            lnZ += lnZ_analytical(n_k1, n_k2, N1, N2, N_f, beta_s, mu_s, alpha);
        }
    }
    
    return lnZ / (MP_Real(N1) * MP_Real(N2));
}

#endif // LNZ_CALCULATION_H