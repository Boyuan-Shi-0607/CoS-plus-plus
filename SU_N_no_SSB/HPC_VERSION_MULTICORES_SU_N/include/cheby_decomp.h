#ifndef CHEBY_DECOMP_H
#define CHEBY_DECOMP_H

#include "../utility/types.h"
#include "../utility/my_funcs.h"
#include <vector>
#include <functional>
#include <iostream>
#include <complex>
#include <cmath>
#include <cassert>

// -----------------------------------------------------------------------------
// Real versions for scalar functions
// -----------------------------------------------------------------------------

inline void compute_T_polys(const Real& t, const int deg, Real1DVector &out_T) {
    out_T.resize(deg+1);
    if(deg >= 0) out_T[0] = 1.0;
    if(deg >= 1) out_T[1] = t;

    for(int m = 1; m < deg; m++) {
        out_T[m+1] = 2.0 * t * out_T[m] - out_T[m-1];
    }
}

// Discrete Chebyshev fit for scalar Real values
static Real1DVector discrete_chebfit_scalar(const Real1DVector& t_nodes,
                                            const Real1DVector& y_data,
                                            const int deg) {
    const int n_nodes = deg + 1;
    Real1DVector c(deg+1, 0.0);

    for(int m=0; m <= deg; m++) {
        Real factor = (m==0) ? 1.0 / Real(n_nodes)
                             : 2.0 / Real(n_nodes);

        Real sum_val = 0.0;
        for(int k=0; k < n_nodes; k++){
            if (m == 0) {
                sum_val += y_data[k];
            } else if (m == 1) {
                sum_val += y_data[k] * t_nodes[k];
            }
            else {
                const Real t = t_nodes[k];
                Real T0 = 1.0;
                Real T1 = t;
                for(int r=1; r < m; r++){
                    const Real T2 = 2.0 * t * T1 - T0;
                    T0 = T1;
                    T1 = T2;
                }
                sum_val += y_data[k] * T1;
            }
        }
        c[m] = factor * sum_val;
    }
    return c;
}

// Chebyshev decomposition for scalar Real functions
inline Real1DVector chebyshev_decomposition_scalar(
    const std::function<Real(Real)> &f,
    const Real& lower,
    const Real& upper,
    const int deg) {

    const int n_nodes = deg + 1;

    Real1DVector y_values(n_nodes, 0.0);
    Real1DVector t_nodes(n_nodes, 0.0);
    Real1DVector x_nodes(n_nodes, 0.0);

    const Real my_pi = M_PI;
    for (int k=0; k < n_nodes; k++) {
        const Real numerator = 2.0 * Real(k) + 1.0;
        const Real denominator = 2.0 * Real(n_nodes);
        const Real arg = my_pi * (numerator / denominator);
        t_nodes[k] = std::cos(arg);
    }

    for (int k=0; k < n_nodes; k++) {
        x_nodes[k] = lower + (upper - lower) * ((t_nodes[k] + 1.0) * 0.5);
        y_values[k] = f(x_nodes[k]);
    }

    Real1DVector t_nodes_mapped(n_nodes);
    for (int k=0; k < n_nodes; k++){
        t_nodes_mapped[k] = 2.0 * (x_nodes[k] - lower) / (upper - lower) - 1.0;
    }

    return discrete_chebfit_scalar(t_nodes_mapped, y_values, deg);
}

#endif //CHEBY_DECOMP_H