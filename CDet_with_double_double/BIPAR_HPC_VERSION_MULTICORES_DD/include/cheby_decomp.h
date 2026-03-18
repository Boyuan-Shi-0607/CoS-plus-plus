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
// MP_Real versions
// -----------------------------------------------------------------------------

inline void compute_T_polys_mp(const MP_Real& t, const int deg, MP_Real1DVector &out_T) {
    out_T.resize(deg+1);
    if(deg >= 0) out_T[0] = MP_Real(1);
    if(deg >= 1) out_T[1] = t;

    for(int m = 1; m < deg; m++) {
        out_T[m+1] = MP_Real(2) * t * out_T[m] - out_T[m-1];
    }
}

static MP_Complex1DVector discrete_chebfit_mp(const MP_Real1DVector& t_nodes,
                                              const MP_Complex1DVector& y_data,
                                              const int deg) {
    const int n_nodes = deg + 1;
    MP_Complex1DVector c(deg+1, MP_Complex(0, 0));

    for(int m=0; m <= deg; m++) {
        MP_Real factor = (m==0) ? MP_Real(1) / MP_Real(n_nodes)
                                : MP_Real(2) / MP_Real(n_nodes);

        MP_Complex sum_val(0, 0);
        for(int k=0; k < n_nodes; k++){
            if (m == 0) {
                sum_val += y_data[k] * MP_Complex(MP_Real(1), MP_Real(0));
            } else if (m == 1) {
                sum_val += y_data[k] * MP_Complex(t_nodes[k], MP_Real(0));
            }
            else {
                const MP_Real t = t_nodes[k];
                MP_Real T0 = MP_Real(1);
                MP_Real T1 = t;
                for(int r=1; r < m; r++){
                    const MP_Real T2 = MP_Real(2) * t * T1 - T0;
                    T0 = T1;
                    T1 = T2;
                }
                sum_val += y_data[k] * MP_Complex(T1, MP_Real(0));
            }
        }
        c[m] = MP_Complex(factor) * MP_Complex(sum_val);
    }
    return c;
}

inline MP_Complex3DVector chebyshev_decomposition_matrix_mp(
    const std::function<MP_Complex2DVector(MP_Real)> &f,
    const MP_Real& lower,
    const MP_Real& upper,
    const int deg) {

    const int n_nodes = deg + 1;
    const MP_Complex2DVector dummy = f(lower);
    const int n = static_cast<int>(dummy.size());
    if (n == 0) {
        std::cerr << "f(x) size is 0.\n" << std::endl;
    }
    const int n_cols = static_cast<int>(dummy[0].size());

    MP_Complex3DVector y_values(n_nodes, MP_Complex2DVector(n, MP_Complex1DVector(n_cols, MP_Complex(0, 0))));
    MP_Real1DVector t_nodes(n_nodes, MP_Real(0));
    MP_Real1DVector x_nodes(n_nodes, MP_Real(0));

    const MP_Real my_pi = boost::multiprecision::acos(MP_Real(-1));
    for (int k=0; k < n_nodes; k++) {
        const MP_Real numerator = MP_Real(2) * MP_Real(k) + MP_Real(1);
        const MP_Real denominator = MP_Real(2) * MP_Real(n_nodes);
        const MP_Real arg = my_pi * (numerator / denominator);
        t_nodes[k] = MP_COS(arg);
    }

    for (int k=0; k < n_nodes; k++) {
        x_nodes[k] = lower + (upper - lower) * ((t_nodes[k] + MP_Real(1)) * MP_Real(0.5));
        y_values[k] = f(x_nodes[k]);
    }

    MP_Real1DVector t_nodes_mapped(n_nodes);
    for (int k=0; k < n_nodes; k++){
        t_nodes_mapped[k] = MP_Real(2) * (x_nodes[k] - lower) / (upper - lower) - MP_Real(1);
    }

    MP_Complex3DVector coeffs(
        deg+1,
        MP_Complex2DVector(n, MP_Complex1DVector(n_cols, MP_Complex(0, 0)))
    );

    MP_Complex1DVector tmp_y(n_nodes);

    for(int i=0; i < n; i++){
        for(int j=0; j < n_cols; j++){
            for(int k=0; k < n_nodes; k++){
                tmp_y[k] = y_values[k][i][j];
            }
            MP_Complex1DVector c = discrete_chebfit_mp(t_nodes_mapped, tmp_y, deg);
            for(int m=0; m <= deg; m++){
                coeffs[m][i][j] = c[m];
            }
        }
    }

    return coeffs;
}

#endif //CHEBY_DECOMP_H