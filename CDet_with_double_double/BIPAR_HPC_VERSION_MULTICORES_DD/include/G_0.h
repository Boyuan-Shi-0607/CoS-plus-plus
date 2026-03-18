#ifndef G_0_H
#define G_0_H

#include "../utility/types.h"
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include "../utility/my_funcs.h"

// MP_Real version of multiply2x2
inline MP_Complex2DVector multiply2x2_mp(const MP_Complex2DVector &A, const MP_Complex2DVector &B) {
    MP_Complex2DVector C(2, MP_Complex1DVector(2, MP_Complex(0, 0)));
    for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
            MP_Complex sum(0, 0);
            for(int k=0; k<2; k++){
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = sum;
        }
    }
    return C;
}

// MP_Real version of dagger2x2
inline MP_Complex2DVector dagger2x2_mp(const MP_Complex2DVector &A) {
    MP_Complex2DVector A_dag(2, MP_Complex1DVector(2, MP_Complex(0, 0)));
    for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
            A_dag[i][j] = conj(A[j][i]);
        }
    }
    return A_dag;
}

// MP_Real version of h_matrix
inline MP_Complex2DVector h_matrix_mp(const int n_k1, const int n_k2,
                                     const MP_Real& Delta_s, const MP_Real& mu_s,
                                     const MP_Real& H_s,
                                     const int N1, const int N2) {
    MP_Complex2DVector h(2, MP_Complex1DVector(2, MP_Complex(0, 0)));

    // Energy dispersion with MP_Real
    const MP_Real r_nk1 = MP_Real(n_k1);
    const MP_Real r_nk2 = MP_Real(n_k2);
    const MP_Real r_N1 = MP_Real(N1);
    const MP_Real r_N2 = MP_Real(N2);

    const MP_Real two_pi = MP_Real(2) * boost::multiprecision::acos(MP_Real(-1));
    const MP_Real alpha = two_pi * r_nk1 / r_N1;
    const MP_Real beta = two_pi * r_nk2 / r_N2;

    MP_Complex E(-1, 0);
    E -= MP_Complex(MP_COS(alpha), MP_SIN(alpha));
    E -= MP_Complex(MP_COS(beta), MP_SIN(beta));

    const MP_Real a = Delta_s - mu_s + H_s;
    const MP_Real d = -Delta_s - mu_s - H_s;

    h[0][0] = MP_Complex(a, 0);
    h[0][1] = E;
    h[1][0] = conj(E);
    h[1][1] = MP_Complex(d, 0);

    return h;
}

// MP_Real version of eigh2x2
struct EigenResult2x2_mp {
    std::vector<MP_Real> eigenvalues;
    MP_Complex2DVector eigenvectors;
};

inline EigenResult2x2_mp eigh2x2_mp(const MP_Complex2DVector &h) {
    const MP_Real a = real(h[0][0]);
    const MP_Real d = real(h[1][1]);
    const MP_Complex b = h[0][1];
    const MP_Real bmag = MP_SQRT(real(b)*real(b) + imag(b)*imag(b));

    const MP_Real sum = a + d;
    const MP_Real diff = a - d;
    const MP_Real R = MP_SQRT(diff*diff + MP_Real(4) * bmag*bmag);

    const MP_Real w1 = MP_Real(0.5) * (sum - R);
    const MP_Real w2 = MP_Real(0.5) * (sum + R);

    const MP_Real phi = boost::multiprecision::atan2(imag(b), real(b));
    const MP_Real cos2theta = (d - a) / R;
    const MP_Real sin2theta = (MP_Real(2) * bmag) / R;

    const MP_Real angle2 = MP_ATAN2(sin2theta, cos2theta);
    const MP_Real theta = MP_Real(0.5) * angle2;

    const MP_Real cth = MP_COS(theta);
    const MP_Real sth = MP_SIN(theta);

    MP_Complex2DVector K(2, MP_Complex1DVector(2, MP_Complex(0, 0)));

    K[0][0] = MP_Complex(cth, 0);
    K[0][1] = MP_Complex(MP_COS(phi)*sth, MP_SIN(phi)*sth);
    K[1][0] = -MP_Complex(MP_COS(-phi)*sth, MP_SIN(-phi)*sth);
    K[1][1] = MP_Complex(cth, 0);

    EigenResult2x2_mp result;
    result.eigenvalues.resize(2);
    if(w1 < w2) {
        result.eigenvalues[0] = w1;
        result.eigenvalues[1] = w2;
    } else {
        result.eigenvalues[0] = w2;
        result.eigenvalues[1] = w1;
    }

    result.eigenvectors = K;
    return result;
}

// MP_Real version of G_0
inline MP_Complex2DVector G_0_mp(const MP_Real& beta_s, const MP_Real& H_s, const MP_Real& Delta_s, const MP_Real& mu_s,
                                const int N1, const int N2, const int n_k1, const int n_k2, const MP_Real& cut,
                                const MP_Real& tau) {

    const auto h_matrix_s = h_matrix_mp(n_k1, n_k2, Delta_s, mu_s, H_s, N1, N2);
    const auto eigh_s = eigh2x2_mp(h_matrix_s);
    const auto xi_s = eigh_s.eigenvalues;
    const auto K_s = eigh_s.eigenvectors;
    const auto K_dagger_s = dagger2x2_mp(K_s);

    const MP_Real xi0 = xi_s[0];
    const MP_Real xi1 = xi_s[1];

    const MP_Real n_F0 = MP_Real(1) / (MP_Real(1) + MP_EXP(beta_s * xi0));
    const MP_Real n_F1 = MP_Real(1) / (MP_Real(1) + MP_EXP(beta_s * xi1));

    MP_Complex2DVector M(2, MP_Complex1DVector(2, MP_Complex(0, 0)));

    if (tau > 0) {
        MP_Complex b00, b11;
        if (beta_s * xi0 < cut) {
            b00 = MP_Complex(MP_EXP((beta_s - tau) * xi0) * n_F0);
        } else {
            b00 = MP_Complex(MP_EXP(-tau * xi0));
        }
        if (beta_s * xi1 < cut) {
            b11 = MP_Complex(MP_EXP((beta_s - tau) * xi1) * n_F1);
        } else {
            b11 = MP_Complex(MP_EXP(-tau * xi1));
        }

        M[0][0] = b00;
        M[0][1] = MP_Complex(0, 0);
        M[1][0] = MP_Complex(0, 0);
        M[1][1] = b11;

        const MP_Complex2DVector temp = multiply2x2_mp(K_s, M);
        MP_Complex2DVector result = multiply2x2_mp(temp, K_dagger_s);
        return result;
    } else {
        MP_Complex c00, c11;
        if (beta_s * xi0 < cut) {
            c00 = MP_Complex(MP_EXP(-tau * xi0) * n_F0);
        } else {
            c00 = MP_Complex(MP_EXP(-(tau + beta_s) * xi0));
        }
        if (beta_s * xi1 < cut) {
            c11 = MP_Complex(MP_EXP(-tau * xi1) * n_F1);
        } else {
            c11 = MP_Complex(MP_EXP(-(tau + beta_s) * xi1));
        }
        M[0][0] = c00;
        M[0][1] = MP_Complex(0, 0);
        M[1][0] = MP_Complex(0, 0);
        M[1][1] = c11;

        const MP_Complex2DVector temp = multiply2x2_mp(K_s, M);
        MP_Complex2DVector result = multiply2x2_mp(temp, K_dagger_s);

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++){
                result[i][j] = MP_Complex(-1, 0) * result[i][j];
            }
        }
        return result;
    }
}


#endif