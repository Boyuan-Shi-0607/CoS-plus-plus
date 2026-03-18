#ifndef FFT_H
#define FFT_H

#include <cmath>
#include <complex>
#include <cstddef>
#include "types.h"
#include "my_funcs.h"

// -----------------------------------------------------------------------------
// MP_Real FFT implementation (DFT for arbitrary precision)
// -----------------------------------------------------------------------------

inline MP_Complex1DVector fft_mp(const MP_Complex1DVector& x) {
    std::size_t N = x.size();
    if (N <= 1) return x;

    MP_Complex1DVector X(N);
    MP_Real two_pi = MP_Real(2) * boost::multiprecision::acos(MP_Real(-1)); // 2*pi

    for (std::size_t k = 0; k < N; ++k) {
        X[k] = MP_Complex(0, 0);
        for (std::size_t n = 0; n < N; ++n) {
            MP_Real angle = -two_pi * MP_Real(k) * MP_Real(n) / MP_Real(N);
            MP_Complex w(MP_COS(angle), MP_SIN(angle));
            X[k] += x[n] * w;
        }
    }
    return X;
}

inline MP_Complex1DVector ifft_mp(const MP_Complex1DVector& X) {
    std::size_t N = X.size();
    if (N <= 1) return X;

    MP_Complex1DVector x(N);
    MP_Real two_pi = 2 * boost::multiprecision::acos(MP_Real(-1)); // 2*pi
    MP_Real scale = MP_Real(1) / MP_Real(N);

    for (std::size_t n = 0; n < N; ++n) {
        x[n] = MP_Complex(0, 0);
        for (std::size_t k = 0; k < N; ++k) {
            MP_Real angle = two_pi * MP_Real(k) * MP_Real(n) / MP_Real(N);
            MP_Complex w(MP_COS(angle), MP_SIN(angle));
            x[n] += X[k] * w;
        }
        x[n] *= MP_Complex(scale, 0);  // Fix: multiply by MP_Complex
    }
    return x;
}

// MP_Real 2D FFT
inline void fft2d_mp(MP_Complex2DVector& data) {
    std::size_t M = data.size();
    if (M == 0) return;
    std::size_t N = data[0].size();
    if (N == 0) return;

    for (std::size_t m = 0; m < M; ++m) {
        data[m] = fft_mp(data[m]);
    }

    for (std::size_t n = 0; n < N; ++n) {
        MP_Complex1DVector col(M);
        for (std::size_t m = 0; m < M; ++m) {
            col[m] = data[m][n];
        }

        MP_Complex1DVector colFFT = fft_mp(col);

        for (std::size_t m = 0; m < M; ++m) {
            data[m][n] = colFFT[m];
        }
    }
}

inline void ifft2d_mp(MP_Complex2DVector& data) {
    std::size_t M = data.size();
    if (M == 0) return;
    std::size_t N = data[0].size();
    if (N == 0) return;

    for (std::size_t m = 0; m < M; ++m) {
        data[m] = ifft_mp(data[m]);
    }

    for (std::size_t n = 0; n < N; ++n) {
        MP_Complex1DVector col(M);
        for (std::size_t m = 0; m < M; ++m) {
            col[m] = data[m][n];
        }

        MP_Complex1DVector colIFFT = ifft_mp(col);

        for (std::size_t m = 0; m < M; ++m) {
            data[m][n] = colIFFT[m];
        }
    }
}

// Transform function for MP_Real
inline void transform_vec_mp(const MP_Complex5DVector &vec2, MP_Real5DVector &vec1) {
    size_t N1 = vec2.size();
    if (N1 == 0) return;
    size_t N2 = vec2[0].size();
    if (N2 == 0) return;
    size_t M  = vec2[0][0].size();
    if (M == 0) return;

    for (size_t m = 0; m < M; m++) {
        for (size_t i = 0; i < 2; i++) {
            for (size_t j = 0; j < 2; j++) {
                MP_Complex2DVector temp2D(N1, MP_Complex1DVector(N2));
                for (size_t n1 = 0; n1 < N1; n1++) {
                    for (size_t n2 = 0; n2 < N2; n2++) {
                        temp2D[n1][n2] = vec2[n1][n2][m][i][j];
                    }
                }

                fft2d_mp(temp2D);

                for (size_t n1 = 0; n1 < N1; n1++) {
                    for (size_t n2 = 0; n2 < N2; n2++) {
                        vec1[m][i][j][n1][n2] = real(temp2D[n1][n2]) / (MP_Real(N1) * MP_Real(N2));
                    }
                }
            }
        }
    }
}

#endif // FFT_H