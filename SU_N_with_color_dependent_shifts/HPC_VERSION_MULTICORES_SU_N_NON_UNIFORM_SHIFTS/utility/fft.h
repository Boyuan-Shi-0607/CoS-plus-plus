#ifndef FFT_H
#define FFT_H

#include <cmath>
#include <complex>
#include <cstddef>
#include "types.h"
#include "my_funcs.h"

// -----------------------------------------------------------------------------
// Real FFT implementation (using standard complex FFT)
// -----------------------------------------------------------------------------

inline Complex1DVector fft(const Complex1DVector& x) {
    std::size_t N = x.size();
    if (N <= 1) return x;

    Complex1DVector X(N);
    Real two_pi = 2.0 * M_PI;

    for (std::size_t k = 0; k < N; ++k) {
        X[k] = Complex(0, 0);
        for (std::size_t n = 0; n < N; ++n) {
            Real angle = -two_pi * Real(k) * Real(n) / Real(N);
            Complex w(std::cos(angle), std::sin(angle));
            X[k] += x[n] * w;
        }
    }
    return X;
}

inline Complex1DVector ifft(const Complex1DVector& X) {
    std::size_t N = X.size();
    if (N <= 1) return X;

    Complex1DVector x(N);
    Real two_pi = 2.0 * M_PI;
    Real scale = 1.0 / Real(N);

    for (std::size_t n = 0; n < N; ++n) {
        x[n] = Complex(0, 0);
        for (std::size_t k = 0; k < N; ++k) {
            Real angle = two_pi * Real(k) * Real(n) / Real(N);
            Complex w(std::cos(angle), std::sin(angle));
            x[n] += X[k] * w;
        }
        x[n] *= scale;
    }
    return x;
}

// Real 2D FFT
inline void fft2d(Complex2DVector& data) {
    std::size_t M = data.size();
    if (M == 0) return;
    std::size_t N = data[0].size();
    if (N == 0) return;

    for (std::size_t m = 0; m < M; ++m) {
        data[m] = fft(data[m]);
    }

    for (std::size_t n = 0; n < N; ++n) {
        Complex1DVector col(M);
        for (std::size_t m = 0; m < M; ++m) {
            col[m] = data[m][n];
        }

        Complex1DVector colFFT = fft(col);

        for (std::size_t m = 0; m < M; ++m) {
            data[m][n] = colFFT[m];
        }
    }
}

inline void ifft2d(Complex2DVector& data) {
    std::size_t M = data.size();
    if (M == 0) return;
    std::size_t N = data[0].size();
    if (N == 0) return;

    for (std::size_t m = 0; m < M; ++m) {
        data[m] = ifft(data[m]);
    }

    for (std::size_t n = 0; n < N; ++n) {
        Complex1DVector col(M);
        for (std::size_t m = 0; m < M; ++m) {
            col[m] = data[m][n];
        }

        Complex1DVector colIFFT = ifft(col);

        for (std::size_t m = 0; m < M; ++m) {
            data[m][n] = colIFFT[m];
        }
    }
}

// Transform function for Real - modified for scalar values (3D instead of 5D)
inline void transform_vec(const Complex3DVector &vec2, Real3DVector &vec1) {
    size_t N1 = vec2.size();
    if (N1 == 0) return;
    size_t N2 = vec2[0].size();
    if (N2 == 0) return;
    size_t M  = vec2[0][0].size();
    if (M == 0) return;

    // For each Chebyshev coefficient
    for (size_t m = 0; m < M; m++) {
        Complex2DVector temp2D(N1, Complex1DVector(N2));

        // Gather the m-th coefficient from all k-points
        for (size_t n1 = 0; n1 < N1; n1++) {
            for (size_t n2 = 0; n2 < N2; n2++) {
                temp2D[n1][n2] = vec2[n1][n2][m];
            }
        }

        // Perform 2D FFT
        fft2d(temp2D);

        // Store the real part back
        for (size_t n1 = 0; n1 < N1; n1++) {
            for (size_t n2 = 0; n2 < N2; n2++) {
                vec1[m][n1][n2] = temp2D[n1][n2].real() / (Real(N1) * Real(N2));
            }
        }
    }
}

#endif // FFT_H