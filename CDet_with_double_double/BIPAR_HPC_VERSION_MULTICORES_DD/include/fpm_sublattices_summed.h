#ifndef FPM_SUBLATTICES_SUMMED_H
#define FPM_SUBLATTICES_SUMMED_H

#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>
#include <chrono>
#include <random>
#include <iomanip>
#include <algorithm>
#include <stdexcept>
#include <type_traits>
#include <mutex>
#include <numeric>

#include "../utility/types.h"
#include "../utility/my_funcs.h"

// Constants
constexpr int FPM_MAXN = 12;
constexpr int FPM_MAXDIM = 2 * FPM_MAXN;

// Precomputed powers
constexpr int pow2_table[FPM_MAXN + 1] = {
    1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096
};
constexpr int pow3_table[FPM_MAXN + 1] = {
    1, 3, 9, 27, 81, 243, 729, 2187, 6561, 19683, 59049, 177147, 531441
};

// Calculate required workspace size at compile time
constexpr size_t calculate_fpm_workspace_size() {
    size_t total_size = 0;
    for (int n = 1; n <= FPM_MAXN; ++n) {
        size_t path_size = 0;
        for (int i = 0; i < n; ++i) {
            size_t dim = 2 * (n - i) - 1;
            if (dim > 0) {
                path_size += dim * dim;
            }
        }
        if (path_size > total_size) {
            total_size = path_size;
        }
    }
    return total_size;
}

constexpr size_t FPM_WORKSPACE_SIZE = calculate_fpm_workspace_size();
constexpr size_t FPM_OUT_RAW_SIZE = pow3_table[FPM_MAXN];

// Type traits for math operations
template<typename T>
struct MathOps;

template<>
struct MathOps<double> {
    static double abs(const double& x) { return std::abs(x); }
    static double epsilon() { return 1e-15; }
};

template<>
struct MathOps<DoubleDouble> {
    static DoubleDouble abs(const DoubleDouble& x) { return DD_ABS(x); }
    static DoubleDouble epsilon() { return DoubleDouble(1e-25); }
};

// Matrix view structure
template<typename T>
struct MatrixView {
    T* data;
    int dim;
    int stride;

    T& at(int r, int c) { return data[r * stride + c]; }
    const T& at(int r, int c) const { return data[r * stride + c]; }

    MatrixView<T> subView(int row_start, int col_start, int new_dim) const {
        return {data + row_start * stride + col_start, new_dim, stride};
    }
};

// Function declarations
template<typename T>
void schur_complement(const MatrixView<T>& M, T pivot, MatrixView<T> dst, T thresh) {
    if (MathOps<T>::abs(pivot) < thresh) {
        for (int i = 0; i < dst.dim; ++i) {
            std::fill(&dst.at(i, 0), &dst.at(i, dst.dim), T(0));
        }
        return;
    }

    const T inv_pivot = T(1) / pivot;
    for (int i = 0; i < dst.dim; ++i) {
        const T s_val = M.at(i + 1, 0) * inv_pivot;
        for (int j = 0; j < dst.dim; ++j) {
            if constexpr (std::is_same_v<T, double>) {
                dst.at(i, j) = std::fma(-s_val, M.at(0, j + 1), M.at(i + 1, j + 1));
            } else {
                dst.at(i, j) = M.at(i + 1, j + 1) - s_val * M.at(0, j + 1);
            }
        }
    }
}

template<typename T>
void dfs_recursive(int depth, int n, const MatrixView<T>& M, T det, int pat_prefix,
                   T thresh, T C_shift, T* workspace_ptr, T* out_raw,
                   std::vector<std::tuple<int, int, int>>& zero_pivots_info) {
    if (depth == n) {
        out_raw[pat_prefix] = det;
        return;
    }

    const int d = 2 * (n - depth);
    const int next_pat_prefix = pat_prefix * 3;

    if (d >= 2) {
        dfs_recursive(depth + 1, n, M.subView(2, 2, d - 2), det, next_pat_prefix + 0,
                     thresh, C_shift, workspace_ptr, out_raw, zero_pivots_info);
    } else {
        dfs_recursive(depth + 1, n, {}, det, next_pat_prefix + 0,
                     thresh, C_shift, workspace_ptr, out_raw, zero_pivots_info);
    }

    T p1 = M.at(0, 0);
    if (MathOps<T>::abs(p1) < thresh) {
        zero_pivots_info.emplace_back(depth, 1, pat_prefix);
        p1 = C_shift;
    }
    if (d > 1) {
        MatrixView<T> schur_mat = {workspace_ptr, d - 1, d - 1};
        schur_complement(M, p1, schur_mat, thresh);
        T* next_workspace_ptr = workspace_ptr + (d - 1) * (d - 1);
        dfs_recursive(depth + 1, n, schur_mat.subView(1, 1, d - 2), det * p1,
                     next_pat_prefix + 1, thresh, C_shift, next_workspace_ptr,
                     out_raw, zero_pivots_info);
    } else {
        dfs_recursive(depth + 1, n, {}, det * p1, next_pat_prefix + 1,
                     thresh, C_shift, workspace_ptr, out_raw, zero_pivots_info);
    }

    if (d >= 2) {
        MatrixView<T> sub_mat = M.subView(1, 1, d - 1);
        T p2 = sub_mat.at(0, 0);
        if (MathOps<T>::abs(p2) < thresh) {
            zero_pivots_info.emplace_back(depth, 2, pat_prefix);
            p2 = C_shift;
        }
        if (d > 2) {
            MatrixView<T> schur_mat = {workspace_ptr, d - 2, d - 2};
            schur_complement(sub_mat, p2, schur_mat, thresh);
            T* next_workspace_ptr = workspace_ptr + (d - 2) * (d - 2);
            dfs_recursive(depth + 1, n, schur_mat, det * p2, next_pat_prefix + 2,
                         thresh, C_shift, next_workspace_ptr, out_raw, zero_pivots_info);
        } else {
            dfs_recursive(depth + 1, n, {}, det * p2, next_pat_prefix + 2,
                         thresh, C_shift, workspace_ptr, out_raw, zero_pivots_info);
        }
    }
}

template<typename T>
void apply_multilinearity_corrections(int n, T* out_raw, T C_shift,
                                    const std::vector<std::tuple<int, int, int>>& zero_pivots_info) {
    auto sorted_info = zero_pivots_info;
    std::sort(sorted_info.begin(), sorted_info.end(), [](const auto& a, const auto& b) {
        return std::get<0>(a) > std::get<0>(b);
    });

    for (const auto& pivot_info : sorted_info) {
        auto [depth, choice, parent_prefix] = pivot_info;
        int n_suffix = n - (depth + 1);
        int suffix_count = pow3_table[n_suffix];
        for (int s = 0; s < suffix_count; ++s) {
            int base_pattern_minor = (parent_prefix * 3 + 0) * suffix_count + s;
            int base_pattern_affected = (parent_prefix * 3 + choice) * suffix_count + s;
            out_raw[base_pattern_affected] -= C_shift * out_raw[base_pattern_minor];
        }
    }
}

template<typename T>
void compute_3n_minors(const std::vector<T>& A, int n, T* workspace, T* out_raw, T tol = T(1e-8)) {
    if (n < 0 || n > FPM_MAXN) throw std::runtime_error("n is out of supported range [0, 12]");
    const int N = 2 * n;
    if (A.size() != (size_t)N * N) throw std::runtime_error("Invalid matrix size");

    // Dynamic threshold and shift calculation
    T sumAbs = T(0);
    for (const T& v : A) {
        sumAbs += MathOps<T>::abs(v);
    }

    T scale = T(1);
    if (N > 0) {
        T calculated_scale = sumAbs / (N * N);
        if (MathOps<T>::abs(calculated_scale) > MathOps<T>::epsilon()) {
            scale = calculated_scale;
        }
    }

    T thresh = tol * scale;
    T C_shift = scale;

    if (thresh < MathOps<T>::epsilon()) {
        thresh = MathOps<T>::epsilon();
    }

    std::vector<T> M_initial = A;
    std::vector<std::tuple<int, int, int>> zero_pivots_info;
    zero_pivots_info.reserve(64); // Reserve some space to avoid reallocations

    MatrixView<T> initial_view = {M_initial.data(), N, N};
    dfs_recursive(0, n, initial_view, T(1), 0, thresh, C_shift, workspace, out_raw, zero_pivots_info);

    if (!zero_pivots_info.empty()) {
        apply_multilinearity_corrections(n, out_raw, C_shift, zero_pivots_info);
    }
}

template<typename T>
void postprocess_mat2pm(const T* up_3n, const T* down_3n, int n, T* output_mat2pm) {
    if (n < 0 || n > FPM_MAXN) throw std::runtime_error("n is out of supported range");
    const unsigned num_minors = pow2_table[n];

    // Zero initialize
    std::fill(output_mat2pm, output_mat2pm + (num_minors - 1), T(0));

    for (int pat = 0; pat < pow3_table[n]; ++pat) {
        unsigned mask = 0;
        int temp_pat = pat;
        for (int i = 0; i < n; ++i) {
            if (temp_pat % 3 != 0) mask |= (1u << i);
            temp_pat /= 3;
        }

        if (mask != 0) {
            output_mat2pm[mask - 1] += up_3n[pat] * down_3n[pat];
        }
    }
}

template<typename T>
void compute_postprocessed_minors(const std::vector<T>& A_up, const std::vector<T>& A_down,
                                 int n, T* workspace, T* out_raw_up, T* out_raw_down, T* result) {
    compute_3n_minors(A_up, n, workspace, out_raw_up);
    compute_3n_minors(A_down, n, workspace + FPM_WORKSPACE_SIZE, out_raw_down);
    postprocess_mat2pm(out_raw_up, out_raw_down, n, result);
}

// Convenience function for legacy code compatibility
template<typename T>
void brute_force_postprocessed(const std::vector<T>& A_up, const std::vector<T>& A_down,
                              int n, std::vector<T>& result) {
    std::vector<T> workspace(2 * FPM_WORKSPACE_SIZE);
    std::vector<T> out_raw(2 * FPM_OUT_RAW_SIZE);

    if (pow2_table[n] > 0) {
        result.resize(pow2_table[n] - 1);
    } else {
        result.clear();
        return;
    }

    compute_postprocessed_minors(A_up, A_down, n, workspace.data(),
                                out_raw.data(), out_raw.data() + FPM_OUT_RAW_SIZE,
                                result.data());
}

#endif //FPM_SUBLATTICES_SUMMED_H