#ifndef VEGAS_MAP_NAMESPACE_H
#define VEGAS_MAP_NAMESPACE_H

#include "../utility/types.h"
#include <vector>
#include <string>
#include <cmath>
#include <numeric>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <limits>

namespace VEGAS {

    // --- Core container: holds only what is needed to apply a histogram-built map ---
    struct MapData {
        int N_DIM = 0;
        int N_INTERVALS = 0;         // N
        int N_EDGES = 0;             // N + 1

        std::vector<Real> min;       // per-dim domain min (x-space)
        std::vector<Real> max;       // per-dim domain max (x-space)

        // Per-dimension x-grid derived from histogram
        std::vector<std::vector<Real>> x_edges;  // size: N_DIM x (N+1)
        std::vector<std::vector<Real>> dx_steps; // size: N_DIM x N
    };

    // --- Init a blank map (uniform edges). You will call BuildFromHistogram(...) next. ---
    inline MapData InitializeMap(int NDIM, int Intervals, const Real* min_vals, const Real* max_vals) {
        MapData map;
        map.N_DIM = NDIM;
        map.N_INTERVALS = Intervals;
        map.N_EDGES = Intervals + 1;

        map.min.assign(min_vals, min_vals + NDIM);
        map.max.assign(max_vals, max_vals + NDIM);

        map.x_edges.assign(NDIM, std::vector<Real>(map.N_EDGES, 0.0));
        map.dx_steps.assign(NDIM, std::vector<Real>(map.N_INTERVALS, 0.0));

        // Start with uniform grid (used if a dimension gets a zero histogram)
        for (int d = 0; d < NDIM; ++d) {
            Real a = map.min[d], b = map.max[d];
            Real step = (b - a) / map.N_INTERVALS;
            for (int k = 0; k <= map.N_INTERVALS; ++k) map.x_edges[d][k] = a + step * k;
            for (int k = 0; k <  map.N_INTERVALS; ++k) map.dx_steps[d][k] = map.x_edges[d][k+1] - map.x_edges[d][k];
        }
        return map;
    }

    // --- Build a single dimension's x-edges from a histogram h (length N_INTERVALS) ---
    // h_i is the (pseudo-)histogram value on the uniform source bin [b_{i-1}, b_i] in x-space.
    // Robust to zeros: if Z == 0 or a local h_j == 0, we fall back locally (or globally) to uniform.
    inline void BuildFromHistogram1D(MapData& map, int dim, const std::vector<Real>& h) {
        if (dim < 0 || dim >= map.N_DIM) throw std::out_of_range("BuildFromHistogram1D: bad dim");
        if ((int)h.size() != map.N_INTERVALS) throw std::runtime_error("histogram length != N_INTERVALS");

        const int N = map.N_INTERVALS;
        const Real a = map.min[dim], b = map.max[dim];
        const Real DeltaB = (b - a) / N;

        // 1) cumulative mass C (with C_0 = 0)
        std::vector<Real> C(N+1, 0.0);
        for (int i = 1; i <= N; ++i) C[i] = C[i-1] + h[i-1] * DeltaB;
        const Real Z = C[N];

        // If the whole histogram is zero, keep uniform and return
        if (Z == Real(0)) {
            // uniform already set by InitializeMap
            return;
        }

        // Source uniform edges b_i
        auto b_edge = [&](int i)->Real { return a + (b - a) * (static_cast<Real>(i) / N); };

        // 2) targets T_k = (k/N) * Z, then invert C to get x_k
        map.x_edges[dim][0] = a;
        map.x_edges[dim][N] = b;

        int j = 1; // index of source bin carrying the current target
        for (int k = 1; k < N; ++k) {
            const Real T = (static_cast<Real>(k) / N) * Z;

            // advance j so that C_{j-1} <= T <= C_j
            while (j <= N && C[j] < T) ++j;
            if (j > N) j = N;

            const Real Cjm1 = C[j-1];
            const Real hj   = h[j-1];

            if (hj > Real(0)) {
                // within that uniform source bin, invert linearly
                const Real frac = (T - Cjm1) / hj; // length inside the bin
                map.x_edges[dim][k] = b_edge(j-1) + frac; // since DeltaB * hj * (delta/DeltaB)/hj collapses to delta
            } else {
                // hj == 0: target lies in an empty bin; pin to its left edge
                map.x_edges[dim][k] = b_edge(j-1);
            }
        }

        // 3) finalize dx_k
        for (int k = 0; k < N; ++k) {
            map.dx_steps[dim][k] = map.x_edges[dim][k+1] - map.x_edges[dim][k];
            // basic degeneracy guard: ensure strictly positive with tiny epsilon fallback
            if (!(map.dx_steps[dim][k] > Real(0))) {
                // replace with a tiny positive width so J != 0; preserves monotonicity
                map.dx_steps[dim][k] = std::max((b-a) * Real(1e-12), std::numeric_limits<Real>::epsilon());
                map.x_edges[dim][k+1] = map.x_edges[dim][k] + map.dx_steps[dim][k];
            }
        }
        // ensure last edge is exactly b (accumulated rounding)
        map.x_edges[dim][N] = b;
        map.dx_steps[dim][N-1] = map.x_edges[dim][N] - map.x_edges[dim][N-1];
    }

    // --- Convenience: build all dimensions at once (hists[d] is size N_INTERVALS) ---
    inline void BuildFromHistograms(MapData& map, const std::vector<std::vector<Real>>& hists) {
        if ((int)hists.size() != map.N_DIM) throw std::runtime_error("BuildFromHistograms: wrong outer size");
        for (int d = 0; d < map.N_DIM; ++d) BuildFromHistogram1D(map, d, hists[d]);
    }

    // --- Interval helpers in y-space (y in [0,1]^N_DIM) ---
    inline std::vector<int> GetIntervalID(const MapData& map, const std::vector<Real>& y) {
        std::vector<int> res; res.reserve(map.N_DIM);
        for (int i = 0; i < map.N_DIM; ++i) {
            Real yi = y[i];
            if (yi < Real(0)) yi = Real(0);
            if (yi >= Real(1)) yi = std::nextafter(Real(1), Real(0)); // clamp to [0,1)
            int id = static_cast<int>(std::floor(yi * map.N_INTERVALS));
            if (id < 0) id = 0;
            if (id >= map.N_INTERVALS) id = map.N_INTERVALS - 1;
            res.push_back(id);
        }
        return res;
    }

    inline std::vector<Real> GetIntervalOffset(const MapData& map, const std::vector<Real>& y) {
        std::vector<int> ID = GetIntervalID(map, y);
        std::vector<Real> off; off.reserve(map.N_DIM);
        for (int i = 0; i < map.N_DIM; ++i) {
            Real yi = std::min(std::max(y[i], Real(0)), std::nextafter(Real(1), Real(0)));
            off.push_back(yi * map.N_INTERVALS - ID[i]);
        }
        return off;
    }

    // --- Map y -> x (piecewise-linear, per your slide's Eq. (within-bin)) ---
    inline std::vector<Real> GetX(const MapData& map, const std::vector<Real>& y) {
        std::vector<Real> x; x.reserve(map.N_DIM);
        std::vector<int> ID = GetIntervalID(map, y);
        std::vector<Real> t = GetIntervalOffset(map, y);
        for (int i = 0; i < map.N_DIM; ++i) {
            const int k = ID[i];
            x.push_back(map.x_edges[i][k] + map.dx_steps[i][k] * t[i]);
        }
        return x;
    }

    // --- Jacobian and log10 Jacobian: J = Π_i [N * Δx_k^(i)] ---
    inline Real GetJac(const MapData& map, const std::vector<Real>& y) {
        std::vector<int> ID = GetIntervalID(map, y);
        Real jac = Real(1);
        for (int i = 0; i < map.N_DIM; ++i) {
            const int k = ID[i];
            jac *= map.N_INTERVALS * map.dx_steps[i][k];
        }
        return jac;
    }

    inline Real GetLogJac(const MapData& map, const std::vector<Real>& y) {
        std::vector<int> ID = GetIntervalID(map, y);
        Real log10jac = Real(0);
        for (int i = 0; i < map.N_DIM; ++i) {
            const int k = ID[i];
            // log10(N * Δx_k)
            log10jac += std::log10(map.N_INTERVALS * map.dx_steps[i][k]);
        }
        return log10jac;
    }

    // --- Introspection / I/O ---
    inline void PrintEdges(const MapData& map) {
        std::cout << "Mapped x-edges per dimension:\n";
        for (int d = 0; d < map.N_DIM; ++d) {
            std::cout << "  dim " << d << " edges:";
            for (int k = 0; k <= map.N_INTERVALS; ++k) std::cout << "\t" << map.x_edges[d][k];
            std::cout << "\n  dim " << d << " dx:   ";
            for (int k = 0; k <  map.N_INTERVALS; ++k) std::cout << "\t" << map.dx_steps[d][k];
            std::cout << "\n";
        }
    }

    // Prints y_k = k/N, x_k, and J_k = N * Δx_k for each k
    inline void PrintMap(const MapData& map) {
        for (int d = 0; d < map.N_DIM; ++d) {
            std::cout << "dim " << d << " : y_k, x_k, J_k\n";
            for (int k = 0; k < map.N_INTERVALS; ++k) {
                Real yk   = static_cast<Real>(k) / map.N_INTERVALS;
                Real xk   = map.x_edges[d][k];
                Real Jk   = map.N_INTERVALS * map.dx_steps[d][k];
                std::cout << "  " << yk << "\t" << xk << "\t" << Jk << "\n";
            }
            // last edge
            std::cout << "  " << 1.0 << "\t" << map.x_edges[d][map.N_INTERVALS] << "\t-\n";
        }
    }

    inline void SaveEdges(const MapData& map, const std::string& basename) {
        for (int d = 0; d < map.N_DIM; ++d) {
            std::ofstream out(basename + std::to_string(d) + ".txt");
            for (int k = 0; k <= map.N_INTERVALS; ++k) out << map.x_edges[d][k] << "\n";
        }
    }

} // namespace VEGAS

#endif // VEGAS_MAP_NAMESPACE_H
