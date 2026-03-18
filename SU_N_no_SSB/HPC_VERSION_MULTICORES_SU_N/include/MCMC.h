#ifndef MCMC_VEGAS_H
#define MCMC_VEGAS_H

#include "../utility/types.h"
#include "online_stats.h"
#include "../utility/my_funcs.h"

#include <vector>
#include <random>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <iomanip>
#include <limits>
#include <tuple>
#include <chrono>
#include <algorithm>

// ===== VEGAS mapping API (refactored header you provided) =====
#include "vegas_map_namespace.h" // MapData, InitializeMap, BuildFromHistogram1D, GetX, GetLogJac

// ===== Single compile-time switch: auto-train & use VEGAS (no main.cpp changes) =====
#define USE_VEGAS_TRAINING 1

namespace MC {

    // =========================
    // Tunables / constants
    // =========================
    constexpr Real P_GAUSSIAN  = 0.4;
    constexpr Real P_UNIFORM   = 0.2;
    constexpr Real P_HISTOGRAM = 0.4;

    constexpr int  TRACKED_DIMENSION         = 8;
    constexpr unsigned long long TRACK_EVERY = 1000;

    constexpr int  DEFAULT_BINS = 32;     // 32 bins everywhere (proposal + VEGAS)
    const     Real EPSILON_THRESHOLD = 5e-15;

#if USE_VEGAS_TRAINING
    // Pre-training controls (can be adjusted by the app prior to calling run_*)
    inline int VEGAS_ORDER_FOR_SIGN = 4;                // which order to learn from
    inline std::vector<int> VEGAS_SIGN_HIST_DIMS = {1, 5, 9}; // which dims to record
    inline unsigned long long VEGAS_TRAIN_STEPS  = (1ULL << 21); // 2^21 MCMC steps
    inline Real VEGAS_TRAIN_SIGMA = 0.3;               // proposal scale in pre-train
#endif

    // =========================
    // Histogram (shared by warm-up & VEGAS pre-training)
    // =========================
    struct Histogram1D {
        int B;
        Real bin_width;
        std::vector<unsigned long long> counts;
        std::vector<Real> probs;      // normalized after finalize()
        bool finalized{false};

        explicit Histogram1D(int bins = DEFAULT_BINS)
            : B(bins),
              bin_width(Real(1) / Real(bins)),
              counts(bins, 0ULL),
              probs(bins, 0.0) {}

        inline int bin_index(Real y) const {
            if (y <= 0.0) return 0;
            if (y >= 1.0) return B - 1;
            int i = static_cast<int>(y / bin_width);
            if (i < 0) i = 0;
            if (i >= B) i = B - 1;
            return i;
        }

        inline void observe(Real y) {
            counts[bin_index(y)]++;
        }

        inline void finalize(Real add_frac = 0.01) {
            unsigned long long T = 0;
            for (auto c : counts) T += c;
            const Real add = add_frac * static_cast<Real>(T);
            probs.assign(B, 0.0);
            Real S = 0.0;
            for (int b = 0; b < B; ++b) {
                probs[b] = static_cast<Real>(counts[b]) + add;
                S += probs[b];
            }
            if (S <= 0.0) {
                for (int b = 0; b < B; ++b) probs[b] = Real(1) / Real(B);
            } else {
                for (int b = 0; b < B; ++b) probs[b] /= S;
            }
            finalized = true;
        }

        inline Real density(Real y) const {
            if (!finalized) return 0.0;
            return probs[bin_index(y)] * Real(B); // piecewise-constant density
        }

        inline Real sample(std::mt19937 &g) const {
            std::uniform_real_distribution<Real> U(0.0, 1.0);
            if (!finalized) return U(g);
            Real u = U(g), cum = 0.0;
            int b = 0;
            for (; b < B; ++b) {
                cum += probs[b];
                if (u <= cum) break;
            }
            if (b >= B) b = B - 1;
            Real left = b * bin_width;
            Real y = left + U(g) * bin_width;
            if (y < 0) y = 0;
            if (y > 1) y = 1;
            return y;
        }
    };

    // Save histogram to CSV-like text (used by main::run_all)
    inline void save_histogram(const Histogram1D &hist, const std::string &filename) {
        std::ofstream f(filename);
        if (!f) { std::cerr << "Cannot open " << filename << " for write\n"; return; }
        f << std::fixed << std::setprecision(15);
        f << "bin_index,count,probability\n";
        for (int i = 0; i < hist.B; ++i) {
            f << i << "," << hist.counts[i] << "," << hist.probs[i] << "\n";
        }
    }

    // =========================
    // Helpers
    // =========================
    inline Real Uniform01(std::mt19937 &g) {
        std::uniform_real_distribution<Real> U(0.0, 1.0);
        return U(g);
    }

    inline int RandomInt(int a, int b, std::mt19937 &g) {
        std::uniform_int_distribution<int> D(a, b);
        return D(g);
    }

    inline Real sign_func(Real x) {
        return (x > 0) - (x < 0);
    }

    inline bool is_even_support(Real s) {
        if (s <= 0.0) return true;
        if (s >= 1.0) return false;
        int i = static_cast<int>(std::floor(s * 100.0)); // 100 bins of 0.01
        if (i < 0) i = 0;
        if (i > 99) i = 99;
        return (i % 2) == 0;
    }

    inline Real block_variance(const std::vector<Real> &means) {
        if (means.size() <= 1) return 0.0;
        Real mu = 0.0;
        for (auto v : means) mu += v;
        mu /= Real(means.size());
        Real v = 0.0;
        for (auto x : means) {
            Real d = x - mu;
            v += d * d;
        }
        v /= Real(means.size());
        v /= Real(means.size() - 1);
        return v;
    }

    // =========================
    // Function types + evaluators (nullable VEGAS map inside)
    // =========================
    using RealFunc = std::function<Real(const Real *, int)>;

    // Single-function evaluator (sign path)
    inline Real eval_f_n(const RealFunc &f,
                         const std::vector<Real> &Y,
                         int n,
                         VEGAS::MapData *vm /*nullable*/) {
        if (!vm) return f(Y.data(), n);
        auto x = VEGAS::GetX(*vm, Y);
        Real v = f(x.data(), n);
        if (v == 0) return 0;
        Real lj = VEGAS::GetLogJac(*vm, Y);  // base-10
        Real L  = std::log10(std::fabs(v)) + lj;
        if (L > 50.) throw std::runtime_error("func*Jac > 1e50");
        return (v > 0 ? 1 : -1) * std::pow(Real(10), L);
    }

    // Two-function evaluator (unified path). Maps only the first DIM coords for a branch.
    inline Real eval_unified(const RealFunc &measure,
                             const RealFunc &ref,
                             int DIM,
                             const std::vector<Real> &Y,
                             int n,
                             VEGAS::MapData *vm_measure /*nullable*/,
                             VEGAS::MapData *vm_ref /*nullable*/) {
        const bool even = is_even_support(Y.back());

        if (even) {
            if (!vm_measure) return measure(Y.data(), n);
            std::vector<Real> yfirst(Y.begin(), Y.begin() + DIM);
            auto  x  = VEGAS::GetX(*vm_measure, yfirst);
            Real  v  = measure(x.data(), n);
            if (v == 0) return 0;
            Real  lj = VEGAS::GetLogJac(*vm_measure, yfirst);
            Real  L  = std::log10(std::fabs(v)) + lj;
            if (L > 50.) throw std::runtime_error("func*Jac > 1e50");
            return (v > 0 ? 1 : -1) * std::pow(Real(10), L);
        } else {
            if (!vm_ref) return ref(Y.data(), n);
            std::vector<Real> yfirst(Y.begin(), Y.begin() + DIM);
            auto  x  = VEGAS::GetX(*vm_ref, yfirst);
            Real  v  = ref(x.data(), n);
            if (v == 0) return 0;
            Real  lj = VEGAS::GetLogJac(*vm_ref, yfirst);
            Real  L  = std::log10(std::fabs(v)) + lj;
            if (L > 50.) throw std::runtime_error("func*Jac > 1e50");
            return (v > 0 ? 1 : -1) * std::pow(Real(10), L);
        }
    }

    // =========================
    // MCMC steps (VEGAS-aware via nullable vm / vm_measure / vm_ref)
    // =========================
    inline auto MC_step_sign(
        const std::vector<Real> &Y_old,
        Real w_old,
        Real o_old,
        Real sigma,
        int order_num,
        const RealFunc &func,
        std::mt19937 &gen,
        int &selected_direction,
        const std::vector<Histogram1D> *warm_hist = nullptr,
        bool use_hist_after_warmup = false,
        VEGAS::MapData *vm /*nullable*/ = nullptr) {

        auto Y_new = Y_old;
        const int DIM = (int)Y_old.size();
        selected_direction = RandomInt(0, DIM - 1, gen);

        enum class K { GAUSS, UNIF, HIST };
        K k = K::GAUSS;
        Real u = Uniform01(gen);
        if (!use_hist_after_warmup || !warm_hist) {
            k = (u < 0.7) ? K::GAUSS : K::UNIF;
        } else {
            if (u < P_GAUSSIAN) k = K::GAUSS;
            else if (u < P_GAUSSIAN + P_UNIFORM) k = K::UNIF;
            else k = K::HIST;
        }

        if (k == K::GAUSS) {
            std::normal_distribution<Real> ND(0.0, sigma);
            Real p = std::fmod(Y_old[selected_direction] + ND(gen), Real(1));
            if (p < 0) p += 1;
            Y_new[selected_direction] = p;
        } else if (k == K::UNIF) {
            std::uniform_real_distribution<Real> U(0.0, 1.0);
            Y_new[selected_direction] = U(gen);
        } else {
            Y_new[selected_direction] = (*warm_hist)[selected_direction].sample(gen);
        }

        const Real fn  = eval_f_n(func, Y_new, order_num, vm);
        const Real w_n = std::fabs(fn);

        Real ratio = (w_old < EPSILON_THRESHOLD) ?
                     ((w_n < EPSILON_THRESHOLD) ? 0.0 : 1.0) : (w_n / w_old);

        if (k == K::HIST) {
            Real hx = (*warm_hist)[selected_direction].density(Y_old[selected_direction]);
            Real hy = (*warm_hist)[selected_direction].density(Y_new[selected_direction]);
            if (hy < EPSILON_THRESHOLD) return std::make_tuple(Y_old, w_old, o_old, false);
            ratio *= (hx < EPSILON_THRESHOLD ? 0.0 : (hx / hy));
        }

        if (ratio > Uniform01(gen)) {
            return std::make_tuple(Y_new, w_n, sign_func(fn), true);
        }
        return std::make_tuple(Y_old, w_old, o_old, false);
    }

    inline auto MC_step_nd_unified(
        const std::vector<Real> &Y_old,   // size DIM+1
        Real w_old,
        Real o_odd_old,                   // ref indicator
        Real o_even_old,                  // measure indicator
        Real sigma,
        int order_num,
        const RealFunc &measure,
        const RealFunc &ref,
        int DIM,                          // visible dims
        std::mt19937 &gen,
        int &selected_direction,
        const std::vector<Histogram1D> *warm_hist = nullptr,
        bool use_hist_after_warmup = false,
        VEGAS::MapData *vm_measure /*nullable*/ = nullptr,
        VEGAS::MapData *vm_ref     /*nullable*/ = nullptr) {

        auto Y_new = Y_old;
        const int Dp1 = (int)Y_old.size();
        selected_direction = RandomInt(0, Dp1 - 1, gen);

        enum class K { GAUSS, UNIF, HIST };
        K k = K::GAUSS;
        Real u = Uniform01(gen);
        if (!use_hist_after_warmup || !warm_hist) {
            k = (u < 0.7) ? K::GAUSS : K::UNIF;
        } else {
            if (u < P_GAUSSIAN) k = K::GAUSS;
            else if (u < P_GAUSSIAN + P_UNIFORM) k = K::UNIF;
            else k = K::HIST;
        }

        if (k == K::GAUSS) {
            std::normal_distribution<Real> ND(0.0, sigma);
            Real p = std::fmod(Y_old[selected_direction] + ND(gen), Real(1));
            if (p < 0) p += 1;
            Y_new[selected_direction] = p;
        } else if (k == K::UNIF) {
            std::uniform_real_distribution<Real> U(0.0, 1.0);
            Y_new[selected_direction] = U(gen);
        } else {
            Y_new[selected_direction] = (*warm_hist)[selected_direction].sample(gen);
        }

        const Real v  = eval_unified(measure, ref, DIM, Y_new, order_num, vm_measure, vm_ref);
        const Real wn = std::fabs(v);

        Real ratio = (w_old < EPSILON_THRESHOLD) ?
                     ((wn < EPSILON_THRESHOLD) ? 0.0 : 1.0) : (wn / w_old);

        if (k == K::HIST) {
            Real hx = (*warm_hist)[selected_direction].density(Y_old[selected_direction]);
            Real hy = (*warm_hist)[selected_direction].density(Y_new[selected_direction]);
            if (hy < EPSILON_THRESHOLD) return std::make_tuple(Y_old, w_old, o_odd_old, o_even_old, false);
            ratio *= (hx < EPSILON_THRESHOLD ? 0.0 : (hx / hy));
        }

        if (ratio > Uniform01(gen)) {
            const bool even = is_even_support(Y_new.back());
            Real o_even_new = even ? 1.0 : 0.0;
            Real o_odd_new  = even ? 0.0 : 1.0;
            return std::make_tuple(Y_new, wn, o_odd_new, o_even_new, true);
        }
        return std::make_tuple(Y_old, w_old, o_odd_old, o_even_old, false);
    }

    // =========================
    // Good starting points (use current evaluators; maps may be nullptr)
    // =========================
    inline bool find_good_starting_point_sign(
        std::vector<Real> &Y0, Real &w0, Real &o0,
        const RealFunc &f, Real sigma, int order_n, int DIM, std::mt19937 &gen,
        int max_attempts = 200, int test_steps = 100, Real min_acc = 0.03,
        VEGAS::MapData *vm /*nullable*/ = nullptr) {

        for (int a = 0; a < max_attempts; ++a) {
            std::uniform_real_distribution<Real> U(0.0, 1.0);
            Y0.assign(DIM, 0.5);
            for (int i = 0; i < DIM; ++i) Y0[i] = U(gen);
            Real v  = eval_f_n(f, Y0, order_n, vm);
            w0 = std::fabs(v);
            o0 = sign_func(v);

            std::vector<Real> Yt = Y0;
            Real wt = w0, ot = o0;
            int acc = 0;
            for (int s = 0; s < test_steps; ++s) {
                int dir = -1;
                auto r = MC_step_sign(Yt, wt, ot, sigma, order_n, f, gen, dir, nullptr, false, vm);
                Yt = std::get<0>(r);
                wt = std::get<1>(r);
                ot = std::get<2>(r);
                if (std::get<3>(r)) ++acc;
            }
            if (Real(acc) / Real(test_steps) > min_acc) return true;
        }
        return false;
    }

    inline bool find_good_starting_point_unified(
        std::vector<Real> &Y0, Real &w0, Real &o_odd0, Real &o_even0,
        const RealFunc &measure, const RealFunc &ref,
        int DIM, Real sigma, int order_n, int Dp1, std::mt19937 &gen,
        int max_attempts = 200, int test_steps = 100, Real min_acc = 0.03,
        VEGAS::MapData *vm_measure = nullptr, VEGAS::MapData *vm_ref = nullptr) {

        for (int a = 0; a < max_attempts; ++a) {
            std::uniform_real_distribution<Real> U(0.0, 1.0);
            Y0.assign(Dp1, 0.5);
            for (int i = 0; i < Dp1; ++i) Y0[i] = U(gen);

            const bool even = is_even_support(Y0.back());
            Real v = eval_unified(measure, ref, DIM, Y0, order_n, vm_measure, vm_ref);
            w0 = std::fabs(v);
            o_even0 = even ? 1.0 : 0.0;
            o_odd0  = even ? 0.0 : 1.0;

            std::vector<Real> Yt = Y0;
            Real wt = w0, o1 = o_odd0, o2 = o_even0;
            int acc = 0;
            for (int s = 0; s < test_steps; ++s) {
                int dir = -1;
                auto r = MC_step_nd_unified(Yt, wt, o1, o2, sigma, order_n,
                                            measure, ref, DIM, gen, dir,
                                            nullptr, false, vm_measure, vm_ref);
                Yt = std::get<0>(r);
                wt = std::get<1>(r);
                o1 = std::get<2>(r);
                o2 = std::get<3>(r);
                if (std::get<4>(r)) ++acc;
            }
            if (Real(acc) / Real(test_steps) > min_acc) return true;
        }
        return false;
    }

    // =========================
    // VEGAS pre-training (MCMC-based, like warm-up)
    // =========================
#if USE_VEGAS_TRAINING
    // MCMC train histograms for selected dims using the *sign* step (no histogram kernel).
    inline std::vector<Histogram1D>
    vegas_train_hist_sign_mcmc(const RealFunc &f,
                               int DIM,
                               int order_n,
                               std::mt19937 &gen,
                               unsigned long long steps,
                               Real sigma,
                               const std::vector<int> &dims,
                               int bins = DEFAULT_BINS) {

        // Good start
        std::vector<Real> Y(DIM);
        Real w = 0.0, o = 0.0;
        if (!find_good_starting_point_sign(Y, w, o, f, sigma, order_n, DIM, gen)) {
            std::uniform_real_distribution<Real> U(0.0, 1.0);
            for (int i = 0; i < DIM; ++i) Y[i] = U(gen);
            Real v = f(Y.data(), order_n);
            w = std::fabs(v);
            o = sign_func(v);
        }

        // Hist containers
        std::vector<Histogram1D> H;
        H.reserve(dims.size());
        for (std::size_t j = 0; j < dims.size(); ++j) H.emplace_back(bins);

        // Run a pure GAUSS/UNIF MCMC (no histogram kernel) and record current state
        for (unsigned long long t = 0; t < steps; ++t) {
            int dir = -1;
            auto r = MC_step_sign(Y, w, o, sigma, order_n, f, gen,
                                  dir, nullptr, /*use_hist_after_warmup*/ false,
                                  /*vm*/ nullptr);
            Y = std::get<0>(r);
            w = std::get<1>(r);
            o = std::get<2>(r);

            for (std::size_t j = 0; j < dims.size(); ++j) {
                int d = dims[j];
                if (0 <= d && d < DIM) H[j].observe(Y[d]);
            }
        }

        for (auto &h : H) h.finalize(0.01);
        return H;
    }

    // Assign map per dim: multiples of ord → uniform; otherwise cycle across trained groups.
    inline VEGAS::MapData make_map_from_hist(int DIM, int ord, const std::vector<Histogram1D> &Hs) {
        std::vector<Real> mins(DIM, 0.0), maxs(DIM, 1.0);
        auto vm = VEGAS::InitializeMap(DIM, DEFAULT_BINS, mins.data(), maxs.data());
        const int K = static_cast<int>(Hs.size());
        if (K == 0) return vm;

        std::vector<std::vector<Real>> groups;
        groups.reserve(K);
        for (const auto &h : Hs) groups.push_back(h.probs);

        auto choose = [&](int d) -> int {
            if (ord <= 0) return -1;
            if (d % ord == 0) return -1;            // uniform at 0, ord, 2*ord, ...
            int block = (d - 1) / ord;              // 0,1,2,...
            return block % K;                       // cycle trained groups
        };

        for (int d = 0; d < DIM; ++d) {
            int g = choose(d);
            if (g >= 0) VEGAS::BuildFromHistogram1D(vm, d, groups[g]);
        }
        return vm;
    }
#endif // USE_VEGAS_TRAINING

    // =========================
    // SIGN driver (now returns dir ratios, tracked, histograms too)
    // =========================
    inline std::tuple<
        Real,                       // sign_mean
        std::vector<Real>,          // block_variances
        Real,                       // final sigma
        Real,                       // overall accept_ratio
        std::vector<Real>,          // per-direction accept ratios
        std::vector<Real>,          // tracked values (Y[TRACKED_DIMENSION])
        std::vector<Histogram1D>    // finalized warm-up histograms
    >
    run_sign(const RealFunc &f,
             Real initial_sigma,
             int order_num,
             int DIM,
             std::mt19937 &gen,
             unsigned int iter_bits,
             unsigned int warm_bits) {

        const unsigned long long N_warm  = (1ULL << warm_bits);
        const unsigned long long N_main  = (1ULL << iter_bits);
        const unsigned long long N_steps = N_warm + N_main;

        // Auto-train VEGAS map before warm-up, if enabled
        VEGAS::MapData vm;
        VEGAS::MapData *vm_ptr = nullptr;
#if USE_VEGAS_TRAINING
        {
            auto Hs = vegas_train_hist_sign_mcmc(
                f, DIM, VEGAS_ORDER_FOR_SIGN, gen,
                VEGAS_TRAIN_STEPS, VEGAS_TRAIN_SIGMA,
                VEGAS_SIGN_HIST_DIMS, DEFAULT_BINS
            );
            vm = make_map_from_hist(DIM, order_num, Hs);
            vm_ptr = &vm;
        }
#endif

        // Good start (uses vm_ptr if non-null)
        std::vector<Real> Y(DIM);
        Real w = 0.0, o = 0.0;
        (void)find_good_starting_point_sign(Y, w, o, f, initial_sigma, order_num, DIM, gen,
                                            200, 100, 0.03, vm_ptr);

        OnlineStatsT<Real> stats(iter_bits > 10 ? (iter_bits - 10 + 1) : 1);
        std::vector<Histogram1D> warm(DIM, Histogram1D(DEFAULT_BINS));

        std::vector<Real> dir_acc(DIM, 0.0), dir_att(DIM, 0.0);
        std::vector<Real> tracked;

        Real sigma = initial_sigma;
        Real warm_acc = 0.0, main_acc = 0.0;
        const unsigned long long adapt_int = std::max(N_warm / 4, 100ULL);

        auto t0 = std::chrono::high_resolution_clock::now();

        for (unsigned long long step = 0; step < N_steps; ++step) {
            int dir = -1;
            bool after = (step >= N_warm);
            auto r = MC_step_sign(Y, w, o, sigma, order_num, f, gen, dir,
                                  &warm, after, vm_ptr);
            Y = std::get<0>(r);
            w = std::get<1>(r);
            o = std::get<2>(r);
            bool acc = std::get<3>(r);

            if (!after) {
                for (int d = 0; d < DIM; ++d) warm[d].observe(Y[d]);
                if (acc) warm_acc += 1.0;
                if ((step + 1) % adapt_int == 0 && step > 0) {
                    Real ar = warm_acc / Real(adapt_int);
                    if (ar < 0.34) sigma /= 1.5;
                    else if (ar > 0.7) sigma *= 1.9;
                    else if (ar > 0.5) sigma *= 1.4;
                    warm_acc = 0.0;
                }
                if (step + 1 == N_warm) {
                    for (auto &h : warm) h.finalize();
                }
            } else {
                stats.update(o);
                if (dir >= 0 && dir < DIM) {
                    dir_att[dir] += 1.0;
                    if (acc) { dir_acc[dir] += 1.0; main_acc += 1.0; }
                }
                if (TRACKED_DIMENSION < DIM && (step - N_warm) % TRACK_EVERY == 0) {
                    tracked.push_back(Y[TRACKED_DIMENSION]);
                }
            }
        }

        auto t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<Real> dt = t1 - t0;

        std::vector<Real> block_vars;
        for (int i = 0; i < stats.getMaxBlocks(); ++i) {
            const auto &bm = stats.getBlockMeans(i);
            block_vars.push_back(bm.empty() ? 0.0 : block_variance(bm));
        }

        std::vector<Real> dir_ratio(DIM, 0.0);
        for (int d = 0; d < DIM; ++d) {
            dir_ratio[d] = (dir_att[d] > 0.0) ? (dir_acc[d] / dir_att[d]) : 0.0;
        }

        Real accept_ratio = main_acc / Real(N_main);
        Real sign_mean = stats.getMean();

        std::cout << "\nMC finished in " << dt.count() << " s"
                  << " | accept=" << accept_ratio
                  << " | sign(f)=" << sign_mean << "\n";

        return std::make_tuple(sign_mean, block_vars, sigma, accept_ratio,
                               dir_ratio, tracked, warm);
    }

    // =========================
    // UNIFIED driver (measure/ref)
    // =========================
    inline std::tuple<
        Real,                               // measure / ref
        std::vector<Real>,                  // ref block variances
        std::vector<Real>,                  // measure block variances
        Real,                               // final sigma
        Real,                               // overall accept ratio
        std::vector<Real>,                  // per-direction accept ratios
        std::vector<Real>,                  // tracked values
        std::vector<Histogram1D>            // warm-up histograms
    >
    run_with_maps(const RealFunc &measure,
                  const RealFunc &ref,
                  Real initial_sigma,
                  int order_num,
                  int DIM,
                  std::mt19937 &gen,
                  unsigned int iter_bits,
                  unsigned int warm_bits) {

        const int Dp1 = DIM + 1;
        const unsigned long long N_warm  = (1ULL << warm_bits);
        const unsigned long long N_main  = (1ULL << iter_bits);
        const unsigned long long N_steps = N_warm + N_main;

        // Auto-train separate measure/ref maps (before warm-up), if enabled
        VEGAS::MapData vm_measure, vm_ref;
        VEGAS::MapData *pM = nullptr, *pR = nullptr;
#if USE_VEGAS_TRAINING
        {
            // Train with sign-style MCMC on each branch independently
            auto Hm = vegas_train_hist_sign_mcmc(
                measure, DIM, VEGAS_ORDER_FOR_SIGN, gen,
                VEGAS_TRAIN_STEPS, VEGAS_TRAIN_SIGMA,
                VEGAS_SIGN_HIST_DIMS, DEFAULT_BINS
            );
            auto Hr = vegas_train_hist_sign_mcmc(
                ref, DIM, VEGAS_ORDER_FOR_SIGN, gen,
                VEGAS_TRAIN_STEPS, VEGAS_TRAIN_SIGMA,
                VEGAS_SIGN_HIST_DIMS, DEFAULT_BINS
            );

            vm_measure = make_map_from_hist(DIM, order_num, Hm);
            vm_ref     = make_map_from_hist(DIM, order_num, Hr);
            pM = &vm_measure;
            pR = &vm_ref;
        }
#endif

        // Good start (uses pM/pR if present)
        std::vector<Real> Y(Dp1);
        Real w = 0.0, o1 = 0.0, o2 = 0.0;
        (void)find_good_starting_point_unified(
            Y, w, o1, o2, measure, ref, DIM, initial_sigma, order_num, Dp1, gen,
            200, 100, 0.03, pM, pR
        );

        OnlineStatsT<Real> Sref(iter_bits > 10 ? (iter_bits - 10 + 1) : 1);
        OnlineStatsT<Real> Smeas(iter_bits > 10 ? (iter_bits - 10 + 1) : 1);
        std::vector<Histogram1D> warm(Dp1, Histogram1D(DEFAULT_BINS));
        std::vector<Real> tracked, dir_acc(Dp1, 0.0), dir_att(Dp1, 0.0);

        Real sigma = initial_sigma, warm_acc = 0.0, main_acc = 0.0;
        const unsigned long long adapt_int = std::max(N_warm / 4, 100ULL);

        auto t0 = std::chrono::high_resolution_clock::now();

        for (unsigned long long step = 0; step < N_steps; ++step) {
            int dir = -1;
            bool after = (step >= N_warm);
            auto r = MC_step_nd_unified(Y, w, o1, o2, sigma, order_num,
                                        measure, ref, DIM, gen, dir,
                                        &warm, after, pM, pR);
            Y  = std::get<0>(r);
            w  = std::get<1>(r);
            o1 = std::get<2>(r);
            o2 = std::get<3>(r);
            bool acc = std::get<4>(r);

            if (!after) {
                for (int d = 0; d < Dp1; ++d) warm[d].observe(Y[d]);
                if (acc) warm_acc += 1.0;
                if ((step + 1) % adapt_int == 0 && step > 0) {
                    Real ar = warm_acc / Real(adapt_int);
                    if (ar < 0.34) sigma /= 1.5;
                    else if (ar > 0.7) sigma *= 1.9;
                    else if (ar > 0.5) sigma *= 1.4;
                    warm_acc = 0.0;
                }
                if (step + 1 == N_warm) {
                    for (auto &h : warm) h.finalize();
                }
            } else {
                Sref.update(o1);
                Smeas.update(o2);
                if (dir >= 0 && dir < Dp1) {
                    dir_att[dir] += 1.0;
                    if (acc) { dir_acc[dir] += 1.0; main_acc += 1.0; }
                }
                if (TRACKED_DIMENSION < Dp1 && (step - N_warm) % TRACK_EVERY == 0)
                    tracked.push_back(Y[TRACKED_DIMENSION]);
            }
        }

        auto t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<Real> dt = t1 - t0;

        std::vector<Real> ref_vars, meas_vars;
        for (int i = 0; i < Sref.getMaxBlocks(); ++i) {
            const auto &bm = Sref.getBlockMeans(i);
            ref_vars.push_back(bm.empty() ? 0.0 : block_variance(bm));
        }
        for (int i = 0; i < Smeas.getMaxBlocks(); ++i) {
            const auto &bm = Smeas.getBlockMeans(i);
            meas_vars.push_back(bm.empty() ? 0.0 : block_variance(bm));
        }

        std::vector<Real> dir_ratio(Dp1, 0.0);
        for (int d = 0; d < Dp1; ++d) {
            dir_ratio[d] = (dir_att[d] > 0.0) ? (dir_acc[d] / dir_att[d]) : 0.0;
        }

        Real accept_ratio = main_acc / Real(N_main);
        Real mean_odd  = Sref.getMean();
        Real mean_even = Smeas.getMean();
        if (std::fabs(mean_odd) < EPSILON_THRESHOLD)
            throw std::runtime_error("mean_odd ≈ 0");
        Real result = mean_even / mean_odd;

        std::cout << "\nUnified MC finished in " << dt.count() << " s"
                  << " | accept=" << accept_ratio
                  << " | measure/ref=" << result << "\n";

        return std::make_tuple(result, ref_vars, meas_vars, sigma, accept_ratio,
                               dir_ratio, tracked, warm);
    }

} // namespace MC

#endif // MCMC_VEGAS_H
