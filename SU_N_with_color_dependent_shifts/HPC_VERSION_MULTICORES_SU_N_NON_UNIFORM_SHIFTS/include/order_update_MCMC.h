#ifndef MCMC_H
#define MCMC_H

#include <vector>
#include <array>
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <functional>
#include <chrono>
#include <map>

#include "parameters.h"

// Minimal standalone typedefs
using Real = double;

// ===== VEGAS mapping API (assumed to be included) =====
#include "vegas_map_namespace.h"

// ===== Single compile-time switch for VEGAS map usage =====
#define USE_VEGAS_MAP 1

namespace ORDER_UPDATE_MCMC {

    // =============================
    // Configuration & Constants
    // =============================

    // Dimensions to learn histograms for VEGAS map (at order 4)
    inline Int1DVector VEGAS_LEARN_DIMS = {6, 7, 8};

    // VEGAS map learning configuration
    constexpr int  VEGAS_LEARN_ORDER = 4;
    constexpr unsigned VEGAS_LEARN_BITS = 21;  // 2^20 learning samples
    constexpr int VEGAS_LEARN_SEGMENTS = 4;

    // Proposal histogram dimensions (for relearned histograms at N_MIN_ORDER)
    constexpr int HIST_LEARN_DIR_A = 6;
    constexpr int HIST_LEARN_DIR_B = 7;
    constexpr int HIST_LEARN_DIR_C = 8;

    constexpr int  HIST_NUM_BINS   = 50;
    constexpr unsigned LEARN_BITS  = 20;  // 2^20 learning samples for proposal histograms
    constexpr int      LEARN_SEGMENTS = 4;

    constexpr Real P_WITHIN_GAUSS   = 0.4;
    constexpr Real P_WITHIN_UNIFORM = 0.2;
    constexpr Real P_WITHIN_HIST    = 0.4;

    constexpr Real P_WITHIN = 0.5;
    constexpr Real P_ADD    = 0.34;
    constexpr Real P_REM    = 0.16;

    constexpr Real EPSILON_THRESHOLD = 5e-15;

    // =============================
    // Types
    // =============================

    using RealFunc     = std::function<Real(const Real*, int)>;

    // =============================
    // Utility
    // =============================

    inline Real sign_func(Real x) { return (x > 0) ? 1.0 : ((x < 0) ? -1.0 : 0.0); }

    inline Real Uniform01(std::mt19937& gen) {
        static thread_local std::uniform_real_distribution<Real> U(0.0, 1.0);
        return U(gen);
    }

    inline int RandomInt(int min, int max, std::mt19937& gen) {
        std::uniform_int_distribution<int> dist(min, max);
        return dist(gen);
    }

    // =============================
    // Histogram
    // =============================

    struct Histogram1D {
        int B{HIST_NUM_BINS};
        Real bin_width{1.0 / Real(HIST_NUM_BINS)};
        std::vector<unsigned long long> counts;
        std::vector<Real> probs;
        bool finalized{false};

        Histogram1D() : counts(B, 0ULL), probs(B, 0.0) {}
        explicit Histogram1D(int bins) : B(bins), bin_width(1.0/Real(bins)),
                                         counts(bins, 0ULL), probs(bins, 0.0) {}

        inline int bin_index(Real y) const {
            if (y <= 0.0) return 0;
            if (y >= 1.0) return B - 1;
            int idx = static_cast<int>(y / bin_width);
            if (idx < 0) idx = 0;
            if (idx >= B) idx = B - 1;
            return idx;
        }

        void observe(Real y) { counts[bin_index(y)]++; }

        void finalize() {
            unsigned long long T = 0ULL;
            for (auto c : counts) T += c;
            const Real add = 0.01 * static_cast<Real>(T);
            Real sum_plus = 0.0;
            probs.assign(B, 0.0);
            for (int b = 0; b < B; ++b) {
                Real cp = static_cast<Real>(counts[b]) + add;
                probs[b] = cp;
                sum_plus += cp;
            }
            if (sum_plus <= 0.0) {
                for (int b = 0; b < B; ++b) probs[b] = 1.0 / Real(B);
            } else {
                for (int b = 0; b < B; ++b) probs[b] /= sum_plus;
            }
            finalized = true;
        }

        Real density(Real y) const {
            if (!finalized) return 1.0;
            int b = bin_index(y);
            return probs[b] * Real(B);
        }

        Real sample(std::mt19937& gen) const {
            std::uniform_real_distribution<Real> U(0.0, 1.0);
            if (!finalized) return U(gen);
            Real u = U(gen), cum = 0.0;
            int b = 0;
            for (; b < B; ++b) { cum += probs[b]; if (u <= cum) break; }
            if (b >= B) b = B - 1;
            Real left = b * bin_width;
            return left + U(gen) * bin_width;
        }
    };

    // =============================
    // State
    // =============================

    struct VState {
        int n;
        Real1DVector Y; // size = 3*n
    };

    inline int dims_from_order(int n) { return 3 * n; }
    inline int first_dim_of_vertex(int v) { return 3 * v; }

    // =============================
    // VEGAS Maps Container
    // =============================

    struct VegasMapsContainer {
        std::map<int, VEGAS::MapData> maps;  // order -> map
        std::vector<Histogram1D> base_histograms;  // The histograms learned at order 4

        bool has_map(int order) const {
            return maps.find(order) != maps.end();
        }

        VEGAS::MapData* get_map(int order) {
            auto it = maps.find(order);
            return (it != maps.end()) ? &(it->second) : nullptr;
        }
    };

    // =============================
    // Function evaluation with VEGAS map
    // =============================

    inline Real eval_f_n(const RealFunc& f,
                         const Real1DVector& Y,
                         int n,
                         VegasMapsContainer* vmc = nullptr) {
        if (!vmc || !vmc->has_map(n)) {
            return f(Y.data(), n);
        }

        // Get the appropriate map for this order
        VEGAS::MapData* vm = vmc->get_map(n);
        if (!vm) {
            return f(Y.data(), n);
        }

        // Apply VEGAS transformation
        auto x = VEGAS::GetX(*vm, Y);
        Real v = f(x.data(), n);
        if (v == 0) return 0;

        Real lj = VEGAS::GetLogJac(*vm, Y);  // base-10
        Real L = std::log10(std::fabs(v)) + lj;
        if (L > 50.) throw std::runtime_error("func*Jac > 1e50");

        return (v > 0 ? 1 : -1) * std::pow(Real(10), L);
    }

    // --- Helper to evaluate abs(|f|) (clamped) and sign in ONE call ---
    inline void eval_abs_and_sign(const RealFunc& f,
                                  const Real1DVector& Y,
                                  int n,
                                  VegasMapsContainer* vmc,
                                  Real& abs_w,
                                  Real& sgn) {
        Real raw = eval_f_n(f, Y, n, vmc);
        abs_w = std::fabs(raw);
        if (abs_w < EPSILON_THRESHOLD) abs_w = EPSILON_THRESHOLD;
        sgn = sign_func(raw);
    }

    // =============================
    // VEGAS map learning at order 4
    // =============================

    inline std::vector<Histogram1D> learn_vegas_histograms_order4(
        const RealFunc& f,
        std::mt19937& gen,
        Real initial_sigma = 0.25,
        const Real1DVector* Y0_order4 = nullptr) {

        const unsigned long long N_learn = (1ULL << VEGAS_LEARN_BITS);
        const unsigned long long seg_len = N_learn / VEGAS_LEARN_SEGMENTS;

        VState S;
        S.n = VEGAS_LEARN_ORDER;
        S.Y.assign(dims_from_order(S.n), 0.0);
        if (Y0_order4 && static_cast<int>(Y0_order4->size()) == dims_from_order(VEGAS_LEARN_ORDER)) {
            S.Y = *Y0_order4;
        } else {
            for (auto& y : S.Y) y = Uniform01(gen);
        }

        std::vector<Histogram1D> H;
        H.reserve(VEGAS_LEARN_DIMS.size());
        for (size_t j = 0; j < VEGAS_LEARN_DIMS.size(); ++j) {
            H.emplace_back(HIST_NUM_BINS);
        }

        Real fval0 = std::fabs(f(S.Y.data(), S.n));
        if (fval0 < EPSILON_THRESHOLD) fval0 = EPSILON_THRESHOLD;

        for (int seg = 0; seg < VEGAS_LEARN_SEGMENTS; ++seg) {
            Real sigma = initial_sigma;
            std::normal_distribution<Real> ND(0.0, sigma);

            for (unsigned long long t = 0; t < seg_len; ++t) {
                int k = RandomInt(0, static_cast<int>(S.Y.size()) - 1, gen);

                Real1DVector Yp = S.Y;
                Real prop = std::fmod(S.Y[k] + ND(gen), Real(1.0));
                if (prop < Real(0.0)) prop += Real(1.0);
                Yp[k] = prop;

                Real f_new = std::fabs(f(Yp.data(), S.n));
                if (f_new < EPSILON_THRESHOLD) f_new = EPSILON_THRESHOLD;

                Real ratio = f_new / fval0;
                if (ratio > Uniform01(gen)) {
                    S.Y.swap(Yp);
                    fval0 = f_new;
                }

                for (size_t j = 0; j < VEGAS_LEARN_DIMS.size(); ++j) {
                    int d = VEGAS_LEARN_DIMS[j];
                    if (d >= 0 && d < static_cast<int>(S.Y.size())) {
                        H[j].observe(S.Y[d]);
                    }
                }
            }
        }

        for (auto& h : H) h.finalize();
        return H;
    }

    // =============================
    // Create VEGAS map for a specific order
    // =============================

    inline VEGAS::MapData make_vegas_map_for_order(int order, const std::vector<Histogram1D>& base_hists) {
        int DIM = dims_from_order(order);
        std::vector<Real> mins(DIM, 0.0), maxs(DIM, 1.0);
        auto vm = VEGAS::InitializeMap(DIM, HIST_NUM_BINS, mins.data(), maxs.data());

        const int K = static_cast<int>(base_hists.size());
        if (K == 0) return vm;

        // Apply stretching based on dimension index modulo 3
        // Args 0, 1, 2 are uniform (not stretched)
        // Args i (for i >= 3) are stretched based on histogram[i % 3]
        for (int d = 0; d < DIM; ++d) {
            if (d < 3) {
                // Keep uniform for first 3 dimensions
                continue;
            }
            int hist_idx = d % 3;
            if (hist_idx < K) {
                VEGAS::BuildFromHistogram1D(vm, d, base_hists[hist_idx].probs);
            }
        }

        return vm;
    }

    // =============================
    // Create all VEGAS maps for orders N_MIN to N_MAX
    // =============================

    inline VegasMapsContainer create_all_vegas_maps(const std::vector<Histogram1D>& base_hists) {
        VegasMapsContainer vmc;
        vmc.base_histograms = base_hists;

        // Create a map for each possible order
        for (int order = N_MIN_ORDER; order <= N_MAX_ORDER; ++order) {
            vmc.maps[order] = make_vegas_map_for_order(order, base_hists);
        }

        return vmc;
    }

    // =============================
    // Learn proposal histograms with stretched function
    // =============================

    inline std::array<Histogram1D, 3> learn_proposal_histograms(
        const RealFunc& f,
        std::mt19937& gen,
        int order,
        VegasMapsContainer* vmc = nullptr,
        Real initial_sigma = 0.25,
        const Real1DVector* Y0 = nullptr) {

        const unsigned long long N_learn = (1ULL << LEARN_BITS);
        const unsigned long long seg_len = N_learn / LEARN_SEGMENTS;

        VState S;
        S.n = order;
        S.Y.assign(dims_from_order(S.n), 0.0);
        if (Y0 && static_cast<int>(Y0->size()) == dims_from_order(order)) {
            S.Y = *Y0;
        } else {
            for (auto& y : S.Y) y = Uniform01(gen);
        }

        std::array<Histogram1D, 3> H;
        Real fval0 = std::fabs(eval_f_n(f, S.Y, S.n, vmc));
        if (fval0 < EPSILON_THRESHOLD) fval0 = EPSILON_THRESHOLD;

        for (int seg = 0; seg < LEARN_SEGMENTS; ++seg) {
            Real sigma = initial_sigma;
            std::normal_distribution<Real> ND(0.0, sigma);

            for (unsigned long long t = 0; t < seg_len; ++t) {
                int k = RandomInt(0, static_cast<int>(S.Y.size()) - 1, gen);

                Real1DVector Yp = S.Y;
                Real prop = std::fmod(S.Y[k] + ND(gen), Real(1.0));
                if (prop < Real(0.0)) prop += Real(1.0);
                Yp[k] = prop;

                Real f_new = std::fabs(eval_f_n(f, Yp, S.n, vmc));
                if (f_new < EPSILON_THRESHOLD) f_new = EPSILON_THRESHOLD;

                Real ratio = f_new / fval0;
                if (ratio > Uniform01(gen)) {
                    S.Y.swap(Yp);
                    fval0 = f_new;
                }

                const int dirs[3] = {HIST_LEARN_DIR_A, HIST_LEARN_DIR_B, HIST_LEARN_DIR_C};
                for (int j = 0; j < 3; ++j) {
                    int d = dirs[j];
                    if (d >= 0 && d < static_cast<int>(S.Y.size())) {
                        H[j].observe(S.Y[d]);
                    }
                }
            }
        }

        for (auto& h : H) h.finalize();
        return H;
    }

    // =============================
    // Seed finder
    // =============================

    inline bool find_good_starting_point_eval(Real1DVector& Y_start,
                                              Real& w_start,
                                              Real& o_start,
                                              const RealFunc& func,
                                              Real sigma,
                                              int order_num,
                                              std::mt19937& gen,
                                              VegasMapsContainer* vmc = nullptr,
                                              int max_attempts = 200,
                                              int test_steps = 100,
                                              Real min_accept_ratio = 0.03) {
        const int DIM = dims_from_order(order_num);
        Y_start.assign(DIM, 0.0);

        for (int attempt = 0; attempt < max_attempts; ++attempt) {
            for (int i = 0; i < DIM; ++i) Y_start[i] = Uniform01(gen);

            Real func_val = eval_f_n(func, Y_start, order_num, vmc);
            w_start = std::fabs(func_val);
            if (w_start < EPSILON_THRESHOLD) w_start = EPSILON_THRESHOLD;
            o_start = sign_func(func_val);

            VState S{order_num, Y_start};

            int accept_count = 0;
            for (int step = 0; step < test_steps; ++step) {
                const int D = static_cast<int>(S.Y.size());
                int k = RandomInt(0, D - 1, gen);

                Real1DVector Yp = S.Y;
                if (Uniform01(gen) < Real(0.7)) {
                    std::normal_distribution<Real> ND(0.0, sigma);
                    Real prop = std::fmod(S.Y[k] + ND(gen), Real(1.0));
                    if (prop < Real(0.0)) prop += Real(1.0);
                    Yp[k] = prop;
                } else {
                    Yp[k] = Uniform01(gen);
                }

                Real w_old = std::fabs(eval_f_n(func, S.Y, S.n, vmc));
                Real w_new = std::fabs(eval_f_n(func, Yp, S.n, vmc));
                if (w_old < EPSILON_THRESHOLD) w_old = EPSILON_THRESHOLD;
                if (w_new < EPSILON_THRESHOLD) w_new = EPSILON_THRESHOLD;

                Real ratio = (w_new / w_old);
                if (ratio > Uniform01(gen)) {
                    S.Y.swap(Yp);
                    accept_count++;
                }
            }

            Real test_accept_ratio = static_cast<Real>(accept_count) / static_cast<Real>(test_steps);
            if (test_accept_ratio > min_accept_ratio) return true;
        }

        return false;
    }

    // =============================
    // Proposals
    // =============================

    enum class Kernel { GAUSS, UNIF, HIST };

    struct StepResult {
        VState state;
        Real abs_weight;
        Real sign_val;
        bool accepted;
        int  last_direction;
    };

    struct StepProbs{
        Real p_add  = P_ADD;
        Real p_rem  = P_REM;
        Real p_with = P_WITHIN;
    };

    inline void clamp_probs_at_boundaries(StepProbs& P, int n){
        if(n<=N_MIN_ORDER){
            P.p_rem = 0.0;
            P.p_with = 1.0 - P.p_add;
        }
        if(n>=N_MAX_ORDER){
            P.p_add = 0.0;
            P.p_with = 1.0 - P.p_rem;
        }
    }

    // NOTE: All proposals now receive the *current state's* (w_old, sign_old)
    // and never recompute it. They evaluate the candidate ONCE and reuse it.

    inline StepResult propose_within(const VState& S,
                                     const RealFunc& f,
                                     std::mt19937& gen,
                                     const std::array<Histogram1D, 3>* learned_hists,
                                     bool use_histograms,
                                     Real sigma,
                                     VegasMapsContainer* vmc,
                                     Real w_old,       // passed-in current abs weight
                                     Real sign_old) {  // passed-in current sign
        StepResult R{S, w_old, sign_old, false, -1};

        const int D = static_cast<int>(S.Y.size());
        int k = RandomInt(0, D - 1, gen);
        R.last_direction = k;

        Kernel ker;
        Real u = Uniform01(gen);
        if (!use_histograms || learned_hists == nullptr) {
            ker = (u < Real(0.7)) ? Kernel::GAUSS : Kernel::UNIF;
        } else {
            if (u < P_WITHIN_GAUSS) ker = Kernel::GAUSS;
            else if (u < P_WITHIN_GAUSS + P_WITHIN_UNIFORM) ker = Kernel::UNIF;
            else ker = Kernel::HIST;
        }

        Real1DVector Yp = S.Y;
        if (ker == Kernel::GAUSS) {
            std::normal_distribution<Real> ND(0.0, sigma);
            Real prop = std::fmod(S.Y[k] + ND(gen), Real(1.0));
            if (prop < Real(0.0)) prop += Real(1.0);
            Yp[k] = prop;
        } else if (ker == Kernel::UNIF) {
            Yp[k] = Uniform01(gen);
        } else {
            const auto& hist = (*learned_hists)[k % 3];
            Yp[k] = hist.sample(gen);
        }

        Real w_new, sign_new;
        eval_abs_and_sign(f, Yp, S.n, vmc, w_new, sign_new);

        Real hast = 1.0;
        if (learned_hists != nullptr && (ker == Kernel::HIST)) {
            const auto& hist = (*learned_hists)[k % 3];
            Real h_x = hist.density(S.Y[k]);
            Real h_y = hist.density(Yp[k]);
            if (h_y < EPSILON_THRESHOLD) {
                // reject without recomputing old state
                return R;
            }
            hast = (h_x < EPSILON_THRESHOLD) ? 0.0 : (h_x / h_y);
        }

        Real ratio = (w_new / w_old) * hast;

        if (ratio > Uniform01(gen)) {
            R.state.Y.swap(Yp);
            R.abs_weight = w_new;
            R.sign_val   = sign_new;
            R.accepted   = true;
        }
        // else: keep R as old (already set)

        return R;
    }

    inline StepResult propose_add_at(const VState& S,
                                     const RealFunc& f,
                                     std::mt19937& gen,
                                     const std::array<Histogram1D,3>* learned_hists,
                                     int insert_vertex,
                                     Real /*sigma*/,
                                     VegasMapsContainer* vmc,
                                     Real w_old, Real sign_old) {
        StepResult R{S, w_old, sign_old, false, -1};
        if(S.n >= N_MAX_ORDER){
            return R; // cannot add; keep old evals
        }

        // sample new vertex
        Real y_new[3];
        Real h_density_prod=1.0;
        for(int j=0;j<3;++j){
            y_new[j]=(*learned_hists)[j].sample(gen);
            h_density_prod *= (*learned_hists)[j].density(y_new[j]);
        }

        const int n=S.n;
        int v = std::max(0,std::min(insert_vertex,n));
        int insert_pos = first_dim_of_vertex(v);

        VState Sp;
        Sp.n=n+1;
        Sp.Y.reserve(S.Y.size()+3);
        Sp.Y.insert(Sp.Y.end(),S.Y.begin(),S.Y.begin()+insert_pos);
        Sp.Y.insert(Sp.Y.end(),y_new,y_new+3);
        Sp.Y.insert(Sp.Y.end(),S.Y.begin()+insert_pos,S.Y.end());

        Real w_new, sign_new;
        eval_abs_and_sign(f, Sp.Y, Sp.n, vmc, w_new, sign_new);

        StepProbs P_old; clamp_probs_at_boundaries(P_old,n);
        StepProbs P_new; clamp_probs_at_boundaries(P_new,Sp.n);

        Real q_fwd=P_old.p_add*h_density_prod;
        Real q_rev=P_new.p_rem;

        Real ratio=(w_new/w_old)*(q_rev/q_fwd);

        if (ratio > Uniform01(gen)) {
            R.state=std::move(Sp);
            R.abs_weight=w_new;
            R.sign_val=sign_new;
            R.accepted=true;
        }
        // else keep old (already in R)

        return R;
    }

    inline StepResult propose_remove_at(const VState& S,
                                        const RealFunc& f,
                                        std::mt19937& gen,
                                        const std::array<Histogram1D,3>* learned_hists,
                                        int rem_vertex,
                                        Real /*sigma*/,
                                        VegasMapsContainer* vmc,
                                        Real w_old, Real sign_old) {
        StepResult R{S, w_old, sign_old, false, -1};
        if(S.n<=N_MIN_ORDER){
            return R; // cannot remove; keep old evals
        }

        const int n=S.n;
        int v=std::max(0,std::min(rem_vertex,n-1));
        int i0=first_dim_of_vertex(v);

        VState Sp;
        Sp.n=n-1;
        Sp.Y.reserve(S.Y.size()-3);
        Sp.Y.insert(Sp.Y.end(),S.Y.begin(),S.Y.begin()+i0);
        Sp.Y.insert(Sp.Y.end(),S.Y.begin()+i0+3,S.Y.end());

        Real w_new, sign_new;
        eval_abs_and_sign(f, Sp.Y, Sp.n, vmc, w_new, sign_new);

        StepProbs P_old; clamp_probs_at_boundaries(P_old,n);
        StepProbs P_new; clamp_probs_at_boundaries(P_new,Sp.n);

        Real q_fwd=P_old.p_rem;

        Real h_prod=1.0;
        h_prod *= (*learned_hists)[0].density(S.Y[i0+0]);
        h_prod *= (*learned_hists)[1].density(S.Y[i0+1]);
        h_prod *= (*learned_hists)[2].density(S.Y[i0+2]);

        Real q_rev=P_new.p_add*h_prod;

        Real ratio=(w_new/w_old)*(q_rev/q_fwd);

        if(ratio>Uniform01(gen)){
            R.state=std::move(Sp);
            R.abs_weight=w_new;
            R.sign_val=sign_new;
            R.accepted=true;
        }
        // else keep old

        return R;
    }

    // =============================
    // Main variable-order run (sign-only)
    // =============================

    struct RunSummary {
        Real sign_mean;
        Real accept_ratio;
        Real final_sigma;
        std::array<Histogram1D, 3> hists;
        std::vector<Real> dir_accept_ratio;

        std::vector<Real>                sign_mean_per_order;
        std::vector<unsigned long long>  sign_count_per_order;

        unsigned long long prod_attempt_within{0};
        unsigned long long prod_accept_within{0};
        unsigned long long prod_attempt_add{0};
        unsigned long long prod_accept_add{0};
        unsigned long long prod_attempt_remove{0};
        unsigned long long prod_accept_remove{0};

        unsigned long long warm_attempt_within{0};
        unsigned long long warm_accept_within{0};
        unsigned long long warm_attempt_add{0};
        unsigned long long warm_accept_add{0};
        unsigned long long warm_attempt_remove{0};
        unsigned long long warm_accept_remove{0};
    };

    inline RunSummary run_sign_variable_order(const RealFunc& f,
                                              Real initial_sigma,
                                              unsigned iter_bits,
                                              unsigned warm_bits,
                                              std::mt19937& gen) {

        auto t_start = std::chrono::high_resolution_clock::now();

        // Step 1: Create VEGAS maps for all orders if enabled
        VegasMapsContainer vegas_container;
        VegasMapsContainer* vmc_ptr = nullptr;

    #if USE_VEGAS_MAP
        std::cout << "Learning VEGAS histograms at order " << VEGAS_LEARN_ORDER << "...\n";

        // Learn histograms at order 4 for VEGAS map
        Real1DVector Y0_vegas(dims_from_order(VEGAS_LEARN_ORDER));
        for (auto& y : Y0_vegas) y = Uniform01(gen);

        auto vegas_hists = learn_vegas_histograms_order4(f, gen, 0.25, &Y0_vegas);

        // Create VEGAS maps for all orders
        vegas_container = create_all_vegas_maps(vegas_hists);
        vmc_ptr = &vegas_container;

        std::cout << "VEGAS maps created for orders " << N_MIN_ORDER << " to " << N_MAX_ORDER << ".\n";
    #endif

        // Step 2: Find good starting point with stretched function
        Real1DVector Y_seed;
        Real w_seed = 0.0, o_seed = 0.0;
        const int seed_order = N_MIN_ORDER;

        bool have_seed = find_good_starting_point_eval(
            Y_seed, w_seed, o_seed, f, initial_sigma, seed_order, gen, vmc_ptr,
            200, 100, 0.03
        );

        if (!have_seed) {
            Y_seed.assign(dims_from_order(seed_order), 0.0);
            for (auto& y : Y_seed) y = Uniform01(gen);
            eval_abs_and_sign(f, Y_seed, seed_order, vmc_ptr, w_seed, o_seed);
        }

        // Step 3: Learn proposal histograms with stretched function at N_MIN_ORDER
        std::cout << "Learning proposal histograms at order " << N_MIN_ORDER
                  << " with stretched function...\n";

        Real1DVector Y0_for_learn = Y_seed;
        if (static_cast<int>(Y0_for_learn.size()) != dims_from_order(N_MIN_ORDER)) {
            Y0_for_learn.assign(dims_from_order(N_MIN_ORDER), 0.0);
            for (auto &y : Y0_for_learn) y = Uniform01(gen);
        }

        auto learned = learn_proposal_histograms(f, gen, N_MIN_ORDER, vmc_ptr, 0.25, &Y0_for_learn);

        std::cout << "Proposal histograms learned.\n";

        // Step 4: Run MCMC with stretched function and proposal histograms
        const unsigned long long N_warm = (1ULL << warm_bits);
        const unsigned long long N_prod = (1ULL << iter_bits);

        VState S{seed_order, Y_seed};

        Real w_cur = 0.0, sign_cur = 0.0;
        eval_abs_and_sign(f, S.Y, S.n, vmc_ptr, w_cur, sign_cur);
        Real sigma = initial_sigma;

        unsigned long long warm_attempt_within = 0, warm_accept_within = 0;
        unsigned long long warm_attempt_add    = 0, warm_accept_add    = 0;
        unsigned long long warm_attempt_remove = 0, warm_accept_remove = 0;

        Real accept_warm = 0.0;
        const unsigned long long adapt_iv = std::max(100ULL, N_warm / 8);

        std::cout << "Running warm-up phase (" << N_warm << " steps)...\n";

        for (unsigned long long t = 0; t < N_warm; ++t) {
            StepProbs Pstep; clamp_probs_at_boundaries(Pstep, S.n);
            Real u = Uniform01(gen);
            StepResult R;

            if (u < Pstep.p_with) {
                ++warm_attempt_within;
                R = propose_within(S, f, gen, &learned, true, sigma, vmc_ptr, w_cur, sign_cur);
                if (R.accepted) ++warm_accept_within;
            } else if (u < Pstep.p_with + Pstep.p_add) {
                ++warm_attempt_add;
                int choice = RandomInt(0, S.n, gen);
                R = propose_add_at(S, f, gen, &learned, choice, sigma, vmc_ptr, w_cur, sign_cur);
                if (R.accepted) ++warm_accept_add;
            } else {
                ++warm_attempt_remove;
                int choice = RandomInt(0, S.n-1, gen);
                R = propose_remove_at(S, f, gen, &learned, choice, sigma, vmc_ptr, w_cur, sign_cur);
                if (R.accepted) ++warm_accept_remove;
            }

            if (R.accepted) {
                S = R.state;
                w_cur   = R.abs_weight;
                sign_cur= R.sign_val;
                accept_warm += 1.0;
            }

            if ((t + 1) % adapt_iv == 0) {
                Real a = accept_warm / Real(adapt_iv);
                if (a < 0.25)      sigma /= 1.6;
                else if (a > 0.65) sigma *= 1.6;
                accept_warm = 0.0;
            }
        }

        std::cout << "Running production phase (" << N_prod << " steps)...\n";

        Real sign_accum = 0.0;
        unsigned long long accepted = 0, attempted = 0;

        std::vector<Real> dir_accept_cnt(dims_from_order(N_MAX_ORDER), 0.0);
        std::vector<Real> dir_attempt_cnt(dims_from_order(N_MAX_ORDER), 0.0);

        const int ORD_VEC_SZ = N_MAX_ORDER + 1;
        std::vector<long double>          sign_sum_per_order(ORD_VEC_SZ, 0.0L);
        std::vector<unsigned long long>   sign_cnt_per_order(ORD_VEC_SZ, 0ULL);

        unsigned long long prod_attempt_within = 0, prod_accept_within = 0;
        unsigned long long prod_attempt_add    = 0, prod_accept_add    = 0;
        unsigned long long prod_attempt_remove = 0, prod_accept_remove = 0;

        for (unsigned long long t = 0; t < N_prod; ++t) {
            StepProbs Pstep; clamp_probs_at_boundaries(Pstep, S.n);
            Real u = Uniform01(gen);
            StepResult R;

            if (u < Pstep.p_with) {
                ++prod_attempt_within;
                R = propose_within(S, f, gen, &learned, true, sigma, vmc_ptr, w_cur, sign_cur);
                if (R.last_direction >= 0) {
                    int K = R.last_direction;
                    if (K < static_cast<int>(dir_attempt_cnt.size())) {
                        dir_attempt_cnt[K] += 1.0;
                        if (R.accepted) dir_accept_cnt[K] += 1.0;
                    }
                }
                if (R.accepted) ++prod_accept_within;
            } else if (u < Pstep.p_with + Pstep.p_add) {
                ++prod_attempt_add;
                int choice = RandomInt(0, S.n, gen);
                R = propose_add_at(S, f, gen, &learned, choice, sigma, vmc_ptr, w_cur, sign_cur);
                if (R.accepted) ++prod_accept_add;
            } else {
                ++prod_attempt_remove;
                int choice = RandomInt(0, S.n-1, gen);
                R = propose_remove_at(S, f, gen, &learned, choice, sigma, vmc_ptr, w_cur, sign_cur);
                if (R.accepted) ++prod_accept_remove;
            }

            attempted += 1ULL;
            if (R.accepted) {
                S = R.state;
                w_cur    = R.abs_weight;
                sign_cur = R.sign_val;
                accepted += 1ULL;
            }

            // Use the sign of the state used this step (R already contains correct sign)
            sign_accum += R.sign_val;
            int n_now = std::clamp(S.n, N_MIN_ORDER, N_MAX_ORDER);
            sign_sum_per_order[n_now] += static_cast<long double>(R.sign_val);
            sign_cnt_per_order[n_now] += 1ULL;
        }

        auto t_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<Real> dt = t_end - t_start;

        std::cout << "\nMC finished in " << dt.count() << " seconds.\n";

        RunSummary out;
        out.sign_mean    = sign_accum / Real(N_prod);
        out.accept_ratio = attempted > 0 ? (Real(accepted) / Real(attempted)) : 0.0;
        out.final_sigma  = sigma;
        out.hists        = learned;

        int Dcap = dims_from_order(N_MAX_ORDER);
        out.dir_accept_ratio.assign(Dcap, 0.0);
        for (int k = 0; k < Dcap; ++k) {
            if (dir_attempt_cnt[k] > 0.0) {
                out.dir_accept_ratio[k] = dir_accept_cnt[k] / dir_attempt_cnt[k];
            }
        }

        out.sign_mean_per_order.assign(ORD_VEC_SZ, 0.0);
        out.sign_count_per_order = sign_cnt_per_order;
        for (int n = 0; n < ORD_VEC_SZ; ++n) {
            if (N_prod > 0ULL) {
                out.sign_mean_per_order[n] =
                    static_cast<Real>(sign_sum_per_order[n] / static_cast<long double>(N_prod));
            } else {
                out.sign_mean_per_order[n] = 0.0;
            }
        }

        out.warm_attempt_within = warm_attempt_within;
        out.warm_accept_within  = warm_accept_within;
        out.warm_attempt_add    = warm_attempt_add;
        out.warm_accept_add     = warm_accept_add;
        out.warm_attempt_remove = warm_attempt_remove;
        out.warm_accept_remove  = warm_accept_remove;

        out.prod_attempt_within = prod_attempt_within;
        out.prod_accept_within  = prod_accept_within;
        out.prod_attempt_add    = prod_attempt_add;
        out.prod_accept_add     = prod_accept_add;
        out.prod_attempt_remove = prod_attempt_remove;
        out.prod_accept_remove  = prod_accept_remove;

        std::cout << "Sign mean: " << out.sign_mean << "\n";
        std::cout << "Accept ratio: " << out.accept_ratio << "\n";

        return out;
    }

} // namespace MC

#endif // MCMC_H
