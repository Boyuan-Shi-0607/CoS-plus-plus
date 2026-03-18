#ifndef MCMC_H
#define MCMC_H

#include <vector>
#include <random>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <thread>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <type_traits>
#include <sstream>
#include "online_stats.h"
#include "../utility/my_funcs.h"
#include "../utility/types.h"
#include <tuple>

namespace MC {

    // Uniform random number in [0, 1] for DoubleDouble
    inline DoubleDouble Uniform01_DD(std::mt19937& gen) {
        std::uniform_real_distribution<double> U(0.0, 1.0);
        return DoubleDouble(U(gen));
    }

    // Uniform integer in [min, max]
    inline int RandomInt(const int min, const int max, std::mt19937& gen) {
        std::uniform_int_distribution<int> dist(min, max);
        return dist(gen);
    }

    inline DoubleDouble sign_func_dd(const DoubleDouble& x) {
        if (x > DoubleDouble(0)) {
            return DoubleDouble(1);
        } else if (x < DoubleDouble(0)) {
            return DoubleDouble(-1);
        } else {
            return DoubleDouble(0);
        }
    }

    // Define epsilon threshold for different types
    template<typename T>
    constexpr T get_epsilon_threshold() {
        if constexpr (std::is_same_v<T, float>) {
            return static_cast<T>(1e-6);
        } else if constexpr (std::is_same_v<T, long double>) {
            return static_cast<T>(4e-18L);
        } else if constexpr (std::is_same_v<T, DoubleDouble>) {
            return DoubleDouble(1e-28);
        } else {  // double
            return static_cast<T>(5e-15);
        }
    }

    // Define a threshold for numerical stability
    const Real EPSILON_THRESHOLD = get_epsilon_threshold<Real>();
    const DoubleDouble EPSILON_THRESHOLD_DD = get_epsilon_threshold<DoubleDouble>();

    // MC_step_nd for DoubleDouble
    inline auto MC_step_nd_dd(
        const DD1DVector& X_old,
        DoubleDouble w_old,
        DoubleDouble o_1_old,
        DoubleDouble o_2_old,
        const DoubleDouble& sigma,
        const int order_num,
        const std::function<DoubleDouble(const DoubleDouble*, int)>& measure_func,
        const std::function<DoubleDouble(const DoubleDouble*, int)>& ref_func,
        const bool if_abs_ref,
        const bool if_abs_measure,
        std::mt19937& gen){

        auto X_new = X_old;
        const int ind = RandomInt(0, static_cast<int>(X_old.size()) - 1, gen);
        // std::normal_distribution<double> ND(0.0, sigma.to_double());
        // X_new[ind] = X_old[ind] + DoubleDouble(ND(gen));
        
        DoubleDouble sigma_in_use = sigma;
        if (ind < 3 * order_num) {
            if (ind % 3 == 1 || ind % 3 == 2) {
                sigma_in_use /= DoubleDouble(2.0);
            }
        }

        if (ind > 3 * order_num) {
            sigma_in_use *= DoubleDouble(2.2);
        }

        std::uniform_real_distribution<Real> coin(0.0, 1.0);
        if (coin(gen) < Real(0.7)) {
            // 70%: wrapped Gaussian step
            std::normal_distribution<Real> ND(0.0, sigma_in_use);
            X_new[ind] = DoubleDouble(std::fmod(X_old[ind].to_double() + ND(gen), Real(1.0)));
            if (X_new[ind] < DoubleDouble(0.0)) {
                X_new[ind] += DoubleDouble(1.0);
            }
        } else {
            // 30%: global uniform in [0,1)
            std::uniform_real_distribution<Real> U(Real(0.0), Real(1.0));
            X_new[ind] = DoubleDouble(U(gen));
        }

        DoubleDouble measure_val = measure_func(X_new.data(), order_num);
        DoubleDouble ref_val = ref_func(X_new.data(), order_num);
        DoubleDouble w_new = DD_ABS(measure_val) + DD_ABS(ref_val);

        if (if_abs_ref) {
            ref_val = DD_ABS(ref_val);
        }
        if (if_abs_measure) {
            measure_val = DD_ABS(measure_val);
        }

        bool accept;

        // Handle near-zero values for numerical stability
        const bool w_old_is_zero = (w_old < EPSILON_THRESHOLD_DD);
        const bool w_new_is_zero = (w_new < EPSILON_THRESHOLD_DD);

        DoubleDouble ratio;
        if (w_old_is_zero) {
            if (!w_new_is_zero) {
                ratio = DoubleDouble(1);
            } else {
                ratio = DoubleDouble(0);
            }
        } else {
            ratio = w_new / w_old;
        }

        if (ratio > Uniform01_DD(gen)) {
            DoubleDouble o_1_new, o_2_new;

            const bool ref_is_small = (DD_ABS(ref_val) < EPSILON_THRESHOLD_DD);
            const bool measure_is_small = (DD_ABS(measure_val) < EPSILON_THRESHOLD_DD);

            if (!ref_is_small && !measure_is_small) {
                o_1_new = ref_val / w_new;
                o_2_new = measure_val / w_new;
            }
            else if (ref_is_small && !measure_is_small) {
                o_1_new = DoubleDouble(0);
                o_2_new = DoubleDouble(1);
            }
            else if (!ref_is_small && measure_is_small) {
                o_1_new = DoubleDouble(1);
                o_2_new = DoubleDouble(0);
            }
            else {
                o_1_new = DoubleDouble(0.5);
                o_2_new = DoubleDouble(0.5);
            }

            accept = true;
            return std::make_tuple(X_new, w_new, o_1_new, o_2_new, accept);
        }
        accept = false;
        return std::make_tuple(X_old, w_old, o_1_old, o_2_old, accept);
    }

    // MC_step_sign for DoubleDouble
    inline auto MC_step_sign_dd(
        const DD1DVector& X_old,
        DoubleDouble w_old,
        DoubleDouble o_old,
        const DoubleDouble& sigma,
        const int order_num,
        const std::function<DoubleDouble(const DoubleDouble*, int)>& func,
        std::mt19937& gen){

        auto X_new = X_old;
        const int ind = RandomInt(0, static_cast<int>(X_old.size()) - 1, gen);
        // std::normal_distribution<double> ND(0.0, sigma.to_double());
        // X_new[ind] = X_old[ind] + DoubleDouble(ND(gen));
        // Use the same proposal as new MCMC.h

        DoubleDouble sigma_in_use = sigma;
        if (ind < 3 * order_num) {
            if (ind % 3 == 1 || ind % 3 == 2) {
                sigma_in_use /= DoubleDouble(2.0);
            }
        }

        if (ind > 3 * order_num) {
            sigma_in_use *= DoubleDouble(2.2);
        }

        std::uniform_real_distribution<Real> coin(0.0, 1.0);
        if (coin(gen) < Real(0.7)) {
            // 70%: wrapped Gaussian step
            std::normal_distribution<Real> ND(0.0, sigma_in_use);
            X_new[ind] = DoubleDouble(std::fmod(X_old[ind].to_double() + ND(gen), Real(1.0)));
            if (X_new[ind] < DoubleDouble(0.0)) {
                X_new[ind] += DoubleDouble(1.0);
            }
        } else {
            // 30%: global uniform in [0,1)
            std::uniform_real_distribution<Real> U(Real(0.0), Real(1.0));
            X_new[ind] = DoubleDouble(U(gen));
        }

        DoubleDouble func_val = func(X_new.data(), order_num);
        DoubleDouble w_new = DD_ABS(func_val);

        bool accept;

        const bool w_old_is_zero = (w_old < EPSILON_THRESHOLD_DD);
        const bool w_new_is_zero = (w_new < EPSILON_THRESHOLD_DD);

        DoubleDouble ratio;
        if (w_old_is_zero) {
            if (!w_new_is_zero) {
                ratio = DoubleDouble(1);
            } else {
                ratio = DoubleDouble(0);
            }
        } else {
            ratio = w_new / w_old;
        }

        if (ratio > Uniform01_DD(gen)) {
            DoubleDouble o_new = sign_func_dd(func_val);
            accept = true;
            return std::make_tuple(X_new, w_new, o_new, accept);
        }
        accept = false;
        return std::make_tuple(X_old, w_old, o_old, accept);
    }

    // Calculate variance of a vector of block means for DoubleDouble
    inline DoubleDouble calculate_block_variance_dd(const std::vector<DoubleDouble>& block_means) {
        if (block_means.size() <= 1) {
            return DoubleDouble(0);
        }

        DoubleDouble mean(0);
        for (const auto& val : block_means) {
            mean = mean + val;
        }
        mean = mean / DoubleDouble(static_cast<unsigned long long>(block_means.size()));

        DoubleDouble variance(0);
        for (const auto& val : block_means) {
            DoubleDouble diff = val - mean;
            variance = variance + diff * diff;
        }

        // Divide by N*(N-1) as requested
        variance = variance / DoubleDouble(static_cast<unsigned long long>(block_means.size()));
        variance = variance / DoubleDouble(static_cast<unsigned long long>(block_means.size()) - 1);

        return variance;
    }

    // run_sign for DoubleDouble with adaptive sigma - MODIFIED to return sigma and accept ratio
    inline std::tuple<DoubleDouble, std::vector<DoubleDouble>, DoubleDouble, DoubleDouble> run_sign_dd(
        const std::function<DoubleDouble(const DoubleDouble*, int)>& func,
        const DoubleDouble& initial_sigma,
        const int order_num,
        const int DIM,
        std::mt19937& gen,
        const unsigned int iter_bits,
        const unsigned int warm_bits,
        const unsigned long initial_seed = 0) {

        // Warmup and total steps
        const unsigned long long N_warm = (1ULL << warm_bits);
        auto N_steps = (1ULL << iter_bits) + N_warm;

        // Adaptive sigma
        DoubleDouble sigma = initial_sigma;

        // Calculate max blocks based on iter_bits
        const int max_blocks = iter_bits > 10 ? iter_bits - 10 + 1 : 1;
        OnlineStatsT<DoubleDouble> stats(max_blocks);

        // Create vectors for the chain's state and associated variables
        DD1DVector V_old(DIM);
        DoubleDouble w_old(0);
        DoubleDouble o_old(0);
        DoubleDouble tot(0);
        bool accept;

        std::uniform_real_distribution<double> dis(0.0, 1.0);

        // Initialize with a starting point in [0,1]^DIM
        for (int i = 0; i < DIM; ++i) {
            V_old[i] = DoubleDouble(dis(gen));
        }
        DoubleDouble func_old = func(V_old.data(), order_num);
        w_old = DD_ABS(func_old);
        o_old = sign_func_dd(func_old);

        DoubleDouble accept_counts(0);
        DoubleDouble warm_accept_counts(0);

        // Track sigma adaptations
        std::vector<std::pair<unsigned long long, DoubleDouble>> sigma_adaptations;

        // Check for warmup adaptation every this many steps (quarter of warmup)
        const unsigned long long adaptation_interval = std::max(N_warm / 4, 100ULL);

        for (unsigned long long step = 0; step < N_steps; ++step) {
            auto result = MC_step_sign_dd(V_old, w_old, o_old, sigma, order_num, func, gen);

            V_old = std::get<0>(result);
            w_old = std::get<1>(result);
            o_old = std::get<2>(result);
            accept = std::get<3>(result);

            if (step < N_warm) {
                // During warmup
                if (accept) {
                    warm_accept_counts = warm_accept_counts + DoubleDouble(1);
                }

                // Check and adapt sigma periodically during warmup
                if ((step + 1) % adaptation_interval == 0 && step > 0) {
                    DoubleDouble warm_accept_ratio = warm_accept_counts / DoubleDouble(adaptation_interval);

                    if (warm_accept_ratio < DoubleDouble(0.2)) {
                        sigma = sigma / DoubleDouble(1.5);
                        sigma_adaptations.emplace_back(step + 1, sigma);
                    } else if (warm_accept_ratio > DoubleDouble(0.7)) {
                        sigma = sigma * DoubleDouble(1.9);
                        sigma_adaptations.emplace_back(step + 1, sigma);
                    } else if (warm_accept_ratio > DoubleDouble(0.5)) {
                        sigma = sigma * DoubleDouble(1.4);
                        sigma_adaptations.emplace_back(step + 1, sigma);
                    }

                    // Reset warm acceptance counter
                    warm_accept_counts = DoubleDouble(0);
                }
            } else {
                // After warmup
                tot = tot + o_old;
                stats.update(o_old);
                if (accept) {
                    accept_counts = accept_counts + DoubleDouble(1);
                }
            }
        }

        tot = tot / DoubleDouble(N_steps - N_warm);
        DoubleDouble accept_ratio = accept_counts / DoubleDouble(N_steps - N_warm);

        // Calculate variances for each block size
        std::vector<DoubleDouble> block_variances;
        for (int i = 0; i < stats.getMaxBlocks(); ++i) {
            const auto& blockMeans = stats.getBlockMeans(i);
            if (!blockMeans.empty()) {
                block_variances.push_back(calculate_block_variance_dd(blockMeans));
            } else {
                block_variances.push_back(DoubleDouble(0)); // No data for this block size
            }
        }

        return std::make_tuple(tot, block_variances, sigma, accept_ratio);
    }

    // Multi-chain run for DoubleDouble with adaptive sigma - MODIFIED to return sigma and accept ratio
    inline std::tuple<DoubleDouble, std::vector<std::vector<DoubleDouble>>, std::vector<std::vector<DoubleDouble>>, DoubleDouble, std::vector<DoubleDouble>>
    run_dd(const std::vector<std::function<DoubleDouble(const DoubleDouble*, int)>>& measure_funcs,
           const std::vector<std::function<DoubleDouble(const DoubleDouble*, int)>>& ref_funcs,
           const std::vector<bool>& if_abs_ref,
           const std::vector<bool>& if_abs_measure,
           const DoubleDouble& initial_sigma,
           const int order_num,
           const int DIM,
           std::mt19937& gen,
           const unsigned int iter_bits,
           const unsigned int warm_bits,
           const unsigned int n_swap_rounds,
           const unsigned long initial_seed = 0) {

        // Warmup and total steps
        const unsigned long long N_warm = (1ULL << warm_bits);
        auto N_steps = (1ULL << iter_bits) + N_warm;

        // Adaptive sigma
        DoubleDouble sigma = initial_sigma;

        // Calculate max blocks based on iter_bits
        const int max_blocks = iter_bits > 10 ? iter_bits - 10 + 1 : 1;

        auto num_chains = measure_funcs.size();
        std::vector<OnlineStatsT<DoubleDouble>> stats_ref(num_chains, OnlineStatsT<DoubleDouble>(max_blocks));
        std::vector<OnlineStatsT<DoubleDouble>> stats_measure(num_chains, OnlineStatsT<DoubleDouble>(max_blocks));

        // Create vectors for the chains' states and associated variables
        std::vector<DD1DVector> V_old_chains(num_chains, DD1DVector(DIM));
        std::vector<DoubleDouble> w_old_chains(num_chains, DoubleDouble(0));
        std::vector<DoubleDouble> o_1_old_chains(num_chains, DoubleDouble(0));
        std::vector<DoubleDouble> o_2_old_chains(num_chains, DoubleDouble(0));
        std::vector<DoubleDouble> measure_tot(num_chains, DoubleDouble(0));
        std::vector<DoubleDouble> ref_tot(num_chains, DoubleDouble(0));
        std::vector<bool> accept(num_chains);

        std::uniform_real_distribution<double> dis(0.0, 1.0);

        // Initialize each chain with a starting point in [0,1]^DIM
        for (int c = 0; c < num_chains; ++c) {
            for (int i = 0; i < DIM; ++i) {
                V_old_chains[c][i] = DoubleDouble(dis(gen));
            }
            DoubleDouble measure_old = measure_funcs[c](V_old_chains[c].data(), order_num);
            DoubleDouble ref_old = ref_funcs[c](V_old_chains[c].data(), order_num);
            w_old_chains[c] = DD_ABS(measure_old) + DD_ABS(ref_old);

            // Case-by-case handling for initialization
            DoubleDouble ref_val_abs = if_abs_ref[c] ? DD_ABS(ref_old) : ref_old;
            DoubleDouble measure_val_abs = if_abs_measure[c] ? DD_ABS(measure_old) : measure_old;

            const bool ref_is_small = (DD_ABS(ref_val_abs) < EPSILON_THRESHOLD_DD);
            const bool measure_is_small = (DD_ABS(measure_val_abs) < EPSILON_THRESHOLD_DD);

            if (!ref_is_small && !measure_is_small) {
                // Both values are numerically significant
                o_1_old_chains[c] = ref_val_abs / w_old_chains[c];
                o_2_old_chains[c] = measure_val_abs / w_old_chains[c];
            }
            else if (ref_is_small && !measure_is_small) {
                // ref is effectively zero, measure is not
                o_1_old_chains[c] = DoubleDouble(0);
                o_2_old_chains[c] = DoubleDouble(1);
            }
            else if (!ref_is_small && measure_is_small) {
                // measure is effectively zero, ref is not
                o_1_old_chains[c] = DoubleDouble(1);
                o_2_old_chains[c] = DoubleDouble(0);
            }
            else {
                o_1_old_chains[c] = DoubleDouble(0.5);
                o_2_old_chains[c] = DoubleDouble(0.5);
            }
        }

        std::vector<DoubleDouble> accept_counts(num_chains, DoubleDouble(0));
        std::vector<DoubleDouble> warm_accept_counts(num_chains, DoubleDouble(0));

        // Track sigma adaptations
        std::vector<std::pair<unsigned long long, DoubleDouble>> sigma_adaptations;

        // Check for warmup adaptation every this many steps (quarter of warmup)
        const unsigned long long adaptation_interval = std::max(N_warm / 4, 100ULL);

        for (unsigned long long step = 0; step < N_steps; ++step) {
            // Update each chain using MC_step_nd_dd
            for (int c = 0; c < num_chains; ++c) {
                auto result = MC_step_nd_dd(V_old_chains[c], w_old_chains[c], o_1_old_chains[c],
                                         o_2_old_chains[c], sigma, order_num, measure_funcs[c],
                                         ref_funcs[c], if_abs_ref[c], if_abs_measure[c], gen);

                V_old_chains[c] = std::get<0>(result);
                w_old_chains[c] = std::get<1>(result);
                o_1_old_chains[c] = std::get<2>(result);
                o_2_old_chains[c] = std::get<3>(result);
                accept[c] = std::get<4>(result);
            }

            if (step < N_warm) {
                // During warmup
                for (int c = 0; c < num_chains; ++c) {
                    if (accept[c]) {
                        warm_accept_counts[c] = warm_accept_counts[c] + DoubleDouble(1);
                    }
                }

                // Check and adapt sigma periodically during warmup
                if ((step + 1) % adaptation_interval == 0 && step > 0) {
                    // Calculate average acceptance ratio across all chains
                    DoubleDouble total_warm_accept = DoubleDouble(0);
                    for (int c = 0; c < num_chains; ++c) {
                        total_warm_accept = total_warm_accept + warm_accept_counts[c];
                    }
                    DoubleDouble avg_warm_accept_ratio = total_warm_accept /
                        (DoubleDouble(num_chains) * DoubleDouble(adaptation_interval));

                    if (avg_warm_accept_ratio < DoubleDouble(0.2)) {
                        sigma = sigma / DoubleDouble(1.5);
                        sigma_adaptations.emplace_back(step + 1, sigma);
                    } else if (avg_warm_accept_ratio > DoubleDouble(0.7)) {
                        sigma = sigma * DoubleDouble(1.9);
                        sigma_adaptations.emplace_back(step + 1, sigma);
                    } else if (avg_warm_accept_ratio > DoubleDouble(0.5)) {
                        sigma = sigma * DoubleDouble(1.4);
                        sigma_adaptations.emplace_back(step + 1, sigma);
                    }

                    // Reset warm acceptance counters
                    for (int c = 0; c < num_chains; ++c) {
                        warm_accept_counts[c] = DoubleDouble(0);
                    }
                }
            } else {
                // After warmup
                for (int c = 0; c < num_chains; ++c) {
                    ref_tot[c] = ref_tot[c] + o_1_old_chains[c];
                    measure_tot[c] = measure_tot[c] + o_2_old_chains[c];
                    stats_ref[c].update(o_1_old_chains[c]);
                    stats_measure[c].update(o_2_old_chains[c]);
                    if (accept[c]) {
                        accept_counts[c] = accept_counts[c] + DoubleDouble(1);
                    }
                }
            }
        }

        // Calculate final accept ratios for each chain
        std::vector<DoubleDouble> accept_ratios(num_chains);
        for (int c = 0; c < num_chains; ++c) {
            ref_tot[c] = ref_tot[c] / DoubleDouble(N_steps - N_warm);
            measure_tot[c] = measure_tot[c] / DoubleDouble(N_steps - N_warm);
            accept_ratios[c] = accept_counts[c] / DoubleDouble(N_steps - N_warm);
        }

        DoubleDouble result(1);
        for (int c = 0; c < num_chains; ++c) {
            // Handle potential divide by zero in final result calculation
            if (DD_ABS(ref_tot[c]) < EPSILON_THRESHOLD_DD) {
                throw std::runtime_error(
                    "Warning: ref_tot[" + std::to_string(c) + "] is effectively zero, cannot compute ratio reliably"
                );
            } else {
                result = result * (measure_tot[c] / ref_tot[c]);
            }
        }

        // Calculate variances for each block size for each chain
        std::vector<std::vector<DoubleDouble>> ref_block_variances(num_chains);
        std::vector<std::vector<DoubleDouble>> measure_block_variances(num_chains);

        for (int c = 0; c < num_chains; ++c) {
            // Process reference function block means
            for (int i = 0; i < stats_ref[c].getMaxBlocks(); ++i) {
                const auto& blockMeans = stats_ref[c].getBlockMeans(i);
                if (!blockMeans.empty()) {
                    ref_block_variances[c].push_back(calculate_block_variance_dd(blockMeans));
                } else {
                    ref_block_variances[c].push_back(DoubleDouble(0));
                }
            }

            // Process measure function block means
            for (int i = 0; i < stats_measure[c].getMaxBlocks(); ++i) {
                const auto& blockMeans = stats_measure[c].getBlockMeans(i);
                if (!blockMeans.empty()) {
                    measure_block_variances[c].push_back(calculate_block_variance_dd(blockMeans));
                } else {
                    measure_block_variances[c].push_back(DoubleDouble(0));
                }
            }
        }

        return std::make_tuple(result, ref_block_variances, measure_block_variances, sigma, accept_ratios);
    }

    // Typedefs for convenience
    using OnlineStatsDD = OnlineStatsT<DoubleDouble>;

} // namespace MC

#endif //MCMC_H