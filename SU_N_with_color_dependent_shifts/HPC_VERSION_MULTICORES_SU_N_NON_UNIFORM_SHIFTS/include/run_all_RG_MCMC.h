#ifndef RUN_ALL_RG_MCMC_H
#define RUN_ALL_RG_MCMC_H

#include "RG_MCMC.h"
#include "parameters.h"
#include "../utility/types.h"
#include "integrand.h"
#include "dag_generator.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <type_traits>
#include <cctype>
#include <future>
#include <iomanip>
#include <future>
#include <filesystem>
#include <regex>
#include <map>
#include <set>
#include <thread>
#include <array>
#include <type_traits>
#include <iostream>
#include <random>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <cctype>
#include <future>
#include <iomanip>
#include <filesystem>
#include <regex>
#include <map>
#include <set>
#include <thread>
#include <array>

// Integrand function wrapper
inline Real integrand(const Real* args, const int order_num, const int i, const int j) {
    // Get pre-built DAG for this order (nullptr for order 1)
    const typed_dag::DAGData* dag = (order_num > 1 && dag_cache.find(order_num) != dag_cache.end())
                         ? &dag_cache[order_num] : nullptr;

    if (lnZ_mode) {
        return integrand_connected_typed(
            beta_vec[j], U, args,
            g_vec_flattened[std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            eps, cheby_degree,
            alpha_shift_vec,
            dag
        );
    }
    if (density_mode) {
        // For density mode, we need to calculate the difference between mu_h and mu_l
        Real result_mu_h = integrand_connected_typed(
            beta_vec[j], U, args,
            g_vec_mu_h_flattened[std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            eps, cheby_degree,
            alpha_shift_vec,
            dag
        );

        Real result_mu_l = integrand_connected_typed(
            beta_vec[j], U, args,
            g_vec_mu_l_flattened[std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            eps, cheby_degree,
            alpha_shift_vec,
            dag
        );

        return (result_mu_h - result_mu_l) / (beta_vec[j] * (mu_h_vec[0] - mu_l_vec[0]));
    }
    else {
        throw std::runtime_error("Invalid mode");
    }
}


// Helper function to select order based on job_index for each case
inline int selectOrderFromJobIndex(int set_num, unsigned long int job_index) {
    switch (set_num) {
        case 1: // Intervals of 40
            if (job_index <= 4) return 1;      // order 1
            else if (job_index <= 24) return 2; // order 2
            else if (job_index <= 42) return 3; // order 3
            else return 4;                       // order 4

        case 2: // Intervals of 50, 60, 50
            if (job_index <= 50) return 5;      // order 5
            else if (job_index <= 110) return 6; // order 6 (50 + 60)
            else return 7;                       // order 7

        case 3: // Intervals of 80
            if (job_index <= 80) return 8;      // order 8
            else return 9;                       // order 9

        default:
            throw std::runtime_error("Invalid set_num: " + std::to_string(set_num));
    }
}

inline void run_all_RG_MCMC(unsigned long int job_index) {
    int actual_iter_bits_thres = 25;
    int order_num;
    if (MCMC_mode == 0) {
        order_num = set_num; // For MCMC_mode = 0, use set_num as order directly
    } else {
        order_num = selectOrderFromJobIndex(set_num, job_index);
    }

    // Get parameters for this specific order
    auto sigma    = sigma_vec[order_num];
    int iter_bits = iter_bits_vec[order_num];
    int warm_bits = warm_bits_vec[order_num];
    int MC_bits   = MC_bits_vec[order_num];
    const int DIM = 3 * order_num;

    // Determine the number of chain transitions required
    const int num_grid_chains = N1_vec.size()   > 1 ? int(N1_vec.size())   - 1 : 0;
    const int num_beta_chains = beta_vec.size() > 1 ? int(beta_vec.size()) - 1 : 0;

    // Setup random number generation
    unsigned long int seed_start = seed_offset;
    for (int k = 1; k < order_num; ++k) {
        seed_start += seed_interval * 20000UL;
    }
    seed_start += (job_index - 1) * seed_interval;
    unsigned long int seed_end = seed_start + seed_interval;
    std::uniform_int_distribution<unsigned long int> dist_seed(seed_start, seed_end - 1);

    std::random_device rd;
    std::mt19937 gen_seed(rd());
    unsigned long int seed = dist_seed(gen_seed);
    std::mt19937 gen(seed);

    using namespace RG_MCMC;

    // ======================
    // Mode 0: Sign calculation
    // ======================
    if (MCMC_mode == 0) {
        // Save seed
        std::string seed_filename = std::string(save_path) + "/output_sign_seed_" + std::to_string(job_index) + ".txt";
        std::ofstream seed_file(seed_filename);
        if (seed_file.is_open()) {
            seed_file << seed << std::endl;
            seed_file.close();
        }

        auto sign_func_lambda = [](const Real* args, int order_num_param) {
            return integrand(args, order_num_param, 0, 0);
        };

        // Updated run_sign returns:
        //   sign_mean, block_variances, final_sigma, accept_ratio,
        //   direction_accept_ratios, tracked_values, histograms
        auto [grid_sign,
              block_variances,
              final_sigma,
              accept_ratio,
              direction_accept_ratios,
              tracked_data,
              histograms] =
            run_sign(sign_func_lambda, sigma, order_num, DIM, gen, iter_bits, warm_bits);

        // Save sign estimate
        {
            std::string filename = std::string(save_path) + "/output_sign_" + std::to_string(job_index) + ".txt";
            std::ofstream file(filename, std::ios::app);
            if (file.is_open()) {
                file << std::fixed << std::setprecision(15) << grid_sign << std::endl;
            }
        }

        // Save block std (sqrt of variance)
        {
            std::string block_std_filename = std::string(save_path) + "/output_sign_block_std_" + std::to_string(job_index) + ".txt";
            std::ofstream block_std_file(block_std_filename);
            if (block_std_file.is_open()) {
                block_std_file << std::fixed << std::setprecision(15);
                for (size_t i = 0; i < block_variances.size(); ++i) {
                    block_std_file << std::sqrt(block_variances[i]) << (i < block_variances.size() - 1 ? "," : "");
                }
                block_std_file << std::endl;
            }
        }

        // Save final sigma
        {
            std::string sigma_filename = std::string(save_path) + "/output_sign_sigma_" + std::to_string(job_index) + ".txt";
            std::ofstream sigma_file(sigma_filename);
            if (sigma_file.is_open()) {
                sigma_file << std::fixed << std::setprecision(15) << final_sigma << std::endl;
            }
        }

        // Save overall accept ratio
        {
            std::string accept_filename = std::string(save_path) + "/output_sign_accept_ratio_" + std::to_string(job_index) + ".txt";
            std::ofstream accept_file(accept_filename);
            if (accept_file.is_open()) {
                accept_file << std::fixed << std::setprecision(15) << accept_ratio << std::endl;
            }
        }

        // NEW: Save per-direction accept ratios for sign
        {
            std::string dir_accept_filename = std::string(save_path) + "/output_sign_direction_accept_ratios_" + std::to_string(job_index) + ".txt";
            std::ofstream dir_accept_file(dir_accept_filename);
            if (dir_accept_file.is_open()) {
                dir_accept_file << std::fixed << std::setprecision(15);
                for (size_t i = 0; i < direction_accept_ratios.size(); ++i) {
                    dir_accept_file << direction_accept_ratios[i];
                    if (i < direction_accept_ratios.size() - 1) dir_accept_file << ",";
                }
                dir_accept_file << std::endl;
                dir_accept_file.close();
            }
        }

        // NEW: Save tracked argument data for sign
        {
            std::string tracked_filename = std::string(save_path) + "/output_sign_tracked_arg_" + std::to_string(job_index) + ".txt";
            std::ofstream tracked_file(tracked_filename);
            if (tracked_file.is_open()) {
                tracked_file << std::fixed << std::setprecision(15);
                for (const auto &val : tracked_data) {
                    tracked_file << val << std::endl;
                }
                tracked_file.close();
            }
        }

        // NEW: Save warm-up histograms for sign
        for (size_t dir_idx = 0; dir_idx < histograms.size(); ++dir_idx) {
            std::string hist_filename = std::string(save_path) + "/output_sign_histogram_dir_" +
                                        std::to_string(dir_idx) + "_" + std::to_string(job_index) + ".txt";
            RG_MCMC::save_histogram(histograms[dir_idx], hist_filename);
        }
    }
    // ======================
    // Grid reduction modes (single transition)
    // ======================
    else if (MCMC_mode >= 1 && MCMC_mode <= num_grid_chains) {
        int chain_idx = MCMC_mode - 1;

        // Save seed
        {
            std::string seed_filename = std::string(save_path) + "/output_grid_mode_" + std::to_string(MCMC_mode) +
                                       "_order_" + std::to_string(order_num) + "_seed_" + std::to_string(job_index) + ".txt";
            std::ofstream seed_file(seed_filename);
            if (seed_file.is_open()) {
                seed_file << seed << std::endl;
                seed_file.close();
            }
        }

        // Measure/ref functions
        auto measure_func = [grid_idx = chain_idx](const Real* args, int order_num_param) {
            return integrand(args, order_num_param, grid_idx, 0);
        };
        auto ref_func = [grid_idx_next = chain_idx + 1](const Real* args, int order_num_param) {
            return integrand(args, order_num_param, grid_idx_next, 0);
        };

        int actual_iter_bits = iter_bits - 3;
        if (actual_iter_bits < actual_iter_bits_thres) actual_iter_bits = actual_iter_bits_thres;

        auto [result,
              ref_block_variances,
              measure_block_variances,
              final_sigma,
              accept_ratio,
              direction_accept_ratios,
              tracked_data,
              histograms] =
            run_with_maps(measure_func, ref_func,
                          sigma, order_num, DIM, gen,
                          actual_iter_bits, warm_bits);

        // Save histograms
        for (size_t dir_idx = 0; dir_idx < histograms.size(); ++dir_idx) {
            std::string hist_filename = std::string(save_path) + "/output_grid_mode_" + std::to_string(MCMC_mode) +
                                        "_order_" + std::to_string(order_num) + "_histogram_dir_" +
                                        std::to_string(dir_idx) + "_" + std::to_string(job_index) + ".txt";
            RG_MCMC::save_histogram(histograms[dir_idx], hist_filename);
        }

        // Save result
        {
            std::string filename = std::string(save_path) + "/output_grid_mode_" + std::to_string(MCMC_mode) +
                                   "_order_" + std::to_string(order_num) + "_" + std::to_string(job_index) + ".txt";
            std::ofstream file(filename, std::ios::app);
            if (file.is_open()) {
                file << std::fixed << std::setprecision(15) << result << std::endl;
            }
        }

        // Save ref-block std (kept same format "0,O1,<std>")
        {
            std::string block_std_filename = std::string(save_path) + "/output_grid_mode_" + std::to_string(MCMC_mode) +
                                             "_order_" + std::to_string(order_num) + "_block_std_" + std::to_string(job_index) + ".txt";
            std::ofstream block_std_file(block_std_filename);
            if (block_std_file.is_open()) {
                block_std_file << std::fixed << std::setprecision(15);
                for (size_t i = 0; i < ref_block_variances.size(); ++i) {
                    block_std_file << "0,O1," << std::sqrt(ref_block_variances[i]) << std::endl;
                }
            }
        }

        // Save final sigma
        {
            std::string sigma_filename = std::string(save_path) + "/output_grid_mode_" + std::to_string(MCMC_mode) +
                                         "_order_" + std::to_string(order_num) + "_sigma_" + std::to_string(job_index) + ".txt";
            std::ofstream sigma_file(sigma_filename);
            if (sigma_file.is_open()) {
                sigma_file << std::fixed << std::setprecision(15) << final_sigma << std::endl;
            }
        }

        // Save overall accept ratio
        {
            std::string accept_filename = std::string(save_path) + "/output_grid_mode_" + std::to_string(MCMC_mode) +
                                          "_order_" + std::to_string(order_num) + "_accept_ratio_" + std::to_string(job_index) + ".txt";
            std::ofstream accept_file(accept_filename);
            if (accept_file.is_open()) {
                accept_file << std::fixed << std::setprecision(15) << accept_ratio << std::endl;
            }
        }

        // Save per-direction accept ratios
        {
            std::string dir_accept_filename = std::string(save_path) + "/output_grid_mode_" + std::to_string(MCMC_mode) +
                                              "_order_" + std::to_string(order_num) + "_direction_accept_ratios_" + std::to_string(job_index) + ".txt";
            std::ofstream dir_accept_file(dir_accept_filename);
            if (dir_accept_file.is_open()) {
                dir_accept_file << std::fixed << std::setprecision(15);
                for (size_t i = 0; i < direction_accept_ratios.size(); ++i) {
                    dir_accept_file << direction_accept_ratios[i];
                    if (i < direction_accept_ratios.size() - 1) dir_accept_file << ",";
                }
                dir_accept_file << std::endl;
                dir_accept_file.close();
            }
        }

        // Save tracked argument data
        {
            std::string tracked_filename = std::string(save_path) + "/output_grid_mode_" + std::to_string(MCMC_mode) +
                                           "_order_" + std::to_string(order_num) + "_tracked_arg_" + std::to_string(job_index) + ".txt";
            std::ofstream tracked_file(tracked_filename);
            if (tracked_file.is_open()) {
                tracked_file << std::fixed << std::setprecision(15);
                for (const auto& val : tracked_data) {
                    tracked_file << val << std::endl;
                }
                tracked_file.close();
            }
        }
    }
    // ======================
    // Beta reduction modes (single transition)
    // ======================
    else if (MCMC_mode > num_grid_chains && MCMC_mode <= num_grid_chains + num_beta_chains) {
        int chain_idx = MCMC_mode - num_grid_chains - 1;
        int grid_idx  = num_grid_chains;

        // Save seed
        {
            std::string seed_filename = std::string(save_path) + "/output_beta_mode_" + std::to_string(MCMC_mode) +
                                        "_order_" + std::to_string(order_num) + "_seed_" + std::to_string(job_index) + ".txt";
            std::ofstream seed_file(seed_filename);
            if (seed_file.is_open()) {
                seed_file << seed << std::endl;
                seed_file.close();
            }
        }

        // Measure/ref functions
        auto measure_func = [=](const Real* args, int order_num_param) {
            return integrand(args, order_num_param, grid_idx, chain_idx);
        };
        auto ref_func = [=](const Real* args, int order_num_param) {
            return integrand(args, order_num_param, grid_idx, chain_idx + 1);
        };

        int actual_iter_bits = iter_bits - 3;
        if (actual_iter_bits < actual_iter_bits_thres) actual_iter_bits = actual_iter_bits_thres;

        auto [result,
              ref_block_variances,
              measure_block_variances,
              final_sigma,
              accept_ratio,
              direction_accept_ratios,
              tracked_data,
              histograms] =
            run_with_maps(measure_func, ref_func,
                          sigma, order_num, DIM, gen,
                          actual_iter_bits, warm_bits);

        // Save histograms
        for (size_t dir_idx = 0; dir_idx < histograms.size(); ++dir_idx) {
            std::string hist_filename = std::string(save_path) + "/output_beta_mode_" + std::to_string(MCMC_mode) +
                                        "_order_" + std::to_string(order_num) + "_histogram_dir_" +
                                        std::to_string(dir_idx) + "_" + std::to_string(job_index) + ".txt";
            RG_MCMC::save_histogram(histograms[dir_idx], hist_filename);
        }

        // Save result
        {
            std::string filename = std::string(save_path) + "/output_beta_mode_" + std::to_string(MCMC_mode) +
                                   "_order_" + std::to_string(order_num) + "_" + std::to_string(job_index) + ".txt";
            std::ofstream file(filename, std::ios::app);
            if (file.is_open()) {
                file << std::fixed << std::setprecision(15) << result << std::endl;
            }
        }

        // Save ref-block std
        {
            std::string block_std_filename = std::string(save_path) + "/output_beta_mode_" + std::to_string(MCMC_mode) +
                                             "_order_" + std::to_string(order_num) + "_block_std_" + std::to_string(job_index) + ".txt";
            std::ofstream block_std_file(block_std_filename);
            if (block_std_file.is_open()) {
                block_std_file << std::fixed << std::setprecision(15);
                for (size_t i = 0; i < ref_block_variances.size(); ++i) {
                    block_std_file << "0,O1," << std::sqrt(ref_block_variances[i]) << std::endl;
                }
            }
        }

        // Save final sigma
        {
            std::string sigma_filename = std::string(save_path) + "/output_beta_mode_" + std::to_string(MCMC_mode) +
                                         "_order_" + std::to_string(order_num) + "_sigma_" + std::to_string(job_index) + ".txt";
            std::ofstream sigma_file(sigma_filename);
            if (sigma_file.is_open()) {
                sigma_file << std::fixed << std::setprecision(15) << final_sigma << std::endl;
            }
        }

        // Save overall accept ratio
        {
            std::string accept_filename = std::string(save_path) + "/output_beta_mode_" + std::to_string(MCMC_mode) +
                                          "_order_" + std::to_string(order_num) + "_accept_ratio_" + std::to_string(job_index) + ".txt";
            std::ofstream accept_file(accept_filename);
            if (accept_file.is_open()) {
                accept_file << std::fixed << std::setprecision(15) << accept_ratio << std::endl;
            }
        }

        // Save per-direction accept ratios
        {
            std::string dir_accept_filename = std::string(save_path) + "/output_beta_mode_" + std::to_string(MCMC_mode) +
                                              "_order_" + std::to_string(order_num) + "_direction_accept_ratios_" + std::to_string(job_index) + ".txt";
            std::ofstream dir_accept_file(dir_accept_filename);
            if (dir_accept_file.is_open()) {
                dir_accept_file << std::fixed << std::setprecision(15);
                for (size_t i = 0; i < direction_accept_ratios.size(); ++i) {
                    dir_accept_file << direction_accept_ratios[i];
                    if (i < direction_accept_ratios.size() - 1) dir_accept_file << ",";
                }
                dir_accept_file << std::endl;
                dir_accept_file.close();
            }
        }

        // Save tracked argument data
        {
            std::string tracked_filename = std::string(save_path) + "/output_beta_mode_" + std::to_string(MCMC_mode) +
                                           "_order_" + std::to_string(order_num) + "_tracked_arg_" + std::to_string(job_index) + ".txt";
            std::ofstream tracked_file(tracked_filename);
            if (tracked_file.is_open()) {
                tracked_file << std::fixed << std::setprecision(15);
                for (const auto& val : tracked_data) {
                    tracked_file << val << std::endl;
                }
                tracked_file.close();
            }
        }
    }
    // ======================
    // Direct Monte Carlo integration (VEGAS-free)
    // ======================
    else if (MCMC_mode == num_grid_chains + num_beta_chains + 1) {
        // Save seed
        {
            std::string seed_filename = std::string(save_path) + "/output_MC_order_" + std::to_string(order_num) +
                                        "_seed_" + std::to_string(job_index) + ".txt";
            std::ofstream seed_file(seed_filename);
            if (seed_file.is_open()) {
                seed_file << seed << std::endl;
                seed_file.close();
            }
        }

        const int dim = 3 * order_num;
        unsigned long numSamples = 1UL << MC_bits;

        int largest_grid_idx = num_grid_chains;
        int largest_beta_idx = num_beta_chains;

        auto mc_func = [=](const Real* args, int order_num_param) {
            return std::abs(integrand(args, order_num_param, largest_grid_idx, largest_beta_idx));
        };

        // Plain Monte Carlo over [0,1]^dim (no VEGAS)
        Real sum = 0.0;
        Real1DVector x(dim);
        std::uniform_real_distribution<Real> dist(0.0, 1.0);

        for (unsigned long i = 0; i < numSamples; i++) {
            for (int j = 0; j < dim; j++) {
                x[j] = dist(gen);
            }
            sum += mc_func(x.data(), order_num);
        }

        Real integral = sum / Real(numSamples);

        std::string filename = std::string(save_path) + "/output_MC_order_" + std::to_string(order_num) +
                               "_" + std::to_string(job_index) + ".txt";
        std::ofstream file(filename, std::ios::app);
        if (file.is_open()) {
            file << std::fixed << std::setprecision(15) << integral << std::endl;
        }
    }
    else {
        throw std::runtime_error("Invalid MCMC_mode: " + std::to_string(MCMC_mode));
    }
}
#endif //RUN_ALL_RG_MCMC_H
