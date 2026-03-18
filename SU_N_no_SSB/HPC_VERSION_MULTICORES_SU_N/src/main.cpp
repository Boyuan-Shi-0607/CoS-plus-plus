#include "../utility/resize_and_initialize.h"
#include "MCMC.h"  // Include the VEGAS-integrated MCMC header
#include "integrand.h"
#include "load_data.h"
#include "../utility/types.h"
#include "../utility/check_range.h"
#include "lnZ_free_part.h"
#include "dag_generator.h"
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

// For CPUID on Linux
#include <cpuid.h>

// Global storage for pre-trained VEGAS maps
std::map<std::string, VEGAS::MapData> global_vegas_maps;
std::mutex vegas_map_mutex;

//================================================================================
// 1. CPU INFORMATION CHECK
//================================================================================

void cpuid(int function_id, std::array<int, 4>& registers) {
    __cpuid(function_id, registers[0], registers[1], registers[2], registers[3]);
}

void display_cpu_info() {
    std::cout << "--- CPU Information ---\n";
    std::array<int, 4> regs;
    std::string vendor_id;
    std::string brand_string;

    // Get Vendor ID
    cpuid(0, regs);
    vendor_id += std::string(reinterpret_cast<char*>(&regs[1]), 4);
    vendor_id += std::string(reinterpret_cast<char*>(&regs[3]), 4);
    vendor_id += std::string(reinterpret_cast<char*>(&regs[2]), 4);
    std::cout << "CPU Vendor: " << vendor_id << std::endl;

    // Get Brand String
    for (int i = 0x80000002; i <= 0x80000004; ++i) {
        cpuid(i, regs);
        brand_string += std::string(reinterpret_cast<char*>(&regs[0]), 4);
        brand_string += std::string(reinterpret_cast<char*>(&regs[1]), 4);
        brand_string += std::string(reinterpret_cast<char*>(&regs[2]), 4);
        brand_string += std::string(reinterpret_cast<char*>(&regs[3]), 4);
    }
    brand_string = brand_string.substr(0, brand_string.find('\0'));
    std::cout << "CPU Brand: " << brand_string << std::endl;

    unsigned int core_count = std::thread::hardware_concurrency();
    std::cout << "Logical Cores Detected: " << core_count << std::endl;
}

//================================================================================
// ORIGINAL FUNCTIONS
//================================================================================

// Convert string to Real
Real string_to_real(const std::string& s) {
    return std::stod(s);
}

// Integrand function wrapper
Real integrand(const Real* args, const int order_num, const int i, const int j) {
    // Get pre-built DAG for this order (nullptr for order 1)
    const DAGData* dag = (order_num > 1 && dag_cache.find(order_num) != dag_cache.end())
                         ? &dag_cache[order_num] : nullptr;

    if (lnZ_mode) {
        return integrand_connected(
            beta_vec[j], U, args,
            g_vec_flattened[std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            eps, cheby_degree,
            alpha_shift,
            dag
        );
    }
    if (density_mode) {
        // For density mode, we need to calculate the difference between mu_h and mu_l
        Real result_mu_h = integrand_connected(
            beta_vec[j], U, args,
            g_vec_mu_h_flattened[std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            eps, cheby_degree,
            alpha_shift,
            dag
        );

        Real result_mu_l = integrand_connected(
            beta_vec[j], U, args,
            g_vec_mu_l_flattened[std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            eps, cheby_degree,
            alpha_shift,
            dag
        );

        return (result_mu_h - result_mu_l) / (beta_vec[j] * (mu_h - mu_l));
    }
    else {
        throw std::runtime_error("Invalid mode");
    }
}

// Helper function to select order based on job_index for each case
int selectOrderFromJobIndex(int set_num, unsigned long int job_index) {
    switch (set_num) {
        case 1: // Intervals of 40
            if (job_index <= 8) return 1;      // order 1
            else if (job_index <= 24) return 2; // order 2
            else if (job_index <= 44) return 3; // order 3
            else return 4;                       // order 4

        case 2: // Intervals of 50, 60, 50
            if (job_index <= 20) return 5;      // order 5
            else if (job_index <= 40) return 6; // order 6 (50 + 60)
            else return 7;                       // order 7

        case 3: // Intervals of 80
            return 8;      // order 8

        default:
            throw std::runtime_error("Invalid set_num: " + std::to_string(set_num));
    }
}

//================================================================================
// MODIFIED MCMC RUN FUNCTION
//================================================================================

void run_all(unsigned long int job_index) {

    int actual_iter_bits_thres = 25;

    int order_num;
    if (MCMC_mode == 0) {
        order_num = set_num; // For MCMC_mode = 0, use set_num as order directly
    } else {
        order_num = selectOrderFromJobIndex(set_num, job_index);
    }

    // Get parameters for this specific order
    auto sigma = sigma_vec[order_num];
    int iter_bits = iter_bits_vec[order_num];
    int warm_bits = warm_bits_vec[order_num];
    int MC_bits = MC_bits_vec[order_num];
    const int DIM = 3 * order_num;

    // Determine the number of chain transitions required
    const int num_grid_chains = N1_vec.size() > 1 ? N1_vec.size() - 1 : 0;
    const int num_beta_chains = beta_vec.size() > 1 ? beta_vec.size() - 1 : 0;

    /*
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
    */
    
    std::random_device rd;
    std::mt19937 gen_seed(rd());
    std::uniform_int_distribution<unsigned long int> dist_seed(0, 10000000000000000ULL);

    unsigned long int seed = dist_seed(gen_seed);
    // *** IMPORTANT ***
    // Seed mt19937 via seed_seq with both halves of seed64
    std::seed_seq seq{
        static_cast<uint32_t>(seed),
        static_cast<uint32_t>(seed >> 32)
    };

    std::mt19937 gen(seq);

    using namespace MC;

    // Mode 0: Sign calculation
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
            MC::save_histogram(histograms[dir_idx], hist_filename);
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
            MC::save_histogram(histograms[dir_idx], hist_filename);
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
            MC::save_histogram(histograms[dir_idx], hist_filename);
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

//================================================================================
// 2. INTEGRAND SPEED TEST
//================================================================================

void run_integrand_speed_test_with_loaded_data() {
   std::cout << "\n--- Speed Test: integrand (using loaded data) ---\n";
   std::cout << "Testing for i=0, j=0\n";
   std::cout << " Order |      Avg Time/Run |   Runs/Sec" << std::endl;
   std::cout << "-------|-------------------|-----------" << std::endl;

   std::mt19937 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count());
   std::uniform_real_distribution<double> dist(0.0, 1.0);
   Real result_sink = 0.0;

   for (int order_num = 4; order_num <= 9; ++order_num) {
       const int dim = 3 * order_num;
       std::vector<Real> args(dim);
       long long runs = 0;
       const auto benchmark_duration = std::chrono::seconds(10);

       auto start_time = std::chrono::steady_clock::now();
       std::chrono::duration<double> elapsed_time;

       do {
           for (int k = 0; k < dim; ++k) {
               args[k] = dist(rng);
           }

           result_sink += integrand(args.data(), order_num, 0, 0);
           runs++;
           elapsed_time = std::chrono::steady_clock::now() - start_time;
       } while (elapsed_time < benchmark_duration);

       double total_seconds = elapsed_time.count();
       double time_per_run = total_seconds / runs;
       double runs_per_second = runs / total_seconds;

       std::cout << "result_sink: " << result_sink << std::endl;

       std::cout << std::setw(6) << order_num << " | ";

       if (time_per_run < 1e-3) {
           std::cout << std::fixed << std::setprecision(3) << std::setw(13) << time_per_run * 1e6 << " us | ";
       } else {
           std::cout << std::fixed << std::setprecision(3) << std::setw(13) << time_per_run * 1e3 << " ms | ";
       }

       std::cout << std::fixed << std::setprecision(1) << std::setw(9) << runs_per_second << std::endl;
   }
}

// Utility function to trim whitespace from both ends of a string
void trim(std::string &s) {
   s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
       return !std::isspace(ch);
   }));
   s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
       return !std::isspace(ch);
   }).base(), s.end());
}

bool stringToBool(const std::string &str) {
   std::string lowerStr = str;
   std::transform(lowerStr.begin(), lowerStr.end(), lowerStr.begin(),
                  [](unsigned char c){ return std::tolower(c); });

   if(lowerStr == "true" || lowerStr == "1")
       return true;
   else if(lowerStr == "false" || lowerStr == "0")
       return false;
   else
       throw std::invalid_argument("Invalid bool string: " + str);
}

//================================================================================
// MAIN FUNCTION
//================================================================================

int main() {
   // 1. Run CPU diagnostic check first
   display_cpu_info();

   // 2. Proceed with the original application logic to load data
   std::cout << "\n--- Starting Main Application Logic ---\n";

   check_range();

   // Initialize parameter vectors
   initialize_parameter_vectors();

    // Dictionary to hold parameters
    std::unordered_map<std::string, std::vector<std::string>> params;
    std::ifstream infile("../params.txt");

    if (!infile) {
        std::cerr << "Could not open param.txt" << std::endl;
        return 1;
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty()) continue;

        size_t pos = line.find('=');
        if (pos == std::string::npos)
            continue;

        std::string key = line.substr(0, pos);
        std::string valueStr = line.substr(pos + 1);
        trim(key);
        trim(valueStr);

        std::vector<std::string> tokens;
        std::stringstream ss(valueStr);
        std::string token;
        while (std::getline(ss, token, ',')) {
            trim(token);
            if (!token.empty())
                tokens.push_back(token);
        }

        params[key] = tokens;
    }

    // Load integer parameters
    for (const auto &entry : params["N1"]) {
        N1_vec.push_back(std::stoi(entry));
    }

    for (const auto &entry : params["N2"]) {
        N2_vec.push_back(std::stoi(entry));
    }

    // Load Real parameters
    for (const auto &entry : params["BETA"]) {
        beta_vec.push_back(string_to_real(entry));
    }

    cheby_degree = std::stoi(params["CHEBY_DEGREE"][0]);
    cores_cpp = std::stoi(params["CORES_CPP"][0]);
    eps = string_to_real(params["EPS"][0]);
    alpha_shift = string_to_real(params["ALPHA_SHIFT"][0]);
    mu = string_to_real(params["MU"][0]);
    mu_h = string_to_real(params["MU_h"][0]);
    mu_l = string_to_real(params["MU_l"][0]);

    U = string_to_real(params["U"][0]);

    lnZ_mode = stringToBool(params["LNZ_MODE"][0]);
    density_mode = stringToBool(params["DENSITY_MODE"][0]);

    // Read parameters for all orders 1-9
    for (int order = 1; order <= 9; ++order) {
        sigma_vec[order] = string_to_real(params["SIGMA_ORDER_" + std::to_string(order)][0]);
        iter_bits_vec[order] = std::stoi(params["ITER_BITS_ORDER_" + std::to_string(order)][0]);
        warm_bits_vec[order] = std::stoi(params["WARM_BITS_ORDER_" + std::to_string(order)][0]);
        MC_bits_vec[order] = std::stoi(params["MC_BITS_ORDER_" + std::to_string(order)][0]);
    }

    seed_offset = std::stoul(params["SEED_OFFSET"][0]);
    seed_interval = std::stoul(params["SEED_INTERVAL"][0]);
    MCMC_mode = std::stoi(params["MCMC_mode"][0]);
    // set_num = std::stoi(params["SET_NUM"][0]);
    MCMC_cut = string_to_real(params["MCMC_CUT"][0]);
    greens_func_exp_cut = string_to_real(params["GREENS_FUNC_EXP_CUT"][0]);

    // Load data
    std::cout << "Loading data...\n";
    load_data();

    // Build all DAGs once at startup
    initialize_all_dags();

    std::cout << "\n--- Loaded Parameters ---\n";
    std::cout << "Using Real (double) for all computations\n";

    // Print parameters
    std::cout << "cheby_degree: " << cheby_degree << "\n";
    std::cout << "cores_cpp: " << cores_cpp << "\n";
    std::cout << "eps: " << eps << "\n";
    std::cout << "mu: " << mu << "\n";
    std::cout << "U: " << U << "\n";
    std::cout << "N1: ";
    for (auto n : N1_vec) std::cout << n << " ";
    std::cout << "\nN2: ";
    for (auto n : N2_vec) std::cout << n << " ";
    std::cout << "\nbeta: ";
    for (auto b : beta_vec) std::cout << b << " ";
    std::cout << "\n";
    std::cout << "MCMC_mode: " << MCMC_mode << "\n";
    std::cout << "set_num: " << set_num << "\n";

    // 4. RUN THE SPEED TEST
    run_integrand_speed_test_with_loaded_data();

    // 5. Continue with the rest of your main application logic
    std::cout << "\n--- Continuing with MCMC simulation ---\n";

    int threads_per_job = cores_cpp;

    // Time the parallel execution of std::async
    auto start_time = std::chrono::high_resolution_clock::now();

    // Create a vector to store the futures
    std::vector<std::future<void>> futures;

    // Launch parallel threads
    for (int thread_num = 1; thread_num <= threads_per_job; ++thread_num) {
        unsigned long int run_index = thread_num;

        futures.push_back(std::async(std::launch::async, [run_index]() {
            run_all(run_index);
        }));
    }

    // Wait for all threads to complete and propagate exceptions
    for (auto& future : futures) {
        try {
            future.get();
        } catch (const std::exception& e) {
            std::cerr << "Thread exception: " << e.what() << std::endl;
            throw;
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end_time - start_time;
    std::cout << "Parallel execution time: " << duration.count() << " seconds\n";

    // Combine output files by order number
   // Wait 3 minutes to ensure all files are written
   std::cout << "Waiting 3 minutes for all files to be written..." << std::endl;
   std::this_thread::sleep_for(std::chrono::minutes(3));
   std::cout << "Starting file combination..." << std::endl;

    // Map to store files by order number with their job indices
    std::map<int, std::vector<std::pair<std::string, int>>> regular_files_by_order;  // path, job_index
    std::map<int, std::vector<std::pair<std::string, int>>> block_std_files_by_order;
    std::map<int, std::vector<std::pair<std::string, int>>> seed_files_by_order;
    std::map<int, std::vector<std::pair<std::string, int>>> sigma_files_by_order;
    std::map<int, std::vector<std::pair<std::string, int>>> accept_ratio_files_by_order;
    std::map<int, std::vector<std::pair<std::string, int>>> direction_accept_ratio_files_by_order;
    std::map<int, std::vector<std::pair<std::string, int>>> tracked_arg_files_by_order;
    std::map<int, std::vector<std::tuple<std::string, int, int>>> histogram_files_by_order; // NEW: path, job_idx, dir_idx

    // Separate vectors for sign files (no order number) with job indices
    std::vector<std::pair<std::string, int>> sign_files;
    std::vector<std::pair<std::string, int>> sign_block_std_files;
    std::vector<std::pair<std::string, int>> sign_seed_files;
    std::vector<std::pair<std::string, int>> sign_sigma_files;
    std::vector<std::pair<std::string, int>> sign_accept_ratio_files;
    std::vector<std::pair<std::string, int>> sign_direction_accept_ratio_files;
    std::vector<std::pair<std::string, int>> sign_tracked_arg_files;
    std::vector<std::tuple<std::string, int, int>> sign_histogram_files; // NEW: path, job_idx, dir_idx

    // Regular expressions to match file patterns and extract order numbers and job indices
    std::regex sign_pattern(R"(output_sign_(\d+)\.txt)");
    std::regex sign_block_pattern(R"(output_sign_block_std_(\d+)\.txt)");
    std::regex sign_seed_pattern(R"(output_sign_seed_(\d+)\.txt)");
    std::regex sign_sigma_pattern(R"(output_sign_sigma_(\d+)\.txt)");
    std::regex sign_accept_pattern(R"(output_sign_accept_ratio_(\d+)\.txt)");
    std::regex sign_dir_accept_pattern(R"(output_sign_direction_accept_ratios_(\d+)\.txt)");
    std::regex sign_tracked_arg_pattern(R"(output_sign_tracked_arg_(\d+)\.txt)");
    std::regex sign_hist_pattern(R"(output_sign_histogram_dir_(\d+)_(\d+)\.txt)"); // NEW

    std::regex grid_pattern(R"(output_grid_mode_\d+_order_(\d+)_(\d+)\.txt)");
    std::regex grid_block_pattern(R"(output_grid_mode_\d+_order_(\d+)_block_std_(\d+)\.txt)");
    std::regex grid_seed_pattern(R"(output_grid_mode_\d+_order_(\d+)_seed_(\d+)\.txt)");
    std::regex grid_sigma_pattern(R"(output_grid_mode_\d+_order_(\d+)_sigma_(\d+)\.txt)");
    std::regex grid_accept_pattern(R"(output_grid_mode_\d+_order_(\d+)_accept_ratio_(\d+)\.txt)");
    std::regex grid_dir_accept_pattern(R"(output_grid_mode_\d+_order_(\d+)_direction_accept_ratios_(\d+)\.txt)");
    std::regex grid_tracked_arg_pattern(R"(output_grid_mode_\d+_order_(\d+)_tracked_arg_(\d+)\.txt)");
    std::regex grid_hist_pattern(R"(output_grid_mode_\d+_order_(\d+)_histogram_dir_(\d+)_(\d+)\.txt)"); // NEW

    std::regex beta_pattern(R"(output_beta_mode_\d+_order_(\d+)_(\d+)\.txt)");
    std::regex beta_block_pattern(R"(output_beta_mode_\d+_order_(\d+)_block_std_(\d+)\.txt)");
    std::regex beta_seed_pattern(R"(output_beta_mode_\d+_order_(\d+)_seed_(\d+)\.txt)");
    std::regex beta_sigma_pattern(R"(output_beta_mode_\d+_order_(\d+)_sigma_(\d+)\.txt)");
    std::regex beta_accept_pattern(R"(output_beta_mode_\d+_order_(\d+)_accept_ratio_(\d+)\.txt)");
    std::regex beta_dir_accept_pattern(R"(output_beta_mode_\d+_order_(\d+)_direction_accept_ratios_(\d+)\.txt)");
    std::regex beta_tracked_arg_pattern(R"(output_beta_mode_\d+_order_(\d+)_tracked_arg_(\d+)\.txt)");
    std::regex beta_hist_pattern(R"(output_beta_mode_\d+_order_(\d+)_histogram_dir_(\d+)_(\d+)\.txt)"); // NEW

    std::regex mc_pattern(R"(output_MC_order_(\d+)_(\d+)\.txt)");
    std::regex mc_seed_pattern(R"(output_MC_order_(\d+)_seed_(\d+)\.txt)");

    // Set to track processed files
    std::set<std::string> files_to_delete;

    // Iterate through all files in save_path
    for (const auto& entry : std::filesystem::directory_iterator(save_path)) {
        if (entry.is_regular_file() && entry.path().extension() == ".txt") {
            std::string filename = entry.path().filename().string();
            std::string full_path = entry.path().string();

            std::smatch match;
            int order_num = -1;
            int job_index = -1;
            bool is_block_std = false;
            bool is_seed_file = false;
            bool is_sigma_file = false;
            bool is_accept_file = false;
            bool is_dir_accept_file = false;
            bool is_tracked_arg_file = false;
            bool is_histogram_file = false; // NEW
            bool is_sign_file = false;

            // Check each pattern and extract job index
            if (std::regex_match(filename, match, sign_pattern)) {
                job_index = std::stoi(match[1]);
                sign_files.push_back({full_path, job_index});
                is_sign_file = true;
            } else if (std::regex_match(filename, match, sign_block_pattern)) {
                job_index = std::stoi(match[1]);
                sign_block_std_files.push_back({full_path, job_index});
                is_sign_file = true;
            } else if (std::regex_match(filename, match, sign_seed_pattern)) {
                job_index = std::stoi(match[1]);
                sign_seed_files.push_back({full_path, job_index});
                is_sign_file = true;
            } else if (std::regex_match(filename, match, sign_sigma_pattern)) {
                job_index = std::stoi(match[1]);
                sign_sigma_files.push_back({full_path, job_index});
                is_sign_file = true;
            } else if (std::regex_match(filename, match, sign_accept_pattern)) {
                job_index = std::stoi(match[1]);
                sign_accept_ratio_files.push_back({full_path, job_index});
                is_sign_file = true;
            } else if (std::regex_match(filename, match, sign_dir_accept_pattern)) {
                job_index = std::stoi(match[1]);
                sign_direction_accept_ratio_files.push_back({full_path, job_index});
                is_sign_file = true;
            } else if (std::regex_match(filename, match, sign_tracked_arg_pattern)) {
                job_index = std::stoi(match[1]);
                sign_tracked_arg_files.push_back({full_path, job_index});
                is_sign_file = true;
            } else if (std::regex_match(filename, match, sign_hist_pattern)) { // NEW
                int dir_index = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                sign_histogram_files.emplace_back(full_path, job_index, dir_index);
                is_sign_file = true; // Still a sign file, but special kind
            } else if (std::regex_match(filename, match, grid_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
            } else if (std::regex_match(filename, match, grid_block_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_block_std = true;
            } else if (std::regex_match(filename, match, grid_seed_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_seed_file = true;
            } else if (std::regex_match(filename, match, grid_sigma_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_sigma_file = true;
            } else if (std::regex_match(filename, match, grid_accept_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_accept_file = true;
            } else if (std::regex_match(filename, match, grid_dir_accept_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_dir_accept_file = true;
            } else if (std::regex_match(filename, match, grid_tracked_arg_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_tracked_arg_file = true;
            } else if (std::regex_match(filename, match, grid_hist_pattern)) { // NEW
                order_num = std::stoi(match[1]);
                int dir_index = std::stoi(match[2]);
                job_index = std::stoi(match[3]);
                is_histogram_file = true;
            } else if (std::regex_match(filename, match, beta_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
            } else if (std::regex_match(filename, match, beta_block_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_block_std = true;
            } else if (std::regex_match(filename, match, beta_seed_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_seed_file = true;
            } else if (std::regex_match(filename, match, beta_sigma_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_sigma_file = true;
            } else if (std::regex_match(filename, match, beta_accept_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_accept_file = true;
            } else if (std::regex_match(filename, match, beta_dir_accept_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_dir_accept_file = true;
            } else if (std::regex_match(filename, match, beta_tracked_arg_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_tracked_arg_file = true;
            } else if (std::regex_match(filename, match, beta_hist_pattern)) { // NEW
                order_num = std::stoi(match[1]);
                int dir_index = std::stoi(match[2]);
                job_index = std::stoi(match[3]);
                is_histogram_file = true;
            } else if (std::regex_match(filename, match, mc_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
            } else if (std::regex_match(filename, match, mc_seed_pattern)) {
                order_num = std::stoi(match[1]);
                job_index = std::stoi(match[2]);
                is_seed_file = true;
            }

            // Add to appropriate map (only if not a sign file)
            if (!is_sign_file && order_num != -1 && job_index != -1) {
                if (is_histogram_file) {
                    int dir_index = std::stoi(match[2]); // Re-extract for safety
                    histogram_files_by_order[order_num].emplace_back(full_path, job_index, dir_index);
                } else if (is_tracked_arg_file) {
                    tracked_arg_files_by_order[order_num].push_back({full_path, job_index});
                } else if (is_seed_file) {
                    seed_files_by_order[order_num].push_back({full_path, job_index});
                } else if (is_block_std) {
                    block_std_files_by_order[order_num].push_back({full_path, job_index});
                } else if (is_sigma_file) {
                    sigma_files_by_order[order_num].push_back({full_path, job_index});
                } else if (is_accept_file) {
                    accept_ratio_files_by_order[order_num].push_back({full_path, job_index});
                } else if (is_dir_accept_file) {
                    direction_accept_ratio_files_by_order[order_num].push_back({full_path, job_index});
                } else {
                    regular_files_by_order[order_num].push_back({full_path, job_index});
                }
            }
        }
    }

    // Lambda to sort files by job index
    auto sort_by_job_index = [](auto& files) {
        std::sort(files.begin(), files.end(),
                  [](const auto& a, const auto& b) { return a.second < b.second; });
    };

    // Lambda to sort histogram files by job index, then by direction index
    auto sort_hist_files = [](auto& files) {
        std::sort(files.begin(), files.end(),
                  [](const auto& a, const auto& b) {
                      if (std::get<1>(a) != std::get<1>(b)) {
                          return std::get<1>(a) < std::get<1>(b); // sort by job_index
                      }
                      return std::get<2>(a) < std::get<2>(b); // then by dir_index
                  });
    };

    // Combine sign files (without order number) with job index
    if (!sign_files.empty()) {
        sort_by_job_index(sign_files);
        std::string output_filename = std::string(save_path) + "/output_I.txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : sign_files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine sign block std files (without order number) with job index
    if (!sign_block_std_files.empty()) {
        sort_by_job_index(sign_block_std_files);
        std::string output_filename = std::string(save_path) + "/output_I_block_std.txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : sign_block_std_files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine sign seed files (without order number) with job index
    if (!sign_seed_files.empty()) {
        sort_by_job_index(sign_seed_files);
        std::string output_filename = std::string(save_path) + "/output_seeds.txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : sign_seed_files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine sign sigma files (without order number) with job index
    if (!sign_sigma_files.empty()) {
        sort_by_job_index(sign_sigma_files);
        std::string output_filename = std::string(save_path) + "/output_sigmas.txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : sign_sigma_files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine sign accept ratio files (without order number) with job index
    if (!sign_accept_ratio_files.empty()) {
        sort_by_job_index(sign_accept_ratio_files);
        std::string output_filename = std::string(save_path) + "/output_accept_ratios.txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : sign_accept_ratio_files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine sign direction accept ratio files (without order number) with job index
    if (!sign_direction_accept_ratio_files.empty()) {
        sort_by_job_index(sign_direction_accept_ratio_files);
        std::string output_filename = std::string(save_path) + "/output_direction_accept_ratios.txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : sign_direction_accept_ratio_files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine sign tracked arg files (without order number) with job index
    if (!sign_tracked_arg_files.empty()) {
        sort_by_job_index(sign_tracked_arg_files);
        std::string output_filename = std::string(save_path) + "/output_tracked_arg.txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : sign_tracked_arg_files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine sign histogram files - NEW
    if (!sign_histogram_files.empty()) {
        sort_hist_files(sign_histogram_files);
        std::string output_filename = std::string(save_path) + "/output_histograms.txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            output_file << "job_idx,direction_idx,bin_idx,count,probability" << std::endl;
            for (const auto& [input_file, job_idx, dir_idx] : sign_histogram_files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    std::getline(input, line); // Skip header 1
                    std::getline(input, line); // Skip header 2
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << dir_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine regular files by order with job index
    for (auto& [order_num, files] : regular_files_by_order) {
        sort_by_job_index(files);
        std::string output_filename = std::string(save_path) + "/output_I_order_" + std::to_string(order_num) + ".txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine block std files by order with job index
    for (auto& [order_num, files] : block_std_files_by_order) {
        sort_by_job_index(files);
        std::string output_filename = std::string(save_path) + "/output_block_std_order_" + std::to_string(order_num) + ".txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine seed files by order with job index
    for (auto& [order_num, files] : seed_files_by_order) {
        sort_by_job_index(files);
        std::string output_filename = std::string(save_path) + "/output_seeds_order_" + std::to_string(order_num) + ".txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine sigma files by order with job index
    for (auto& [order_num, files] : sigma_files_by_order) {
        sort_by_job_index(files);
        std::string output_filename = std::string(save_path) + "/output_sigmas_order_" + std::to_string(order_num) + ".txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine accept ratio files by order with job index
    for (auto& [order_num, files] : accept_ratio_files_by_order) {
        sort_by_job_index(files);
        std::string output_filename = std::string(save_path) + "/output_accept_ratios_order_" + std::to_string(order_num) + ".txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine direction accept ratio files by order with job index
    for (auto& [order_num, files] : direction_accept_ratio_files_by_order) {
        sort_by_job_index(files);
        std::string output_filename = std::string(save_path) + "/output_direction_accept_ratios_order_" + std::to_string(order_num) + ".txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine tracked arg files by order with job index
    for (auto& [order_num, files] : tracked_arg_files_by_order) {
        sort_by_job_index(files);
        std::string output_filename = std::string(save_path) + "/output_tracked_arg_order_" + std::to_string(order_num) + ".txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            for (const auto& [input_file, job_idx] : files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }

    // Combine histogram files by order - NEW
    for (auto& [order_num, files] : histogram_files_by_order) {
        sort_hist_files(files);
        std::string output_filename = std::string(save_path) + "/output_histograms_order_" + std::to_string(order_num) + ".txt";
        std::ofstream output_file(output_filename);

        if (output_file.is_open()) {
            output_file << "job_idx,direction_idx,bin_idx,count,probability" << std::endl;
            for (const auto& [input_file, job_idx, dir_idx] : files) {
                std::ifstream input(input_file);
                if (input.is_open()) {
                    std::string line;
                    std::getline(input, line); // Skip header 1
                    std::getline(input, line); // Skip header 2
                    while (std::getline(input, line)) {
                        output_file << job_idx << "," << dir_idx << "," << line << std::endl;
                    }
                    input.close();
                    files_to_delete.insert(input_file);
                }
            }
            output_file.close();
        }
    }


    // Delete the original files
    for (const auto& file : files_to_delete) {
        try {
            std::filesystem::remove(file);
        } catch (const std::exception& e) {
            std::cerr << "Error deleting " << file << ": " << e.what() << std::endl;
        }
    }

    std::cout << "File combination complete." << std::endl;

   return 0;
}