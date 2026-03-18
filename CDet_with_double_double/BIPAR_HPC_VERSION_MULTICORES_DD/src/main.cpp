#include "../utility/resize_and_initialize.h"
#include "MCMC.h"
#include "integrand_sublattice_randomly_sampled.h"
#include "integrand_sublattice_summed.h"
#include "load_data.h"
#include "../utility/types.h"
#include "../utility/check_range.h"
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
// ORIGINAL UNCHANGED FUNCTIONS
//================================================================================

// Create lookup table at compile time using constexpr
std::array<unsigned char, 256> createMSBTable() {
    std::array<unsigned char, 256> table{};
    for (int i = 0; i < 256; ++i) {
        unsigned char msb = 0;
        if (i != 0) {
            unsigned char val = i;
            msb = 128;
            while ((val & msb) == 0) {
                msb >>= 1;
            }
        }
        table[i] = msb;
    }
    return table;
}

// Convert string to MP_Real
MP_Real string_to_mp_real(const std::string& s) {
    return MP_Real(s);
}

// Convert string to Real (double)
Real string_to_real(const std::string& s) {
    return std::stod(s);
}

// DoubleDouble version of integrand that does all computation in DoubleDouble
DoubleDouble integrand_dd(const DoubleDouble* args, const int order_num, const int i, const int j) {
    if (lnZ_mode) {
        return integrand_temp_dd(
            DoubleDouble(beta_vec[j]), DoubleDouble(U), args,
            g_vec_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
            g_vec_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            DoubleDouble(eps), cheby_degree,
            DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
            DoubleDouble(Delta_up), DoubleDouble(Delta_down),
            DoubleDouble(tol_fpm)
        );
    }
    if (density_mode) {
        // For density mode, we need to calculate the difference between mu_h and mu_l
        DoubleDouble result_mu_h = integrand_temp_dd(
            DoubleDouble(beta_vec[j]), DoubleDouble(U), args,
            g_vec_mu_h_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
            g_vec_mu_h_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            DoubleDouble(eps), cheby_degree,
            DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
            DoubleDouble(Delta_up), DoubleDouble(Delta_down),
            DoubleDouble(tol_fpm)
        );

        DoubleDouble result_mu_l = integrand_temp_dd(
            DoubleDouble(beta_vec[j]), DoubleDouble(U), args,
            g_vec_mu_l_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
            g_vec_mu_l_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            DoubleDouble(eps), cheby_degree,
            DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
            DoubleDouble(Delta_up), DoubleDouble(Delta_down),
            DoubleDouble(tol_fpm)
        );

        return (result_mu_h - result_mu_l) / (DoubleDouble(beta_vec[j]) * (DoubleDouble(mu_h) - DoubleDouble(mu_l)));
    }
    if (compressibility_mode) {
        DoubleDouble val_0 = integrand_temp_dd(
            DoubleDouble(beta_vec[j]), DoubleDouble(U), args,
            g_vec_mu_h_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
            g_vec_mu_h_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            DoubleDouble(eps), cheby_degree,
            DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
            DoubleDouble(Delta_up), DoubleDouble(Delta_down),
            DoubleDouble(tol_fpm)
        );
        DoubleDouble val_1 = integrand_temp_dd(
            DoubleDouble(beta_vec[j]), DoubleDouble(U), args,
            g_vec_mu_l_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
            g_vec_mu_l_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            DoubleDouble(eps), cheby_degree,
            DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
            DoubleDouble(Delta_up), DoubleDouble(Delta_down),
            DoubleDouble(tol_fpm)
        );
        DoubleDouble val_2 = integrand_temp_dd(
            DoubleDouble(beta_vec[j]), DoubleDouble(U), args,
            g_vec_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
            g_vec_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            DoubleDouble(eps), cheby_degree,
            DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
            DoubleDouble(Delta_up), DoubleDouble(Delta_down),
            DoubleDouble(tol_fpm)
        );
        return (val_0 + val_1 - DoubleDouble(2.0) * val_2) / (DoubleDouble(beta_vec[j]) * dd_pow((DoubleDouble(mu_h) - DoubleDouble(mu_l))/DoubleDouble(2.0), DoubleDouble(2.0)));
    }
    if (energy_mode) {
        DoubleDouble val_beta_h = integrand_temp_dd(
            DoubleDouble(beta_h_vec[j]), DoubleDouble(U), args,
            g_vec_beta_h_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
            g_vec_beta_h_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            DoubleDouble(eps), cheby_degree,
            DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
            DoubleDouble(Delta_up), DoubleDouble(Delta_down),
            DoubleDouble(tol_fpm)
        );
        DoubleDouble val_beta_l = integrand_temp_dd(
            DoubleDouble(beta_l_vec[j]), DoubleDouble(U), args,
            g_vec_beta_l_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
            g_vec_beta_l_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            DoubleDouble(eps), cheby_degree,
            DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
            DoubleDouble(Delta_up), DoubleDouble(Delta_down),
            DoubleDouble(tol_fpm)
        );
        return -(val_beta_h - val_beta_l) / (DoubleDouble(beta_h_vec[j]) - DoubleDouble(beta_l_vec[j]));
    }
    if (double_occupancy_mode) {
         DoubleDouble val_0 = integrand_temp_dd(
            DoubleDouble(beta_vec[j]), DoubleDouble(U_h), args,
            g_vec_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
            g_vec_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            DoubleDouble(eps), cheby_degree,
            DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
            DoubleDouble(Delta_up), DoubleDouble(Delta_down),
            DoubleDouble(tol_fpm));
        DoubleDouble val_1 = integrand_temp_dd(
           DoubleDouble(beta_vec[j]), DoubleDouble(U_l), args,
           g_vec_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
           g_vec_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
           order_num, N1_vec[i], N2_vec[i],
           DoubleDouble(eps), cheby_degree,
           DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
           DoubleDouble(Delta_up), DoubleDouble(Delta_down),
           DoubleDouble(tol_fpm));
        return -(val_0 - val_1) / (DoubleDouble(beta_vec[j]) * (DoubleDouble(U_h) - DoubleDouble(U_l)));
    }
    if (staggered_magnetization_mode) {
        // (n_A,up - n_A,down - (n_B,up - n_B,down)) / 2
        DoubleDouble val_0 = integrand_temp_dd(
            DoubleDouble(beta_vec[j]), DoubleDouble(U), args,
            g_vec_H_h_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
            g_vec_H_h_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            DoubleDouble(eps), cheby_degree,
            DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
            DoubleDouble(Delta_up), DoubleDouble(Delta_down),
            DoubleDouble(tol_fpm)
        );
        DoubleDouble val_1 = integrand_temp_dd(
            DoubleDouble(beta_vec[j]), DoubleDouble(U), args,
            g_vec_H_l_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
            g_vec_H_l_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            DoubleDouble(eps), cheby_degree,
            DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
            DoubleDouble(Delta_up), DoubleDouble(Delta_down),
            DoubleDouble(tol_fpm)
        );
        return -DoubleDouble(0.5) * (val_0 - val_1) / (DoubleDouble(beta_vec[j]) * (DoubleDouble(H_h) - DoubleDouble(H_l)));
    } else {
        throw std::runtime_error("Invalid mode");
    }
}

// Sublattice-summed version of integrand_dd
DoubleDouble integrand_dd_sublattice_summed(const DoubleDouble* args, const int order_num, const int i, const int j) {
    if (lnZ_mode) {
        return integrand_sublattice_summed_dd(
            DoubleDouble(beta_vec[j]), DoubleDouble(U), args,
            g_vec_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
            g_vec_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            DoubleDouble(eps), cheby_degree,
            DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
            DoubleDouble(Delta_up), DoubleDouble(Delta_down),
            DoubleDouble(tol_fpm)
        );
    }
    if (density_mode) {
        DoubleDouble result_mu_h = integrand_sublattice_summed_dd(
            DoubleDouble(beta_vec[j]), DoubleDouble(U), args,
            g_vec_mu_h_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
            g_vec_mu_h_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            DoubleDouble(eps), cheby_degree,
            DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
            DoubleDouble(Delta_up), DoubleDouble(Delta_down),
            DoubleDouble(tol_fpm)
        );

        DoubleDouble result_mu_l = integrand_sublattice_summed_dd(
            DoubleDouble(beta_vec[j]), DoubleDouble(U), args,
            g_vec_mu_l_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
            g_vec_mu_l_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            DoubleDouble(eps), cheby_degree,
            DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
            DoubleDouble(Delta_up), DoubleDouble(Delta_down),
            DoubleDouble(tol_fpm)
        );

        return (result_mu_h - result_mu_l) / (DoubleDouble(beta_vec[j]) * (DoubleDouble(mu_h) - DoubleDouble(mu_l)));
    }
    if (compressibility_mode) {
        DoubleDouble val_0 = integrand_sublattice_summed_dd(
            DoubleDouble(beta_vec[j]), DoubleDouble(U), args,
            g_vec_mu_h_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
            g_vec_mu_h_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            DoubleDouble(eps), cheby_degree,
            DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
            DoubleDouble(Delta_up), DoubleDouble(Delta_down),
            DoubleDouble(tol_fpm)
        );
        DoubleDouble val_1 = integrand_sublattice_summed_dd(
            DoubleDouble(beta_vec[j]), DoubleDouble(U), args,
            g_vec_mu_l_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
            g_vec_mu_l_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            DoubleDouble(eps), cheby_degree,
            DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
            DoubleDouble(Delta_up), DoubleDouble(Delta_down),
            DoubleDouble(tol_fpm)
        );
        DoubleDouble val_2 = integrand_sublattice_summed_dd(
            DoubleDouble(beta_vec[j]), DoubleDouble(U), args,
            g_vec_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
            g_vec_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            DoubleDouble(eps), cheby_degree,
            DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
            DoubleDouble(Delta_up), DoubleDouble(Delta_down),
            DoubleDouble(tol_fpm)
        );
        return (val_0 + val_1 - DoubleDouble(2.0) * val_2) / (DoubleDouble(beta_vec[j]) * dd_pow((DoubleDouble(mu_h) - DoubleDouble(mu_l))/DoubleDouble(2.0), DoubleDouble(2.0)));
    }
    if (energy_mode) {
        DoubleDouble val_beta_h = integrand_sublattice_summed_dd(
            DoubleDouble(beta_h_vec[j]), DoubleDouble(U), args,
            g_vec_beta_h_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
            g_vec_beta_h_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            DoubleDouble(eps), cheby_degree,
            DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
            DoubleDouble(Delta_up), DoubleDouble(Delta_down),
            DoubleDouble(tol_fpm)
        );
        DoubleDouble val_beta_l = integrand_sublattice_summed_dd(
            DoubleDouble(beta_l_vec[j]), DoubleDouble(U), args,
            g_vec_beta_l_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
            g_vec_beta_l_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            DoubleDouble(eps), cheby_degree,
            DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
            DoubleDouble(Delta_up), DoubleDouble(Delta_down),
            DoubleDouble(tol_fpm)
        );
        return -(val_beta_h - val_beta_l) / (DoubleDouble(beta_h_vec[j]) - DoubleDouble(beta_l_vec[j]));
    }
    if (double_occupancy_mode) {
        DoubleDouble val_0 = integrand_sublattice_summed_dd(
            DoubleDouble(beta_vec[j]), DoubleDouble(U_h), args,
            g_vec_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
            g_vec_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            DoubleDouble(eps), cheby_degree,
            DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
            DoubleDouble(Delta_up), DoubleDouble(Delta_down),
            DoubleDouble(tol_fpm));
        DoubleDouble val_1 = integrand_sublattice_summed_dd(
           DoubleDouble(beta_vec[j]), DoubleDouble(U_l), args,
           g_vec_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
           g_vec_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
           order_num, N1_vec[i], N2_vec[i],
           DoubleDouble(eps), cheby_degree,
           DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
           DoubleDouble(Delta_up), DoubleDouble(Delta_down),
           DoubleDouble(tol_fpm));
        return -(val_0 - val_1) / (DoubleDouble(beta_vec[j]) * (DoubleDouble(U_h) - DoubleDouble(U_l)));
    }
    if (staggered_magnetization_mode) {
        DoubleDouble val_0 = integrand_sublattice_summed_dd(
            DoubleDouble(beta_vec[j]), DoubleDouble(U), args,
            g_vec_H_h_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
            g_vec_H_h_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            DoubleDouble(eps), cheby_degree,
            DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
            DoubleDouble(Delta_up), DoubleDouble(Delta_down),
            DoubleDouble(tol_fpm)
        );
        DoubleDouble val_1 = integrand_sublattice_summed_dd(
            DoubleDouble(beta_vec[j]), DoubleDouble(U), args,
            g_vec_H_l_flattened["up_" + std::to_string(i) + "_" + std::to_string(j)],
            g_vec_H_l_flattened["down_" + std::to_string(i) + "_" + std::to_string(j)],
            order_num, N1_vec[i], N2_vec[i],
            DoubleDouble(eps), cheby_degree,
            DoubleDouble(alpha_shift_up), DoubleDouble(alpha_shift_down),
            DoubleDouble(Delta_up), DoubleDouble(Delta_down),
            DoubleDouble(tol_fpm)
        );
        return -DoubleDouble(0.5) * (val_0 - val_1) / (DoubleDouble(beta_vec[j]) * (DoubleDouble(H_h) - DoubleDouble(H_l)));
    } else {
        throw std::runtime_error("Invalid mode");
    }
}

// New helper function to select order based on job_index for each case
int selectOrderFromJobIndex(int set_num, unsigned long int job_index) {
    switch (set_num) {
        case 1: // Intervals of 40
            if (job_index <= 4) return 1;      // order 1
            else if (job_index <= 20) return 2; // order 2
            else if (job_index <= 36) return 3; // order 3
            else return 4;                       // order 4

        case 2: // Intervals of 50
            if (job_index <= 20) return 5;      // order 5
            else if (job_index <= 40) return 6; // order 6
            else return 7;                       // order 7

        case 3: // Intervals of 80
            return 8;      // order 8

        default:
            throw std::runtime_error("Invalid set_num: " + std::to_string(set_num));
    }
}

//================================================================================
// MODIFIED MCMC RUN FUNCTIONS WITHOUT CHECKPOINT SAVING
//================================================================================

void run_all_dd(unsigned long int job_index) {

    int actual_iter_bits_thres = 26;

    int order_num;
    if (MCMC_mode == 0) {
        order_num = set_num; // For MCMC_mode = 0, use set_num as order directly
    } else {
        order_num = selectOrderFromJobIndex(set_num, job_index);
    }

    // Get parameters for this specific order
    auto sigma = DoubleDouble(sigma_vec[order_num]);
    int iter_bits = iter_bits_vec[order_num];
    int warm_bits = warm_bits_vec[order_num];
    int MC_bits = MC_bits_vec[order_num];
    const int DIM = 4 * order_num;

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
    // ---------------------- RNG setup (drop-in replacement) ----------------------
    unsigned long long seed_start = seed_offset;
    for (int k = 1; k < order_num; ++k) {
        seed_start += static_cast<unsigned long long>(seed_interval) * 20000ULL;
    }
    seed_start += static_cast<unsigned long long>(job_index - 1) * seed_interval;
    unsigned long long seed_end = seed_start + seed_interval;

    // Use a temporary generator from std::random_device
    std::random_device rd;
    std::mt19937 gen_seed(rd());
    std::uniform_int_distribution<unsigned long long> dist_seed(seed_start, seed_end - 1);

    // Draw a 64-bit seed candidate
    unsigned long long seed = dist_seed(gen_seed);

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

        auto sign_func = [](const DoubleDouble* args, int order_num_param) {
            return integrand_dd(args, order_num_param, 0, 0);
        };

        auto [grid_sign_dd, block_variances_dd, final_sigma_dd, accept_ratio_dd] = run_sign_dd(sign_func, sigma, order_num, DIM,
                                                             gen, iter_bits, warm_bits, seed);
        Real grid_sign = grid_sign_dd.to_double();
        Real final_sigma = final_sigma_dd.to_double();
        Real accept_ratio = accept_ratio_dd.to_double();

        std::string filename = std::string(save_path) + "/output_sign_" + std::to_string(job_index) + ".txt";
        std::ofstream file(filename, std::ios::app);
        if (file.is_open()) {
            file << std::fixed << std::setprecision(15) << grid_sign << std::endl;
        }

        std::string block_std_filename = std::string(save_path) + "/output_sign_block_std_" + std::to_string(job_index) + ".txt";
        std::ofstream block_std_file(block_std_filename);
        if (block_std_file.is_open()) {
            block_std_file << std::fixed << std::setprecision(15);
            for (size_t i = 0; i < block_variances_dd.size(); ++i) {
                block_std_file << DD_SQRT(block_variances_dd[i]).to_double() << (i < block_variances_dd.size() - 1 ? "," : "");
            }
            block_std_file << std::endl;
        }

        // Save final sigma
        std::string sigma_filename = std::string(save_path) + "/output_sign_sigma_" + std::to_string(job_index) + ".txt";
        std::ofstream sigma_file(sigma_filename);
        if (sigma_file.is_open()) {
            sigma_file << std::fixed << std::setprecision(15) << final_sigma << std::endl;
        }

        // Save accept ratio
        std::string accept_filename = std::string(save_path) + "/output_sign_accept_ratio_" + std::to_string(job_index) + ".txt";
        std::ofstream accept_file(accept_filename);
        if (accept_file.is_open()) {
            accept_file << std::fixed << std::setprecision(15) << accept_ratio << std::endl;
        }
    }
    // Grid reduction modes: one mode per chain transition
    else if (MCMC_mode >= 1 && MCMC_mode <= num_grid_chains) {
        int chain_idx = MCMC_mode - 1; // 0-based index for the chain

        // Save seed
        std::string seed_filename = std::string(save_path) + "/output_grid_mode_" + std::to_string(MCMC_mode) +
                                   "_order_" + std::to_string(order_num) + "_seed_" + std::to_string(job_index) + ".txt";
        std::ofstream seed_file(seed_filename);
        if (seed_file.is_open()) {
            seed_file << seed << std::endl;
            seed_file.close();
        }

        // Function containers for a single chain
        std::vector<std::function<DoubleDouble(const DoubleDouble*, int)>> measure_funcs(1);
        std::vector<std::function<DoubleDouble(const DoubleDouble*, int)>> ref_funcs(1);

        // Setup measure function (current step in the chain)
        measure_funcs[0] = [grid_idx = chain_idx](const DoubleDouble* args, int order_num_param) {
            return DD_ABS(integrand_dd(args, order_num_param, grid_idx, 0));
        };

        // Setup reference function (next step in the chain)
        ref_funcs[0] = [grid_idx = chain_idx + 1](const DoubleDouble* args, int order_num_param) {
            return DD_ABS(integrand_dd(args, order_num_param, grid_idx, 0));
        };

        std::vector<bool> if_abs_ref(1, true);
        std::vector<bool> if_abs_measure(1, true);

        int actual_iter_bits = iter_bits - 3;
        if (actual_iter_bits < actual_iter_bits_thres) actual_iter_bits = actual_iter_bits_thres;

        auto [result_dd, ref_block_variances_dd, measure_block_variances_dd, final_sigma_dd, accept_ratios_dd] =
            run_dd(measure_funcs, ref_funcs, if_abs_ref, if_abs_measure,
                   sigma, order_num, DIM, gen, actual_iter_bits, warm_bits, 0, seed);

        Real result = result_dd.to_double();
        Real final_sigma = final_sigma_dd.to_double();
        Real accept_ratio = accept_ratios_dd[0].to_double(); // Single chain

        std::string filename = std::string(save_path) + "/output_grid_mode_" + std::to_string(MCMC_mode) +
                               "_order_" + std::to_string(order_num) + "_" + std::to_string(job_index) + ".txt";
        std::ofstream file(filename, std::ios::app);
        if (file.is_open()) {
            file << std::fixed << std::setprecision(15) << result << std::endl;
        }

        std::string block_std_filename = std::string(save_path) + "/output_grid_mode_" + std::to_string(MCMC_mode) +
                                         "_order_" + std::to_string(order_num) + "_block_std_" + std::to_string(job_index) + ".txt";
        std::ofstream block_std_file(block_std_filename);
        if (block_std_file.is_open()) {
            block_std_file << std::fixed << std::setprecision(15);
            for (size_t i = 0; i < ref_block_variances_dd[0].size(); ++i) {
                block_std_file << "0,O1," << DD_SQRT(ref_block_variances_dd[0][i]).to_double() << std::endl;
            }
        }

        // Save final sigma
        std::string sigma_filename = std::string(save_path) + "/output_grid_mode_" + std::to_string(MCMC_mode) +
                                    "_order_" + std::to_string(order_num) + "_sigma_" + std::to_string(job_index) + ".txt";
        std::ofstream sigma_file(sigma_filename);
        if (sigma_file.is_open()) {
            sigma_file << std::fixed << std::setprecision(15) << final_sigma << std::endl;
        }

        // Save accept ratio
        std::string accept_filename = std::string(save_path) + "/output_grid_mode_" + std::to_string(MCMC_mode) +
                                     "_order_" + std::to_string(order_num) + "_accept_ratio_" + std::to_string(job_index) + ".txt";
        std::ofstream accept_file(accept_filename);
        if (accept_file.is_open()) {
            accept_file << std::fixed << std::setprecision(15) << accept_ratio << std::endl;
        }
    }
    // Beta reduction modes: one mode per chain transition
    else if (MCMC_mode > num_grid_chains && MCMC_mode <= num_grid_chains + num_beta_chains) {
        int chain_idx = MCMC_mode - num_grid_chains - 1; // 0-based index for the chain
        int grid_idx = num_grid_chains; // Beta chains run at the final grid size

        // Save seed
        std::string seed_filename = std::string(save_path) + "/output_beta_mode_" + std::to_string(MCMC_mode) +
                                   "_order_" + std::to_string(order_num) + "_seed_" + std::to_string(job_index) + ".txt";
        std::ofstream seed_file(seed_filename);
        if (seed_file.is_open()) {
            seed_file << seed << std::endl;
            seed_file.close();
        }

        std::vector<std::function<DoubleDouble(const DoubleDouble*, int)>> measure_funcs(1);
        std::vector<std::function<DoubleDouble(const DoubleDouble*, int)>> ref_funcs(1);

        measure_funcs[0] = [=](const DoubleDouble* args, int order_num_param) {
            return DD_ABS(integrand_dd(args, order_num_param, grid_idx, chain_idx));
        };

        ref_funcs[0] = [=](const DoubleDouble* args, int order_num_param) {
            return DD_ABS(integrand_dd(args, order_num_param, grid_idx, chain_idx + 1));
        };

        std::vector<bool> if_abs_ref(1, true);
        std::vector<bool> if_abs_measure(1, true);

        int actual_iter_bits = iter_bits - 3;
        if (actual_iter_bits < actual_iter_bits_thres) actual_iter_bits = actual_iter_bits_thres;

        auto [result_dd, ref_block_variances_dd, measure_block_variances_dd, final_sigma_dd, accept_ratios_dd] =
            run_dd(measure_funcs, ref_funcs, if_abs_ref, if_abs_measure,
                   sigma, order_num, DIM, gen, actual_iter_bits, warm_bits, 0, seed);

        Real result = result_dd.to_double();
        Real final_sigma = final_sigma_dd.to_double();
        Real accept_ratio = accept_ratios_dd[0].to_double();

        std::string filename = std::string(save_path) + "/output_beta_mode_" + std::to_string(MCMC_mode) +
                               "_order_" + std::to_string(order_num) + "_" + std::to_string(job_index) + ".txt";
        std::ofstream file(filename, std::ios::app);
        if (file.is_open()) {
            file << std::fixed << std::setprecision(15) << result << std::endl;
        }

        std::string block_std_filename = std::string(save_path) + "/output_beta_mode_" + std::to_string(MCMC_mode) +
                                         "_order_" + std::to_string(order_num) + "_block_std_" + std::to_string(job_index) + ".txt";
        std::ofstream block_std_file(block_std_filename);
        if (block_std_file.is_open()) {
            block_std_file << std::fixed << std::setprecision(15);
            for (size_t i = 0; i < ref_block_variances_dd[0].size(); ++i) {
                block_std_file << "0,O1," << DD_SQRT(ref_block_variances_dd[0][i]).to_double() << std::endl;
            }
        }

        // Save final sigma
        std::string sigma_filename = std::string(save_path) + "/output_beta_mode_" + std::to_string(MCMC_mode) +
                                    "_order_" + std::to_string(order_num) + "_sigma_" + std::to_string(job_index) + ".txt";
        std::ofstream sigma_file(sigma_filename);
        if (sigma_file.is_open()) {
            sigma_file << std::fixed << std::setprecision(15) << final_sigma << std::endl;
        }

        // Save accept ratio
        std::string accept_filename = std::string(save_path) + "/output_beta_mode_" + std::to_string(MCMC_mode) +
                                     "_order_" + std::to_string(order_num) + "_accept_ratio_" + std::to_string(job_index) + ".txt";
        std::ofstream accept_file(accept_filename);
        if (accept_file.is_open()) {
            accept_file << std::fixed << std::setprecision(15) << accept_ratio << std::endl;
        }
    }
    // Direct Monte Carlo integration for the final configuration
    else if (MCMC_mode == num_grid_chains + num_beta_chains + 1) {
        // Save seed
        std::string seed_filename = std::string(save_path) + "/output_MC_order_" + std::to_string(order_num) +
                                   "_seed_" + std::to_string(job_index) + ".txt";
        std::ofstream seed_file(seed_filename);
        if (seed_file.is_open()) {
            seed_file << seed << std::endl;
            seed_file.close();
        }

        const int dim = 4 * order_num;
        unsigned long numSamples = 1 << MC_bits;
        DoubleDouble sum(0);

        DD1DVector x(dim);
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        int largest_grid_idx = num_grid_chains;
        int largest_beta_idx = num_beta_chains;

        for (unsigned long i = 0; i < numSamples; i++) {
            for (int j = 0; j < dim; j++) {
                x[j] = DoubleDouble(dist(gen));
            }
            sum += DD_ABS(integrand_dd(x.data(), order_num, largest_grid_idx, largest_beta_idx));
        }

        DoubleDouble integral_dd = sum / DoubleDouble(numSamples);
        Real integral = integral_dd.to_double();

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

void run_all_dd_sublattices_summed(unsigned long int job_index) {

    int actual_iter_bits_thres = 20;

    int order_num;
    if (MCMC_mode == 0) {
        order_num = set_num; // For MCMC_mode = 0, use set_num as order directly
    } else {
        order_num = selectOrderFromJobIndex(set_num, job_index);
    }

    // Get parameters for this specific order
    auto sigma = DoubleDouble(sigma_vec[order_num]);
    int iter_bits = iter_bits_vec[order_num];
    int warm_bits = warm_bits_vec[order_num];
    int MC_bits = MC_bits_vec[order_num];
    const int DIM = 3 * order_num; // Use 3 * order_num for sublattice summed

    // Determine the number of chain transitions required
    const int num_grid_chains = N1_vec.size() > 1 ? N1_vec.size() - 1 : 0;
    const int num_beta_chains = beta_vec.size() > 1 ? beta_vec.size() - 1 : 0;

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

        auto sign_func = [](const DoubleDouble* args, int order_num_param) {
            return integrand_dd_sublattice_summed(args, order_num_param, 0, 0);
        };

        auto [grid_sign_dd, block_variances_dd, final_sigma_dd, accept_ratio_dd] = run_sign_dd(sign_func, sigma, order_num, DIM,
                                                             gen, iter_bits, warm_bits, seed);
        Real grid_sign = grid_sign_dd.to_double();
        Real final_sigma = final_sigma_dd.to_double();
        Real accept_ratio = accept_ratio_dd.to_double();

        std::string filename = std::string(save_path) + "/output_sign_" + std::to_string(job_index) + ".txt";
        std::ofstream file(filename, std::ios::app);
        if (file.is_open()) {
            file << std::fixed << std::setprecision(15) << grid_sign << std::endl;
        }

        std::string block_std_filename = std::string(save_path) + "/output_sign_block_std_" + std::to_string(job_index) + ".txt";
        std::ofstream block_std_file(block_std_filename);
        if (block_std_file.is_open()) {
            block_std_file << std::fixed << std::setprecision(15);
            for (size_t i = 0; i < block_variances_dd.size(); ++i) {
                 block_std_file << DD_SQRT(block_variances_dd[i]).to_double() << (i < block_variances_dd.size() - 1 ? "," : "");
            }
            block_std_file << std::endl;
        }

        // Save final sigma
        std::string sigma_filename = std::string(save_path) + "/output_sign_sigma_" + std::to_string(job_index) + ".txt";
        std::ofstream sigma_file(sigma_filename);
        if (sigma_file.is_open()) {
            sigma_file << std::fixed << std::setprecision(15) << final_sigma << std::endl;
        }

        // Save accept ratio
        std::string accept_filename = std::string(save_path) + "/output_sign_accept_ratio_" + std::to_string(job_index) + ".txt";
        std::ofstream accept_file(accept_filename);
        if (accept_file.is_open()) {
            accept_file << std::fixed << std::setprecision(15) << accept_ratio << std::endl;
        }
    }
    // Grid reduction modes: one mode per chain transition
    else if (MCMC_mode >= 1 && MCMC_mode <= num_grid_chains) {
        int chain_idx = MCMC_mode - 1;

        // Save seed
        std::string seed_filename = std::string(save_path) + "/output_grid_mode_" + std::to_string(MCMC_mode) +
                                   "_order_" + std::to_string(order_num) + "_seed_" + std::to_string(job_index) + ".txt";
        std::ofstream seed_file(seed_filename);
        if (seed_file.is_open()) {
            seed_file << seed << std::endl;
            seed_file.close();
        }

        std::vector<std::function<DoubleDouble(const DoubleDouble*, int)>> measure_funcs(1);
        std::vector<std::function<DoubleDouble(const DoubleDouble*, int)>> ref_funcs(1);

        measure_funcs[0] = [grid_idx = chain_idx](const DoubleDouble* args, int order_num_param) {
            return DD_ABS(integrand_dd_sublattice_summed(args, order_num_param, grid_idx, 0));
        };

        ref_funcs[0] = [grid_idx = chain_idx + 1](const DoubleDouble* args, int order_num_param) {
            return DD_ABS(integrand_dd_sublattice_summed(args, order_num_param, grid_idx, 0));
        };

        std::vector<bool> if_abs_ref(1, true);
        std::vector<bool> if_abs_measure(1, true);

        int actual_iter_bits = iter_bits - 3;
        if (actual_iter_bits < actual_iter_bits_thres) actual_iter_bits = actual_iter_bits_thres;

        auto [result_dd, ref_block_variances_dd, measure_block_variances_dd, final_sigma_dd, accept_ratios_dd] =
            run_dd(measure_funcs, ref_funcs, if_abs_ref, if_abs_measure,
                  sigma, order_num, DIM, gen, actual_iter_bits, warm_bits, 0, seed);

        Real result = result_dd.to_double();
        Real final_sigma = final_sigma_dd.to_double();
        Real accept_ratio = accept_ratios_dd[0].to_double();

        std::string filename = std::string(save_path) + "/output_grid_mode_" + std::to_string(MCMC_mode) +
                               "_order_" + std::to_string(order_num) + "_" + std::to_string(job_index) + ".txt";
        std::ofstream file(filename, std::ios::app);
        if (file.is_open()) {
            file << std::fixed << std::setprecision(15) << result << std::endl;
        }

        std::string block_std_filename = std::string(save_path) + "/output_grid_mode_" + std::to_string(MCMC_mode) +
                                         "_order_" + std::to_string(order_num) + "_block_std_" + std::to_string(job_index) + ".txt";
        std::ofstream block_std_file(block_std_filename);
        if (block_std_file.is_open()) {
            block_std_file << std::fixed << std::setprecision(15);
            for (size_t i = 0; i < ref_block_variances_dd[0].size(); ++i) {
                block_std_file << "0,O1," << DD_SQRT(ref_block_variances_dd[0][i]).to_double() << std::endl;
            }
        }

        // Save final sigma
        std::string sigma_filename = std::string(save_path) + "/output_grid_mode_" + std::to_string(MCMC_mode) +
                                    "_order_" + std::to_string(order_num) + "_sigma_" + std::to_string(job_index) + ".txt";
        std::ofstream sigma_file(sigma_filename);
        if (sigma_file.is_open()) {
            sigma_file << std::fixed << std::setprecision(15) << final_sigma << std::endl;
        }

        // Save accept ratio
        std::string accept_filename = std::string(save_path) + "/output_grid_mode_" + std::to_string(MCMC_mode) +
                                     "_order_" + std::to_string(order_num) + "_accept_ratio_" + std::to_string(job_index) + ".txt";
        std::ofstream accept_file(accept_filename);
        if (accept_file.is_open()) {
            accept_file << std::fixed << std::setprecision(15) << accept_ratio << std::endl;
        }
    }
    // Beta reduction modes: one mode per chain transition
    else if (MCMC_mode > num_grid_chains && MCMC_mode <= num_grid_chains + num_beta_chains) {
        int chain_idx = MCMC_mode - num_grid_chains - 1;
        int grid_idx = num_grid_chains;

        // Save seed
        std::string seed_filename = std::string(save_path) + "/output_beta_mode_" + std::to_string(MCMC_mode) +
                                   "_order_" + std::to_string(order_num) + "_seed_" + std::to_string(job_index) + ".txt";
        std::ofstream seed_file(seed_filename);
        if (seed_file.is_open()) {
            seed_file << seed << std::endl;
            seed_file.close();
        }

        std::vector<std::function<DoubleDouble(const DoubleDouble*, int)>> measure_funcs(1);
        std::vector<std::function<DoubleDouble(const DoubleDouble*, int)>> ref_funcs(1);

        measure_funcs[0] = [=](const DoubleDouble* args, int order_num_param) {
            return DD_ABS(integrand_dd_sublattice_summed(args, order_num_param, grid_idx, chain_idx));
        };

        ref_funcs[0] = [=](const DoubleDouble* args, int order_num_param) {
            return DD_ABS(integrand_dd_sublattice_summed(args, order_num_param, grid_idx, chain_idx + 1));
        };

        std::vector<bool> if_abs_ref(1, true);
        std::vector<bool> if_abs_measure(1, true);

        int actual_iter_bits = iter_bits - 3;
        if (actual_iter_bits < actual_iter_bits_thres) actual_iter_bits = actual_iter_bits_thres;

        auto [result_dd, ref_block_variances_dd, measure_block_variances_dd, final_sigma_dd, accept_ratios_dd] =
            run_dd(measure_funcs, ref_funcs, if_abs_ref, if_abs_measure,
                   sigma, order_num, DIM, gen, actual_iter_bits, warm_bits, 0, seed);

        Real result = result_dd.to_double();
        Real final_sigma = final_sigma_dd.to_double();
        Real accept_ratio = accept_ratios_dd[0].to_double();

        std::string filename = std::string(save_path) + "/output_beta_mode_" + std::to_string(MCMC_mode) +
                               "_order_" + std::to_string(order_num) + "_" + std::to_string(job_index) + ".txt";
        std::ofstream file(filename, std::ios::app);
        if (file.is_open()) {
            file << std::fixed << std::setprecision(15) << result << std::endl;
        }

        std::string block_std_filename = std::string(save_path) + "/output_beta_mode_" + std::to_string(MCMC_mode) +
                                         "_order_" + std::to_string(order_num) + "_block_std_" + std::to_string(job_index) + ".txt";
        std::ofstream block_std_file(block_std_filename);
        if (block_std_file.is_open()) {
            block_std_file << std::fixed << std::setprecision(15);
            for (size_t i = 0; i < ref_block_variances_dd[0].size(); ++i) {
                block_std_file << "0,O1," << DD_SQRT(ref_block_variances_dd[0][i]).to_double() << std::endl;
            }
        }

        // Save final sigma
        std::string sigma_filename = std::string(save_path) + "/output_beta_mode_" + std::to_string(MCMC_mode) +
                                    "_order_" + std::to_string(order_num) + "_sigma_" + std::to_string(job_index) + ".txt";
        std::ofstream sigma_file(sigma_filename);
        if (sigma_file.is_open()) {
            sigma_file << std::fixed << std::setprecision(15) << final_sigma << std::endl;
        }

        // Save accept ratio
        std::string accept_filename = std::string(save_path) + "/output_beta_mode_" + std::to_string(MCMC_mode) +
                                     "_order_" + std::to_string(order_num) + "_accept_ratio_" + std::to_string(job_index) + ".txt";
        std::ofstream accept_file(accept_filename);
        if (accept_file.is_open()) {
            accept_file << std::fixed << std::setprecision(15) << accept_ratio << std::endl;
        }
    }
    // Direct Monte Carlo integration for the final configuration
    else if (MCMC_mode == num_grid_chains + num_beta_chains + 1) {
        // Save seed
        std::string seed_filename = std::string(save_path) + "/output_MC_order_" + std::to_string(order_num) +
                                   "_seed_" + std::to_string(job_index) + ".txt";
        std::ofstream seed_file(seed_filename);
        if (seed_file.is_open()) {
            seed_file << seed << std::endl;
            seed_file.close();
        }

        const int dim = 3 * order_num;
        unsigned long numSamples = 1 << MC_bits;
        DoubleDouble sum(0);

        DD1DVector x(dim);
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        int largest_grid_idx = num_grid_chains;
        int largest_beta_idx = num_beta_chains;

        for (unsigned long i = 0; i < numSamples; i++) {
            for (int j = 0; j < dim; j++) {
                x[j] = DoubleDouble(dist(gen));
            }
            sum += DD_ABS(integrand_dd_sublattice_summed(x.data(), order_num, largest_grid_idx, largest_beta_idx));
        }

        DoubleDouble integral_dd = sum / DoubleDouble(numSamples);
        Real integral = integral_dd.to_double();

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
// 2. NEW INTEGRAND SPEED TEST
//================================================================================

void run_integrand_speed_test_with_loaded_data() {
    std::cout << "\n--- Speed Test: integrand_dd_sublattice_summed (using loaded data) ---\n";
    std::cout << "Testing for i=0, j=0\n";
    std::cout << " Order |      Avg Time/Run |   Runs/Sec" << std::endl;
    std::cout << "-------|-------------------|-----------" << std::endl;

    std::mt19937 rng(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    DoubleDouble result_sink(0.0); // To prevent compiler optimization

    // Loop through the requested order numbers
    for (int order_num = 4; order_num <= 9; ++order_num) {
        const int dim = 3 * order_num;
        std::vector<DoubleDouble> args(dim);
        long long runs = 0;
        const auto benchmark_duration = std::chrono::seconds(10); // Run for 2 seconds per order

        auto start_time = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_time;

        do {
            // Generate new random args for each call
            for (int k = 0; k < dim; ++k) {
                args[k] = DoubleDouble(dist(rng));
            }

            // Call the target function with loaded global data
            result_sink += integrand_dd_sublattice_summed(args.data(), order_num, 0, 0);
            runs++;
            elapsed_time = std::chrono::steady_clock::now() - start_time;
        } while (elapsed_time < benchmark_duration);

        double total_seconds = elapsed_time.count();
        double time_per_run = total_seconds / runs;
        double runs_per_second = runs / total_seconds;

        std::cout << "result_sink: " << result_sink.to_double() << std::endl;

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

    // Load MP_Real parameters
    for (const auto &entry : params["BETA"]) {
        beta_vec.push_back(string_to_mp_real(entry));
    }

    for (const auto &entry : params["BETA_h"]) {
        beta_h_vec.push_back(string_to_mp_real(entry));
    }

    for (const auto &entry : params["BETA_l"]) {
        beta_l_vec.push_back(string_to_mp_real(entry));
    }

    cheby_degree = std::stoi(params["CHEBY_DEGREE"][0]);
    cores_cpp = std::stoi(params["CORES_CPP"][0]);
    tol_fpm = string_to_mp_real(params["TOL_FPM"][0]);
    eps = string_to_mp_real(params["EPS"][0]);
    alpha_shift_up = string_to_mp_real(params["ALPHA_UP"][0]);
    alpha_shift_down = string_to_mp_real(params["ALPHA_DOWN"][0]);
    Delta_up = string_to_mp_real(params["DELTA_UP"][0]);
    Delta_down = string_to_mp_real(params["DELTA_DOWN"][0]);
    mu = string_to_mp_real(params["MU"][0]);
    mu_h = string_to_mp_real(params["MU_h"][0]);
    mu_l = string_to_mp_real(params["MU_l"][0]);

    U = string_to_mp_real(params["U"][0]);
    U_h = string_to_mp_real(params["U_h"][0]);
    U_l = string_to_mp_real(params["U_l"][0]);

    H = string_to_mp_real(params["H"][0]);
    H_h = string_to_mp_real(params["H_h"][0]);
    H_l = string_to_mp_real(params["H_l"][0]);

    lnZ_mode = stringToBool(params["LNZ_MODE"][0]);
    density_mode = stringToBool(params["DENSITY_MODE"][0]);
    double_occupancy_mode = stringToBool(params["DOUBLE_OCCUPANCY_MODE"][0]);
    energy_mode = stringToBool(params["ENERGY_MODE"][0]);
    staggered_magnetization_mode = stringToBool(params["STAGGERED_MAGNETIZATION_MODE"][0]);
    compressibility_mode = stringToBool(params["COMPRESSIBILITY_MODE"][0]);

    density_direct_mode = stringToBool(params["DENSITY_DIRECT_MODE"][0]);
    double_occupancy_direct_mode = stringToBool(params["DOUBLE_OCCUPANCY_DIRECT_MODE"][0]);

    // Read parameters for all orders 1-9
    for (int order = 1; order <= 9; ++order) {
        sigma_vec[order] = string_to_mp_real(params["SIGMA_ORDER_" + std::to_string(order)][0]);
        iter_bits_vec[order] = std::stoi(params["ITER_BITS_ORDER_" + std::to_string(order)][0]);
        warm_bits_vec[order] = std::stoi(params["WARM_BITS_ORDER_" + std::to_string(order)][0]);
        MC_bits_vec[order] = std::stoi(params["MC_BITS_ORDER_" + std::to_string(order)][0]);
    }

    seed_offset = std::stoul(params["SEED_OFFSET"][0]);
    seed_interval = std::stoul(params["SEED_INTERVAL"][0]);
    MCMC_mode = std::stoi(params["MCMC_mode"][0]);
    MCMC_cut = string_to_mp_real(params["MCMC_CUT"][0]);
    greens_func_exp_cut = string_to_mp_real(params["GREENS_FUNC_EXP_CUT"][0]);

    MSB_table = createMSBTable();

    // Load data with MP_Real precision
    std::cout << "Loading data with MP_Real precision (100 decimal places)...\n";
    load_data();

    // Convert all MP_Real flattened vectors to DoubleDouble
    std::cout << "Converting data to DoubleDouble for fast computation...\n";
    convert_all_mp_to_dd();

    std::cout << "\n--- Loaded Parameters ---\n";
    std::cout << "Using MP_Real with 100 decimal places for data preparation\n";
    std::cout << "Using DoubleDouble for fast computation\n";

    // Print parameters
    std::cout << "cheby_degree: " << cheby_degree << "\n";
    std::cout << "cores_cpp: " << cores_cpp << "\n";
    std::cout << "tol_fpm: " << tol_fpm << "\n";
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

    // 3. RUN THE NEW SPEED TEST USING THE LOADED DATA
    run_integrand_speed_test_with_loaded_data();

    // 4. Continue with the rest of your main application logic
    std::cout << "\n--- Continuing with MCMC simulation ---\n";

    int threads_per_job = cores_cpp;

    // Time the parallel execution of std::async
    auto start_time = std::chrono::high_resolution_clock::now();

    // Create a vector to store the futures
    std::vector<std::future<void>> futures;

    // Launch parallel threads
    for (int thread_num = 1; thread_num <= threads_per_job; ++thread_num) {
        // Calculate the actual run index for this thread
        unsigned long int run_index = thread_num;
        // Launch the thread and store its future
        futures.push_back(std::async(std::launch::async, [run_index]() {
            run_all_dd(run_index);
        }));
    }

    // Wait for all threads to complete and propagate exceptions
    for (auto& future : futures) {
        try {
            future.get();  // This will rethrow any exceptions
        } catch (const std::exception& e) {
            std::cerr << "Thread exception: " << e.what() << std::endl;
            throw;  // Re-throw to maintain original behavior
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

    // Separate vectors for sign files (no order number) with job indices
    std::vector<std::pair<std::string, int>> sign_files;
    std::vector<std::pair<std::string, int>> sign_block_std_files;
    std::vector<std::pair<std::string, int>> sign_seed_files;
    std::vector<std::pair<std::string, int>> sign_sigma_files;
    std::vector<std::pair<std::string, int>> sign_accept_ratio_files;

    // Regular expressions to match file patterns and extract order numbers and job indices
    std::regex sign_pattern(R"(output_sign_(\d+)\.txt)");
    std::regex sign_block_pattern(R"(output_sign_block_std_(\d+)\.txt)");
    std::regex sign_seed_pattern(R"(output_sign_seed_(\d+)\.txt)");
    std::regex sign_sigma_pattern(R"(output_sign_sigma_(\d+)\.txt)");
    std::regex sign_accept_pattern(R"(output_sign_accept_ratio_(\d+)\.txt)");

    std::regex grid_pattern(R"(output_grid_mode_\d+_order_(\d+)_(\d+)\.txt)");
    std::regex grid_block_pattern(R"(output_grid_mode_\d+_order_(\d+)_block_std_(\d+)\.txt)");
    std::regex grid_seed_pattern(R"(output_grid_mode_\d+_order_(\d+)_seed_(\d+)\.txt)");
    std::regex grid_sigma_pattern(R"(output_grid_mode_\d+_order_(\d+)_sigma_(\d+)\.txt)");
    std::regex grid_accept_pattern(R"(output_grid_mode_\d+_order_(\d+)_accept_ratio_(\d+)\.txt)");

    std::regex beta_pattern(R"(output_beta_mode_\d+_order_(\d+)_(\d+)\.txt)");
    std::regex beta_block_pattern(R"(output_beta_mode_\d+_order_(\d+)_block_std_(\d+)\.txt)");
    std::regex beta_seed_pattern(R"(output_beta_mode_\d+_order_(\d+)_seed_(\d+)\.txt)");
    std::regex beta_sigma_pattern(R"(output_beta_mode_\d+_order_(\d+)_sigma_(\d+)\.txt)");
    std::regex beta_accept_pattern(R"(output_beta_mode_\d+_order_(\d+)_accept_ratio_(\d+)\.txt)");

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
                if (is_seed_file) {
                    seed_files_by_order[order_num].push_back({full_path, job_index});
                } else if (is_block_std) {
                    block_std_files_by_order[order_num].push_back({full_path, job_index});
                } else if (is_sigma_file) {
                    sigma_files_by_order[order_num].push_back({full_path, job_index});
                } else if (is_accept_file) {
                    accept_ratio_files_by_order[order_num].push_back({full_path, job_index});
                } else {
                    regular_files_by_order[order_num].push_back({full_path, job_index});
                }
            }
        }
    }

    // Lambda to sort files by job index
    auto sort_by_job_index = [](std::vector<std::pair<std::string, int>>& files) {
        std::sort(files.begin(), files.end(),
                  [](const auto& a, const auto& b) { return a.second < b.second; });
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