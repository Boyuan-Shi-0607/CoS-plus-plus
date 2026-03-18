#ifndef RUN_ALL_ORDER_UPDATE_MCMC_H
#define RUN_ALL_ORDER_UPDATE_MCMC_H

#include "order_update_MCMC.h"
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

// ---------------------- your integrand wrapper (i=j=0) ----------------------
// Calls through to your connected / density modes via your existing globals.
// (i,j) indices are fixed to (0,0).
inline Real integrand_order_update(const Real* args, const int order_num) {
    const typed_dag::DAGData* dag = (order_num > 1 && dag_cache.find(order_num) != dag_cache.end())
                                    ? &dag_cache[order_num] : nullptr;

    if (lnZ_mode) {
        return integrand_connected_typed(
            beta_vec[0], U, args,
            g_vec_flattened["0_0"],
            order_num, N1_vec[0], N2_vec[0],
            eps, cheby_degree,
            alpha_shift_vec,
            dag
        );
    }
    if (density_mode) {
        Real result_mu_h = integrand_connected_typed(
            beta_vec[0], U, args,
            g_vec_mu_h_flattened["0_0"],
            order_num, N1_vec[0], N2_vec[0],
            eps, cheby_degree,
            alpha_shift_vec,
            dag
        );
        Real result_mu_l = integrand_connected_typed(
            beta_vec[0], U, args,
            g_vec_mu_l_flattened["0_0"],
            order_num, N1_vec[0], N2_vec[0],
            eps, cheby_degree,
            alpha_shift_vec,
            dag
        );
        return (result_mu_h - result_mu_l) / (beta_vec[0] * (mu_h_vec[0] - mu_l_vec[0]));
    }
    throw std::runtime_error("Invalid mode: set either LNZ_MODE or DENSITY_MODE in params.txt");
}

inline void run_all_order_update_MCMC(int run_index) {

    // ---------------------- RNG seed ----------------------
    std::random_device rd;
    std::mt19937 gen_seed(rd());
    std::uniform_int_distribution<unsigned long int> dist_seed(0, 100000000000000ULL);

    unsigned long int seed = dist_seed(gen_seed);
    std::mt19937 gen(seed);

    std::cout << "Seed: " << seed << "\n";

    // ---------------------- Bind your integrand (i=j=0) ----------------------
    ORDER_UPDATE_MCMC::RealFunc f = [](const Real* x, int n) {
        return integrand_order_update(x, n); // (i,j) fixed to (0,0)
    };

    const int order_min = N_MIN_ORDER;
    const int order_max = N_MAX_ORDER;
    if (order_max < order_min) {
        std::cerr << "Invalid order range: N_MAX_ORDER < N_MIN_ORDER\n";
        return;
    }
    const int num_orders = order_max - order_min + 1;

    // ---------------------- Main Variable-order MCMC ----------------------
    unsigned WARM_BITS = 20;
    unsigned ITER_BITS = 24;
    const Real SIGMA_INIT = 0.25;

    auto t3 = std::chrono::high_resolution_clock::now();
    ORDER_UPDATE_MCMC::RunSummary S = ORDER_UPDATE_MCMC::run_sign_variable_order(
        f,
        /*initial_sigma*/    SIGMA_INIT,
        /*iter_bits*/        ITER_BITS,
        /*warm_bits*/        WARM_BITS,
        /*gen*/              gen
    );
    auto t4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dt_mcmc = t4 - t3;

    // ---------------------- Report MCMC ----------------------
    auto ratio = [](unsigned long long att, unsigned long long acc){
        return (att==0ULL)?0.0: (double)acc/(double)att;
    };

    std::cout << "\n=== Variable-order MCMC (n ∈ [" << order_min << "," << order_max << "]) — last-vertex add/remove ===\n";
    std::cout << "prod steps = " << (1ULL<<ITER_BITS) << ", warmup = " << (1ULL<<WARM_BITS) << "\n";
    std::cout << "final sigma = " << S.final_sigma << "\n";
    std::cout << "overall acceptance = " << S.accept_ratio << "\n";
    std::cout << "accept(within/add/rem) = "
         << ratio(S.prod_attempt_within,S.prod_accept_within) << " / "
         << ratio(S.prod_attempt_add   ,S.prod_accept_add   ) << " / "
         << ratio(S.prod_attempt_remove,S.prod_accept_remove) << "\n";
    std::cout << "elapsed (MCMC) = " << dt_mcmc.count() << " s\n";

    // Per-order mean sign & visits
    std::cout << "\n--- Per-order mean sign (per TOTAL steps) ---\n";
    std::vector<Real> O_mcmc(num_orders, 0.0);
    std::vector<unsigned long long> visits_vec(num_orders, 0ULL);

    for (int n = order_min; n <= order_max; ++n) {
        Real mean_sign = 0.0;
        unsigned long long visits = 0ULL;

        if ((int)S.sign_mean_per_order.size() > n)  mean_sign = S.sign_mean_per_order[n];
        if ((int)S.sign_count_per_order.size() > n) visits    = S.sign_count_per_order[n];

        O_mcmc[n - order_min] = mean_sign;
        visits_vec[n - order_min] = visits;
        std::cout << "n=" << n << ": mean_sign=" << mean_sign << "  visits=" << visits << "\n";
    }

    // ---------------------- Save all results to text files ----------------------

    // 1. Save visits at all orders
    {
        std::string visits_filename = std::string(save_path) + "/visits_run_" + std::to_string(run_index) + ".txt";
        std::ofstream visits_file(visits_filename);
        if (visits_file.is_open()) {
            visits_file << std::fixed << std::setprecision(15);
            visits_file << "# Order Visits\n";
            for (int i = 0; i < num_orders; ++i) {
                int n = order_min + i;
                visits_file << n << " " << visits_vec[i] << "\n";
            }
            visits_file.close();
        }
    }

    // 2. Save signs at all orders
    {
        std::string signs_filename = std::string(save_path) + "/signs_run_" + std::to_string(run_index) + ".txt";
        std::ofstream signs_file(signs_filename);
        if (signs_file.is_open()) {
            signs_file << std::fixed << std::setprecision(15);
            signs_file << "# Order MeanSign\n";
            for (int i = 0; i < num_orders; ++i) {
                int n = order_min + i;
                signs_file << n << " " << O_mcmc[i] << "\n";
            }
            signs_file.close();
        }
    }

    // 3. Save seed
    {
        std::string seed_filename = std::string(save_path) + "/seed_run_" + std::to_string(run_index) + ".txt";
        std::ofstream seed_file(seed_filename);
        if (seed_file.is_open()) {
            seed_file << seed << "\n";
            seed_file.close();
        }
    }

    // 4. Save acceptance ratios
    {
        std::string accept_filename = std::string(save_path) + "/accept_ratios_run_" + std::to_string(run_index) + ".txt";
        std::ofstream accept_file(accept_filename);
        if (accept_file.is_open()) {
            accept_file << std::fixed << std::setprecision(15);
            accept_file << "# Acceptance ratios\n";
            accept_file << "overall_acceptance " << S.accept_ratio << "\n";
            accept_file << "accept_within " << ratio(S.prod_attempt_within, S.prod_accept_within) << "\n";
            accept_file << "accept_add " << ratio(S.prod_attempt_add, S.prod_accept_add) << "\n";
            accept_file << "accept_remove " << ratio(S.prod_attempt_remove, S.prod_accept_remove) << "\n";
            accept_file << "# Raw counts\n";
            accept_file << "prod_attempt_within " << S.prod_attempt_within << "\n";
            accept_file << "prod_accept_within " << S.prod_accept_within << "\n";
            accept_file << "prod_attempt_add " << S.prod_attempt_add << "\n";
            accept_file << "prod_accept_add " << S.prod_accept_add << "\n";
            accept_file << "prod_attempt_remove " << S.prod_attempt_remove << "\n";
            accept_file << "prod_accept_remove " << S.prod_accept_remove << "\n";
            accept_file.close();
        }
    }

    // 5. Save P_ADD and P_REM
    {
        std::string params_filename = std::string(save_path) + "/params_run_" + std::to_string(run_index) + ".txt";
        std::ofstream params_file(params_filename);
        if (params_file.is_open()) {
            params_file << std::fixed << std::setprecision(15);
            params_file << "# MCMC Parameters\n";
            params_file << "P_ADD " << ORDER_UPDATE_MCMC::P_ADD << "\n";
            params_file << "P_REM " << ORDER_UPDATE_MCMC::P_REM << "\n";
            params_file << "P_WITHIN " << (1.0 - ORDER_UPDATE_MCMC::P_ADD - ORDER_UPDATE_MCMC::P_REM) << "\n";
            params_file.close();
        }
    }

    // 6. Save final_sigma
    {
        std::string sigma_filename = std::string(save_path) + "/sigma_run_" + std::to_string(run_index) + ".txt";
        std::ofstream sigma_file(sigma_filename);
        if (sigma_file.is_open()) {
            sigma_file << std::fixed << std::setprecision(15);
            sigma_file << "# Initial and final sigma\n";
            sigma_file << "sigma_init " << SIGMA_INIT << "\n";
            sigma_file << "final_sigma " << S.final_sigma << "\n";
            sigma_file.close();
        }
    }
}
#endif //RUN_ALL_ORDER_UPDATE_MCMC_H
