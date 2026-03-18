// main.cpp – Variable-order MCMC + Plain Monte Carlo (arbitrary inclusive order range)
// - Uses MC::run_sign_variable_order from MCMC_vegas.h (respects N_MIN_ORDER..N_MAX_ORDER)
// - Binds MC::RealFunc f(x,n) -> integrand_connected(..., order_num, N1_vec[0], N2_vec[0], ...)
//   with (i,j) fixed to (0,0) using your global vectors/data.
// - Prints Plain MC baselines for every order n in [N_MIN_ORDER, N_MAX_ORDER]:
//     I_n, A_n, stderrs, per-order runtime
//   Then forms ΣA = sum_k A_k over the whole range and reports
//     O_n(plain) = I_n / ΣA for each n.
// - Runs variable-order MCMC once and reports acceptance stats, per-order mean sign and visits,
//   and compares O_n(MCMC) = S.sign_mean_per_order[n] vs O_n(plain) for every n.
//
// NOTE: Adjust ITER_BITS/WARM_BITS/MC_SAMPLES_PER_ORDER if runtime is long.

#include "../utility/types.h"
#include "integrand.h"
#include "load_data.h"
#include "../utility/check_range.h"
#include "lnZ_free_part.h"
#include "dag_generator.h"
#include "file_combination.h"
#include "lnZ_free_part.h"
#include "run_all_order_update_MCMC.h"
#include "run_all_RG_MCMC.h"

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

#define USE_RG_MCMC 0

// ---------------------- small helpers ----------------------
static inline Real string_to_real(const std::string& s) { return std::stod(s); }

static inline void trim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch){ return !std::isspace(ch); }));
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch){ return !std::isspace(ch); }).base(), s.end());
}
static inline bool stringToBool(const std::string &str) {
    std::string t = str;
    std::transform(t.begin(), t.end(), t.begin(), [](unsigned char c){ return (char)std::tolower(c); });
    if (t=="true"||t=="1") return true;
    if (t=="false"||t=="0") return false;
    throw std::invalid_argument("Invalid bool string: " + str);
}


// ---------------------- main ----------------------
int main() {
    // run_all();
    std::cout.setf(std::ios::fixed);
    std::cout << std::setprecision(15);
    std::cout << "\n--- Variable-order MCMC + Plain MC using integrand(args, n, 0, 0) ---\n";
    std::cout << "Orders: inclusive range [" << N_MIN_ORDER << ", " << N_MAX_ORDER << "]\n";

    // Basic safety checks / vector sizing
    check_range();
    initialize_parameter_vectors();

    // ---------------------- read params.txt ----------------------
    std::unordered_map<std::string, std::vector<std::string>> params;
    std::ifstream infile("../params.txt");
    if (!infile) {
        std::cerr << "Could not open params.txt\n";
    }
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        size_t pos = line.find('=');
        if (pos == std::string::npos) continue;
        std::string key = line.substr(0, pos);
        std::string val = line.substr(pos + 1);
        trim(key); trim(val);
        std::vector<std::string> toks;
        std::stringstream ss(val);
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            trim(tok);
            if (!tok.empty()) toks.push_back(tok);
        }
        params[key] = toks;
    }

    N_f = std::stoi(params["N_f"][0]);
    N_MIN_ORDER = std::stoi(params["N_MIN_ORDER"][0]);
    N_MAX_ORDER = std::stoi(params["N_MAX_ORDER"][0]);

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

    for (const auto &entry : params["MU"]) {
        mu_vec.push_back(string_to_real(entry));
    }
    for (const auto &entry : params["MU_h"]) {
        mu_h_vec.push_back(string_to_real(entry));
    }
    for (const auto &entry : params["MU_l"]) {
        mu_l_vec.push_back(string_to_real(entry));
    }

    cheby_degree = std::stoi(params["CHEBY_DEGREE"][0]);
    eps          = string_to_real(params["EPS"][0]);
    for (const auto& e : params["ALPHA_SHIFT_VEC"]) alpha_shift_vec.push_back(string_to_real(e));

    U            = string_to_real(params["U"][0]);
    lnZ_mode     = stringToBool(params["LNZ_MODE"][0]);
    density_mode = stringToBool(params["DENSITY_MODE"][0]);
    cores_cpp = std::stoi(params["CORES_CPP"][0]);

    // Optional knobs used in your broader code-paths
    greens_func_exp_cut = string_to_real(params["GREENS_FUNC_EXP_CUT"][0]);
    MCMC_cut           = string_to_real(params["MCMC_CUT"][0]);

    // Optional seeds (not strictly required here)
    seed_offset   = std::stoul(params["SEED_OFFSET"][0]);
    seed_interval = std::stoul(params["SEED_INTERVAL"][0]);
    MCMC_mode     = std::stoi(params["MCMC_mode"][0]);

    // Read parameters for all orders 1-9
    for (int order = 1; order <= 9; ++order) {
        sigma_vec[order] = string_to_real(params["SIGMA_ORDER_" + std::to_string(order)][0]);
        iter_bits_vec[order] = std::stoi(params["ITER_BITS_ORDER_" + std::to_string(order)][0]);
        warm_bits_vec[order] = std::stoi(params["WARM_BITS_ORDER_" + std::to_string(order)][0]);
        MC_bits_vec[order] = std::stoi(params["MC_BITS_ORDER_" + std::to_string(order)][0]);
    }

    // ---------------------- load heavy data and DAGs ----------------------
    std::cout << "Loading data...\n";
    load_data();
    initialize_all_dags();

    std::cout << "\n--- Loaded Parameters ---\n";
    std::cout << "Using Real (double) for all computations\n";

    // Print parameters
    std::cout << "cheby_degree: " << cheby_degree << "\n";
    std::cout << "cores_cpp: " << cores_cpp << "\n";
    std::cout << "eps: " << eps << "\n";
    // std::cout << "mu: " << mu << "\n";
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

    // std::cout << "lnZ free per site: " << calculate_lnZ_per_site(N1, N2, beta, mu, {2.3, 2.3, 2.3}) << std::endl;

    // 5. Continue with the rest of your main application logic
    std::cout << "\n--- Continuing with MCMC simulation ---\n";

#if USE_RG_MCMC
    int threads_per_job = cores_cpp;

    // Time the parallel execution of std::async
    auto start_time = std::chrono::high_resolution_clock::now();

    // Create a vector to store the futures
    std::vector<std::future<void>> futures;

    // Launch parallel threads
    for (int thread_num = 1; thread_num <= threads_per_job; ++thread_num) {
        unsigned long int run_index = thread_num;

        futures.push_back(std::async(std::launch::async, [run_index]() {
            run_all_RG_MCMC(run_index);
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
    combine_files_RG_MCMC(10);
#else
    // ---------------------- Main MCMC runs ----------------------
    int threads_per_job = cores_cpp;

    // Time the parallel execution of std::async
    auto start_time = std::chrono::high_resolution_clock::now();

    // FIX 2 & 3: Use different name for the second futures vector
    std::vector<std::future<void>> mcmc_futures;  // Renamed from 'futures' to 'mcmc_futures'

    // Launch parallel threads for main MCMC
    for (int thread_num = 1; thread_num <= threads_per_job; ++thread_num) {
        unsigned long int run_index = thread_num;

        mcmc_futures.push_back(std::async(std::launch::async, [run_index]() {
            run_all_order_update_MCMC(run_index);
        }));
    }

    // Wait for all MCMC threads to complete and propagate exceptions
    for (auto& future : mcmc_futures) {
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
    combine_order_update_files(10);
#endif

    return 0;
}