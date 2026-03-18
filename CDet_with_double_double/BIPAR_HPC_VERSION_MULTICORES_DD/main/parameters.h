#ifndef BIPARTITE_LATTICE_PRINCIPAL_MINOR_CDET_PARAMETERS_H
#define BIPARTITE_LATTICE_PRINCIPAL_MINOR_CDET_PARAMETERS_H

#include <vector>
#include <unordered_map>
#include <iostream>
#include <array>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "../utility/types.h"

inline int MCMC_mode;

inline Int1DVector N1_vec;
inline Int1DVector N2_vec;

// Use MP_Real for high precision data storage
inline std::vector<MP_Real> beta_vec;
inline std::vector<MP_Real> beta_h_vec;
inline std::vector<MP_Real> beta_l_vec;

inline int cheby_degree;
inline int cores_cpp;
inline MP_Real tol_fpm;

inline MP_Real greens_func_exp_cut;
inline MP_Real MCMC_cut;

inline const char* main_path = "/users";
inline const char* save_path = "/users";

inline MP_Real alpha_shift_up;
inline MP_Real alpha_shift_down;
inline MP_Real Delta_up;
inline MP_Real Delta_down;

// chemical potential
inline MP_Real mu;
inline MP_Real mu_h;
inline MP_Real mu_l;

// interaction
inline MP_Real U;
inline MP_Real U_h;
inline MP_Real U_l;

// magnetic fields
inline MP_Real H;
inline MP_Real H_h;
inline MP_Real H_l;

// Changed from order_num to set_num
inline int set_num = 0;

// Store sigma values for different orders (up to order 9)
inline std::vector<MP_Real> sigma_vec;  // Will be resized to 10 in main
inline MP_Real eps;
inline unsigned long long seed;
inline unsigned int n_swap_rounds;

inline bool lnZ_mode;
inline bool double_occupancy_mode;
inline bool staggered_magnetization_mode;
inline bool density_mode;
inline bool energy_mode;
inline bool compressibility_mode;

// New direct computation modes
inline bool density_direct_mode;
inline bool double_occupancy_direct_mode;

inline unsigned long seed_offset;
inline unsigned long seed_interval;

// Store iter_bits and warm_bits for different orders
inline Int1DVector iter_bits_vec;  // Will be resized to 10 in main
inline Int1DVector warm_bits_vec;  // Will be resized to 10 in main
inline Int1DVector MC_bits_vec;    // Will be resized to 10 in main

// Use MP_Real for data storage
inline std::unordered_map<std::string, MP_Real5DVector> g_vec;
inline std::unordered_map<std::string, MP_Real5DVector> g_vec_mu_h;
inline std::unordered_map<std::string, MP_Real5DVector> g_vec_mu_l;
inline std::unordered_map<std::string, MP_Real5DVector> g_vec_beta_h;
inline std::unordered_map<std::string, MP_Real5DVector> g_vec_beta_l;
inline std::unordered_map<std::string, MP_Real5DVector> g_vec_H_h;
inline std::unordered_map<std::string, MP_Real5DVector> g_vec_H_l;

// Store both MP_Real and DoubleDouble versions
inline std::unordered_map<std::string, MP_Real1DVector> g_vec_flattened_mp;
inline std::unordered_map<std::string, MP_Real1DVector> g_vec_mu_h_flattened_mp;
inline std::unordered_map<std::string, MP_Real1DVector> g_vec_mu_l_flattened_mp;
inline std::unordered_map<std::string, MP_Real1DVector> g_vec_beta_h_flattened_mp;
inline std::unordered_map<std::string, MP_Real1DVector> g_vec_beta_l_flattened_mp;
inline std::unordered_map<std::string, MP_Real1DVector> g_vec_H_h_flattened_mp;
inline std::unordered_map<std::string, MP_Real1DVector> g_vec_H_l_flattened_mp;

// DoubleDouble versions for fast computation
inline std::unordered_map<std::string, DD1DVector> g_vec_flattened;
inline std::unordered_map<std::string, DD1DVector> g_vec_mu_h_flattened;
inline std::unordered_map<std::string, DD1DVector> g_vec_mu_l_flattened;
inline std::unordered_map<std::string, DD1DVector> g_vec_beta_h_flattened;
inline std::unordered_map<std::string, DD1DVector> g_vec_beta_l_flattened;
inline std::unordered_map<std::string, DD1DVector> g_vec_H_h_flattened;
inline std::unordered_map<std::string, DD1DVector> g_vec_H_l_flattened;

inline std::array<unsigned char, 256> MSB_table;

// Conversion helper functions
inline DD1DVector convert_mp_to_dd_vector(const MP_Real1DVector& mp_vec) {
    DD1DVector dd_vec;
    dd_vec.reserve(mp_vec.size());
    for (const auto& val : mp_vec) {
        dd_vec.push_back(DoubleDouble(val));
    }
    return dd_vec;
}

inline void convert_all_mp_to_dd() {
    // Convert all MP_Real flattened vectors to DoubleDouble
    for (const auto& [key, mp_vec] : g_vec_flattened_mp) {
        g_vec_flattened[key] = convert_mp_to_dd_vector(mp_vec);
    }
    for (const auto& [key, mp_vec] : g_vec_mu_h_flattened_mp) {
        g_vec_mu_h_flattened[key] = convert_mp_to_dd_vector(mp_vec);
    }
    for (const auto& [key, mp_vec] : g_vec_mu_l_flattened_mp) {
        g_vec_mu_l_flattened[key] = convert_mp_to_dd_vector(mp_vec);
    }
    for (const auto& [key, mp_vec] : g_vec_beta_h_flattened_mp) {
        g_vec_beta_h_flattened[key] = convert_mp_to_dd_vector(mp_vec);
    }
    for (const auto& [key, mp_vec] : g_vec_beta_l_flattened_mp) {
        g_vec_beta_l_flattened[key] = convert_mp_to_dd_vector(mp_vec);
    }
    for (const auto& [key, mp_vec] : g_vec_H_h_flattened_mp) {
        g_vec_H_h_flattened[key] = convert_mp_to_dd_vector(mp_vec);
    }
    for (const auto& [key, mp_vec] : g_vec_H_l_flattened_mp) {
        g_vec_H_l_flattened[key] = convert_mp_to_dd_vector(mp_vec);
    }
}

// Initialize vectors with proper sizes
inline void initialize_parameter_vectors() {
    sigma_vec.resize(10);  // index 0 unused, 1-9 for orders 1-9
    iter_bits_vec.resize(10);  // index 0 unused, 1-9 for orders 1-9
    warm_bits_vec.resize(10);  // index 0 unused, 1-9 for orders 1-9
    MC_bits_vec.resize(10);    // index 0 unused, 1-9 for orders 1-9
}

#endif