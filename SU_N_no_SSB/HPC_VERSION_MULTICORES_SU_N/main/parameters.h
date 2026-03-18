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
#include "dag_generator.h"

inline int MCMC_mode;

inline Int1DVector N1_vec;
inline Int1DVector N2_vec;

// Use Real for all data
inline std::vector<Real> beta_vec;

inline int cheby_degree;
inline int cores_cpp;

inline Real greens_func_exp_cut;
inline Real MCMC_cut;

inline const char* main_path = "/users";
inline const char* save_path = "/users";

// Single alpha shift for honeycomb lattice
inline Real alpha_shift;

// chemical potential
inline Real mu;
inline Real mu_h;
inline Real mu_l;

// interaction
inline Real U;

// Changed from order_num to set_num
inline int set_num = 0;

// Store sigma values for different orders (up to order 9)
inline std::vector<Real> sigma_vec;  // Will be resized to 10 in main
inline Real eps;
inline unsigned long long seed;

inline bool lnZ_mode;
inline bool density_mode;

inline unsigned long seed_offset;
inline unsigned long seed_interval;

// Store iter_bits and warm_bits for different orders
inline Int1DVector iter_bits_vec;  // Will be resized to 10 in main
inline Int1DVector warm_bits_vec;  // Will be resized to 10 in main
inline Int1DVector MC_bits_vec;    // Will be resized to 10 in main

// Use Real for data storage - now 3D vectors for scalar Green's functions
inline std::unordered_map<std::string, Real3DVector> g_vec;
inline std::unordered_map<std::string, Real3DVector> g_vec_mu_h;
inline std::unordered_map<std::string, Real3DVector> g_vec_mu_l;

// Flattened versions for fast access
inline std::unordered_map<std::string, Real1DVector> g_vec_flattened;
inline std::unordered_map<std::string, Real1DVector> g_vec_mu_h_flattened;
inline std::unordered_map<std::string, Real1DVector> g_vec_mu_l_flattened;

// Pre-built DAGs for each order (built once at startup)
inline std::unordered_map<int, DAGData> dag_cache;

// Helper function to flatten 3D to 1D
inline void flatten_3D_to_1D(const Real3DVector& vec3D, Real1DVector& vec1D) {
    vec1D.clear();
    for (size_t i = 0; i < vec3D.size(); ++i) {
        for (size_t j = 0; j < vec3D[i].size(); ++j) {
            for (size_t k = 0; k < vec3D[i][j].size(); ++k) {
                vec1D.push_back(vec3D[i][j][k]);
            }
        }
    }
}

// Initialize vectors with proper sizes
inline void initialize_parameter_vectors() {
    sigma_vec.resize(10);  // index 0 unused, 1-9 for orders 1-9
    iter_bits_vec.resize(10);  // index 0 unused, 1-9 for orders 1-9
    warm_bits_vec.resize(10);  // index 0 unused, 1-9 for orders 1-9
    MC_bits_vec.resize(10);    // index 0 unused, 1-9 for orders 1-9
}

// Build all DAGs at startup
inline void initialize_all_dags() {
    std::cout << "Building DAGs for all orders..." << std::endl;
    for (int order = 2; order <= 9; ++order) {
        dag_cache[order] = build_dag(order);
        std::cout << "  Order " << order << " DAG built with sink node "
                  << dag_cache[order].sink_node_id
                  << " and " << dag_cache[order].layer_edges.size()
                  << " layers" << std::endl;
    }
}

#endif