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

// -------------------- Global knobs --------------------
inline int cheby_degree;
inline int cores_cpp;

inline int MCMC_mode;
inline int set_num = 0;
inline Real greens_func_exp_cut;
inline Real MCMC_cut;

inline const char* main_path = "/users";
inline const char* save_path = "/users";

// -------------------- Lattice + temperature (scalars) --------------------
inline Int1DVector N1_vec;
inline Int1DVector N2_vec;
inline Real1DVector beta_vec;

// -------------------- Flavours --------------------
inline int N_f;                  // must be set before building DAGs
inline Real1DVector alpha_shift_vec; // size N_f
inline Real1DVector mu_vec;
inline Real1DVector mu_h_vec;
inline Real1DVector mu_l_vec;

// -------------------- Interaction --------------------
inline Real U;

// -------------------- MC scheduling --------------------
inline std::vector<Real> sigma_vec;  // resized to 10 in main
inline Real eps;

inline int N_MIN_ORDER;
inline int N_MAX_ORDER;

inline bool lnZ_mode;
inline bool density_mode;

inline unsigned long seed_offset;
inline unsigned long seed_interval;

inline Int1DVector iter_bits_vec;  // resized to 10 in main
inline Int1DVector warm_bits_vec;  // resized to 10 in main
inline Int1DVector MC_bits_vec;    // resized to 10 in main

// -------------------- Green's data (multi-flavour) --------------------
// g_vec:    [flavour][cheby_degree+1][N1][N2]
inline std::unordered_map<std::string, Real4DVector> g_vec;
inline std::unordered_map<std::string, Real4DVector> g_vec_mu_h;
inline std::unordered_map<std::string, Real4DVector> g_vec_mu_l;

// Flattened versions: [flavour][(cheby_degree+1)*N1*N2]
inline std::unordered_map<std::string, Real2DVector> g_vec_flattened;
inline std::unordered_map<std::string, Real2DVector> g_vec_mu_h_flattened;
inline std::unordered_map<std::string, Real2DVector> g_vec_mu_l_flattened;

// Helper: flatten (cheby_degree+1) x N1 x N2 -> 1D row-major
inline void flatten_3D_to_1D(const Real3DVector& vec3D, Real1DVector& vec1D) {
    const size_t D0 = vec3D.size();             // cheby_degree+1
    const size_t D1 = D0 ? vec3D[0].size() : 0; // N1
    const size_t D2 = (D1 ? vec3D[0][0].size() : 0); // N2
    vec1D.clear();
    vec1D.reserve(D0 * D1 * D2);
    for (size_t i = 0; i < D0; ++i)
        for (size_t j = 0; j < D1; ++j)
            for (size_t k = 0; k < D2; ++k)
                vec1D.push_back(vec3D[i][j][k]);
}

// Initialize vectors with proper sizes
inline void initialize_parameter_vectors() {
    sigma_vec.resize(10);     // index 0 unused, 1-9 for orders 1-9
    iter_bits_vec.resize(10);
    warm_bits_vec.resize(10);
    MC_bits_vec.resize(10);
}

// -------------------- Typed DAG cache --------------------
inline std::unordered_map<int, typed_dag::DAGData> dag_cache;

inline void initialize_all_dags() {
    std::cout << "Building typed DAGs (N_f = " << N_f << ") for orders 2..7..." << std::endl;
    dag_cache.clear();
    for (int order = 2; order <= 7; ++order) {
        dag_cache[order] = typed_dag::build_dag_typed(order, N_f);
        std::cout << "  Order " << order << " DAG built with sink node "
                  << dag_cache[order].sink_node_id
                  << " and " << dag_cache[order].layer_edges.size()
                  << " layers\n";
    }
}

#endif
