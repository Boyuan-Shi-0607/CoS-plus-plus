#ifndef DAG_GENERATOR_H
#define DAG_GENERATOR_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <iomanip> // For printing
#include "../utility/types.h"

#define RESTRICT __restrict__

// ===== STRUCTURES =====
struct RawEdge {
    int src;
    int dst;
    int row;
    int col;
    Real factor;
};

struct DAGData {
    std::vector<std::vector<RawEdge>> layer_edges;
    int sink_node_id;
    int max_node_id;
};

struct State {
    int u;
    int v;
    std::unordered_set<int> R;

    State() : u(0), v(0) {}
    State(int u_, int v_, const std::unordered_set<int>& R_) : u(u_), v(v_), R(R_) {}

    bool operator==(const State& other) const {
        return u == other.u && v == other.v && R == other.R;
    }
};

struct StateHash {
    size_t operator()(State const& s) const noexcept {
        size_t h = 0;
        auto hc = std::hash<int>{};
        h ^= hc(s.u) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= hc(s.v) + 0x9e3779b9 + (h << 6) + (h >> 2);
        std::vector<int> elems(s.R.begin(), s.R.end());
        std::sort(elems.begin(), elems.end());
        for (int r : elems) {
            h ^= hc(r) + 0x9e3779b9 + (h << 6) + (h >> 2);
        }
        return h;
    }
};

// ===== HELPER FUNCTIONS =====

// CORRECTED STATE NORMALIZATION
inline State shift_state(const State& s, int order_num) {
    const int n = order_num;

    // First, normalize the set R to its canonical form. This process is
    // common to all state types.
    std::vector<int> shifted_list;
    shifted_list.reserve(s.R.size());
    for (int r_val : s.R) {
        shifted_list.push_back((r_val - 1) % n + 1);
    }
    std::sort(shifted_list.begin(), shifted_list.end());

    std::unordered_set<int> result_set;
    std::unordered_map<int, int> counts;
    for (int val : shifted_list) {
        counts[val]++;
        if (counts[val] == 1) {
            result_set.insert(val);
        } else {
            result_set.insert(val + n);
        }
    }

    // Second, determine the canonical u and v, ensuring cycle-starters (u=v)
    // are never identified with paths between partners (u=partner(v)).
    int new_u, new_v;

    // Project original u and v to their top-level representatives
    int projected_u = (s.u - 1) % n + 1;
    int projected_v = (s.v - 1) % n + 1;

    if (s.u == s.v) {
        // True cycle-starting state
        new_u = projected_u;
        new_v = projected_u;
    } else {
        if (projected_u != projected_v) {
            new_u = projected_u;
            new_v = projected_v;
        } else {
            // Path between partners (e.g., (x, x+n, R))
            new_u = projected_u;
            new_v = projected_u + n;
        }
    }
    return State(new_u, new_v, result_set);
}

// ===== DAG GENERATION =====
inline DAGData build_dag(int order_num) {
    const int N_VERTICES = 2 * order_num;
    const Real N_f = 3.0;
    std::vector<std::pair<int, int>> INTERACTIONS;
    INTERACTIONS.reserve(order_num);
    for (int i = 1; i <= order_num; i++) {
        INTERACTIONS.push_back({i, i + order_num});
    }

    std::vector<std::unordered_map<State, Real, StateHash>> V_layers(N_VERTICES);
    std::vector<std::unordered_map<State, int, StateHash>> state_to_node_id(N_VERTICES + 1);
    int next_node_id = 0;

    State start_state(1, 1, {1});
    State normalized_start = shift_state(start_state, order_num);
    V_layers[0][normalized_start] = 1.0;
    state_to_node_id[0][normalized_start] = next_node_id++;

    std::vector<std::vector<RawEdge>> layer_edges(N_VERTICES);

    for (int i = 0; i < N_VERTICES - 1; i++) {
        for (const auto& kv : V_layers[i]) {
            const State& current_state = kv.first;
            const int src_node_id = state_to_node_id[i][current_state];
            const int u = current_state.u;
            const int v = current_state.v;
            const auto& R = current_state.R;

            std::unordered_set<int> unvisited;
            unvisited.reserve(N_VERTICES - static_cast<int>(R.size()));
            for (int w = 1; w <= N_VERTICES; w++) {
                if (R.find(w) == R.end()) unvisited.insert(w);
            }

            // Expand to all unvisited vertices
            for (int w : unvisited) {
                std::unordered_set<int> new_R = R;
                new_R.insert(w);
                State new_state(u, w, new_R);
                State normalized_new_state = shift_state(new_state, order_num);
                V_layers[i + 1][normalized_new_state] = 1.0;
                if (state_to_node_id[i + 1].find(normalized_new_state) == state_to_node_id[i + 1].end()) {
                    state_to_node_id[i + 1][normalized_new_state] = next_node_id++;
                }
                int dst_node_id = state_to_node_id[i + 1][normalized_new_state];
                layer_edges[i].push_back({src_node_id, dst_node_id, v - 1, w - 1, 1.0});
            }

            // Partner interaction edges
            std::vector<int> partner_vertices;
            for (const auto& pq : INTERACTIONS) {
                const int p = pq.first;
                const int q = pq.second;
                const bool p_in_R = (R.find(p) != R.end());
                const bool q_in_R = (R.find(q) != R.end());
                if (p_in_R != q_in_R) {
                    if (!p_in_R) partner_vertices.push_back(p);
                    if (!q_in_R) partner_vertices.push_back(q);
                }
            }
            std::sort(partner_vertices.begin(), partner_vertices.end());
            if (!partner_vertices.empty()) {
                int new_head = partner_vertices[0];
                std::unordered_set<int> new_R = R;
                new_R.insert(new_head);
                State new_state(new_head, new_head, new_R);
                State normalized_new_state = shift_state(new_state, order_num);
                V_layers[i + 1][normalized_new_state] = 1.0;
                if (state_to_node_id[i + 1].find(normalized_new_state) == state_to_node_id[i + 1].end()) {
                    state_to_node_id[i + 1][normalized_new_state] = next_node_id++;
                }
                int dst_node_id = state_to_node_id[i + 1][normalized_new_state];
                layer_edges[i].push_back({src_node_id, dst_node_id, v - 1, u - 1, -N_f});
            }
        }
    }

    // Sink connections from last layer
    const int sink_node_id = next_node_id++;
    for (const auto& kv : V_layers[N_VERTICES - 1]) {
        const State& state = kv.first;
        int src_node_id = state_to_node_id[N_VERTICES - 1][state];
        layer_edges[N_VERTICES - 1].push_back({src_node_id, sink_node_id, state.v - 1, state.u - 1, -N_f});
    }

    // ===== Sort edges within each layer for cache-friendly forward pass =====
    for (auto& layer : layer_edges) {
        std::sort(layer.begin(), layer.end(), [](const RawEdge& a, const RawEdge& b) {
            if (a.src != b.src) return a.src < b.src;        // group by src to reuse node_values[src]
            if (a.row != b.row) return a.row < b.row;        // then by row for M's row-major access
            if (a.col != b.col) return a.col < b.col;        // then by col
            return a.dst < b.dst;
        });
    }

    const int max_node_id = sink_node_id;
    return DAGData{std::move(layer_edges), sink_node_id, max_node_id};
}

// ===== FORWARD PROPAGATION (pointer-friendly) =====
thread_local std::vector<Real> tl_node_values;

inline Real compute_sink_value_fast(const DAGData& dag, const Real* RESTRICT M_matrix, int matrix_size) {
    auto& node = tl_node_values;
    if (node.size() < static_cast<size_t>(dag.max_node_id + 1)) node.resize(dag.max_node_id + 1);
    std::fill(node.begin(), node.begin() + dag.max_node_id + 1, 0.0);
    node[0] = 1.0;

    // Iterate layer by layer; edges are sorted primarily by src
    for (const auto& layer : dag.layer_edges) {
        const RawEdge* e = layer.data();
        const RawEdge* e_end = e + layer.size();
        int cur_src = -1;
        Real src_val = 0.0;
        for (; e != e_end; ++e) {
            if (e->src != cur_src) {
                cur_src = e->src;
                src_val = node[cur_src];
            }
            const Real mval = M_matrix[e->row * matrix_size + e->col];
            node[e->dst] += src_val * e->factor * mval;
        }
    }
    return node[dag.sink_node_id];
}

#endif // DAG_GENERATOR_H
