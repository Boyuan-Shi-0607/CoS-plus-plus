#ifndef TYPED_DAG_GENERATOR_H
#define TYPED_DAG_GENERATOR_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cstdint>

// If you already have a common Real typedef, replace this with it.
#ifndef TYPED_DAG_REAL_T
#define TYPED_DAG_REAL_T double
#endif

#if defined(__GNUC__) || defined(__clang__)
  #define RESTRICT __restrict__
#else
  #define RESTRICT
#endif

namespace typed_dag {

using Real = TYPED_DAG_REAL_T;

// ==================== Edge typing ====================
enum class EdgeKind : uint8_t { Mat = 0, Shift = 1 };

struct RawEdge {
    int src;        // node id
    int dst;        // node id
    int row;        // (v-1)%n
    int col;        // (w-1)%n or (u-1)%n when closing
    Real factor;    // +1 or -1 multiplier
    int label;      // flavour index in [0, N_f)
    EdgeKind kind;  // Mat -> A_matrices[label][row*n + col], Shift -> shift_close[label]
};

struct DAGData {
    std::vector<std::vector<RawEdge>> layer_edges; // edges grouped by layer (forward-only)
    int sink_node_id = -1;
    int max_node_id  = -1; // inclusive
};

// ==================== State & hashing ====================
// The state mirrors the Python DP state key: (head_info, v, R, unpaired_info)
// head_info = (u, u_type, u_label)
// R: set of visited vertices
// unpaired: map of vertex in R whose partner is not in R -> (type,label) at time of visit
struct State {
    int u; // head anchor vertex
    int v; // current frontier vertex
    std::unordered_set<int> R; // visited set
    std::unordered_map<int, std::pair<uint8_t,int>> unpaired; // vertex -> (type,label)
    uint8_t u_type; // 0 = normal, 1 = shift
    int u_label;    // current head label [0, N_f)

    State() : u(0), v(0), u_type(0), u_label(0) {}
    State(int u_, int v_, std::unordered_set<int> R_,
          std::unordered_map<int, std::pair<uint8_t,int>> unp_,
          uint8_t t_, int lab_)
        : u(u_), v(v_), R(std::move(R_)), unpaired(std::move(unp_)), u_type(t_), u_label(lab_) {}

    bool operator==(const State& other) const {
        return u == other.u && v == other.v && u_type == other.u_type && u_label == other.u_label
            && R == other.R && unpaired == other.unpaired;
    }
};

struct StateHash {
    size_t operator()(State const& s) const noexcept {
        size_t h = 0;
        auto hi = std::hash<int>{};
        auto hb = std::hash<uint8_t>{};
        auto mix = [&](size_t x){ h ^= x + 0x9e3779b97f4a7c15ULL + (h<<6) + (h>>2); };
        mix(hi(s.u));
        mix(hi(s.v));
        mix(hb(s.u_type));
        mix(hi(s.u_label));
        // hash R deterministically
        std::vector<int> elems(s.R.begin(), s.R.end());
        std::sort(elems.begin(), elems.end());
        for (int r : elems) mix(hi(r));
        // hash unpaired deterministically
        struct Item { int key; uint8_t t; int lab; };
        std::vector<Item> items; items.reserve(s.unpaired.size());
        for (auto &kv : s.unpaired) items.push_back({kv.first, kv.second.first, kv.second.second});
        std::sort(items.begin(), items.end(), [](const Item&a,const Item&b){ return a.key<b.key || (a.key==b.key && (a.t<b.t || (a.t==b.t && a.lab<b.lab))); });
        for (auto &it : items) { mix(hi(it.key)); mix(hb(it.t)); mix(hi(it.lab)); }
        return h;
    }
};

// ==================== Canonicalization ====================
// Equivalent to Python get_canonical_map() + applying to (u,v,R,unpaired).
// Map vertices modulo n and split collisions by +n in the order induced by sorted R.
struct Canonicalized {
    State st;
};

inline Canonicalized canonicalize(const State& s, int n) {
    // Build mapping only over vertices in R
    std::vector<int> sorted_R; sorted_R.reserve(s.R.size());
    for (int x : s.R) sorted_R.push_back(x);
    std::sort(sorted_R.begin(), sorted_R.end());

    std::unordered_map<int,int> vmap; vmap.reserve(sorted_R.size()*2);
    std::unordered_map<int,int> counts; counts.reserve(sorted_R.size()*2);

    for (int old_v : sorted_R) {
        int base = (old_v - 1) % n + 1;
        int& c = counts[base];
        int mapped = base + (c>0 ? n : 0);
        ++c;
        vmap[old_v] = mapped;
    }

    // Map R
    std::unordered_set<int> canonR; canonR.reserve(s.R.size()*2);
    for (int old_v : sorted_R) canonR.insert(vmap[old_v]);

    // Map unpaired keys
    std::unordered_map<int, std::pair<uint8_t,int>> canonUnp; canonUnp.reserve(s.unpaired.size()*2);
    for (auto const& kv : s.unpaired) {
        int old_k = kv.first; auto val = kv.second;
        auto it = vmap.find(old_k);
        if (it != vmap.end()) canonUnp[it->second] = val;
        // else: by construction unpaired keys are subset of R, so should exist
    }

    // Map u,v (both must be in R)
    int nu = vmap.at(s.u);
    int nv = vmap.at(s.v);

    State out(nu, nv, std::move(canonR), std::move(canonUnp), s.u_type, s.u_label);
    return { std::move(out) };
}

// ==================== DAG builder (typed, exact Python parity) ====================
inline DAGData build_dag_typed(int n, int N_f) {
    const int N_VERTICES = 2*n;

    auto partner = [n](int x){ return (x<=n) ? x+n : x-n; };

    // Layers of unique states
    std::vector<std::unordered_map<State, int, StateHash>> state_id(N_VERTICES);
    std::vector<std::vector<RawEdge>> layer_edges(N_VERTICES); // last index will be used for sink edges

    int next_id = 0;

    // ========== Initialize layer 0: for each label, both normal and shift heads at vertex 1 ==========
    for (int lab=0; lab<N_f; ++lab) {
        for (uint8_t t=0; t<2; ++t) {
            std::unordered_set<int> R0 = {1};
            std::unordered_map<int, std::pair<uint8_t,int>> unp0; unp0[1] = {t, lab};
            State s0(1,1,std::move(R0), std::move(unp0), t, lab);
            auto cs = canonicalize(s0, n).st;
            if (!state_id[0].count(cs)) state_id[0][cs] = next_id++;
        }
    }

    // ========== Forward layers (0 .. N_VERTICES-2) ==========
    for (int layer=0; layer<N_VERTICES-1; ++layer) {
        for (auto const& kv : state_id[layer]) {
            const State& st = kv.first;
            int src_node = kv.second;

            // ---- Rule 1: Continue current cycle (only if normal head) ----
            if (st.u_type == 0) {
                // unvisited = [1..2n] \ R
                for (int w=1; w<=N_VERTICES; ++w) if (!st.R.count(w)) {
                    auto newR = st.R; newR.insert(w);
                    auto newUnp = st.unpaired; // copy
                    int pw = partner(w);
                    if (newUnp.count(pw)) newUnp.erase(pw); else newUnp[w] = {st.u_type, st.u_label};

                    State nxt(st.u, w, std::move(newR), std::move(newUnp), st.u_type, st.u_label);
                    auto cn = canonicalize(nxt, n).st;
                    int dst_id;
                    auto it = state_id[layer+1].find(cn);
                    if (it == state_id[layer+1].end()) dst_id = state_id[layer+1][cn] = next_id++;
                    else dst_id = it->second;

                    // Edge pulls A_matrices[st.u_label][(v-1)%n, (w-1)%n]
                    RawEdge e{src_node, dst_id, (st.v-1)%n, (w-1)%n, +1.0, st.u_label, EdgeKind::Mat};
                    layer_edges[layer].push_back(e);
                }
            }

            // ---- Rule 2: Close the current cycle and start a new one ----
            // closed_value: -A[v,u] if normal head, else -shift_close[label]
            // partner_vertices = [partner[p] for p in unpaired.keys()]
            std::vector<int> partner_vertices; partner_vertices.reserve(st.unpaired.size());
            for (auto const& up : st.unpaired) partner_vertices.push_back(partner(up.first));
            if (!partner_vertices.empty()) {
                std::sort(partner_vertices.begin(), partner_vertices.end());
                int new_head_vertex = partner_vertices[0];
                int partner_of_new_head = partner(new_head_vertex);
                auto itup = st.unpaired.find(partner_of_new_head);
                // Should exist by construction
                uint8_t partner_type = itup->second.first;
                int     partner_label = itup->second.second;

                auto newR = st.R; newR.insert(new_head_vertex);
                auto newUnp = st.unpaired; newUnp.erase(partner_of_new_head);

                auto emit_close_edge = [&](int dst_id){
                    RawEdge e; e.src = src_node; e.dst = dst_id; e.factor = -1.0; e.label = st.u_label;
                    if (st.u_type == 0) { e.kind = EdgeKind::Mat; e.row = (st.v-1)%n; e.col = (st.u-1)%n; }
                    else { e.kind = EdgeKind::Shift; e.row = 0; e.col = 0; }
                    layer_edges[layer].push_back(e);
                };

                // (A) Open new NORMAL head
                if (partner_type == 1) {
                    // Must inherit label
                    int new_lab = partner_label;
                    State nxt(new_head_vertex, new_head_vertex, newR, newUnp, /*type*/0, /*label*/new_lab);
                    auto cn = canonicalize(nxt, n).st;
                    int dst_id; auto it = state_id[layer+1].find(cn);
                    if (it == state_id[layer+1].end()) dst_id = state_id[layer+1][cn] = next_id++;
                    else dst_id = it->second;
                    emit_close_edge(dst_id);
                } else {
                    // Can open with any label
                    for (int new_lab=0; new_lab<N_f; ++new_lab) {
                        State nxt(new_head_vertex, new_head_vertex, newR, newUnp, /*type*/0, /*label*/new_lab);
                        auto cn = canonicalize(nxt, n).st;
                        int dst_id; auto it = state_id[layer+1].find(cn);
                        if (it == state_id[layer+1].end()) dst_id = state_id[layer+1][cn] = next_id++;
                        else dst_id = it->second;
                        emit_close_edge(dst_id);
                    }
                }

                // (B) Open new SHIFT head (inherits partner_label)
                {
                    int new_lab = partner_label;
                    State nxt(new_head_vertex, new_head_vertex, newR, newUnp, /*type*/1, /*label*/new_lab);
                    auto cn = canonicalize(nxt, n).st;
                    int dst_id; auto it = state_id[layer+1].find(cn);
                    if (it == state_id[layer+1].end()) dst_id = state_id[layer+1][cn] = next_id++;
                    else dst_id = it->second;
                    emit_close_edge(dst_id);
                }
            }
        }
    }

    // ========== Create sink and connect last layer states with a final closing edge ==========
    int sink_id = next_id++;
    // The final contribution for each state is exactly the "close" value with no new state.
    for (auto const& kv : state_id[N_VERTICES-1]) {
        const State& st = kv.first; int src = kv.second;
        RawEdge e; e.src = src; e.dst = sink_id; e.factor = -1.0; e.label = st.u_label;
        if (st.u_type == 0) { e.kind = EdgeKind::Mat; e.row = (st.v-1)%n; e.col = (st.u-1)%n; }
        else { e.kind = EdgeKind::Shift; e.row = 0; e.col = 0; }
        layer_edges[N_VERTICES-1].push_back(e);
    }

    // sort edges in each layer for cache-friendly forward pass
    for (auto& L : layer_edges) {
        std::sort(L.begin(), L.end(), [](const RawEdge& a, const RawEdge& b){
            if (a.src != b.src) return a.src < b.src;
            if ((int)a.kind != (int)b.kind) return (int)a.kind < (int)b.kind;
            if (a.label != b.label) return a.label < b.label;
            if (a.row != b.row) return a.row < b.row;
            if (a.col != b.col) return a.col < b.col;
            return a.dst < b.dst;
        });
    }

    DAGData dag; dag.layer_edges = std::move(layer_edges); dag.sink_node_id = sink_id; dag.max_node_id = sink_id; return dag;
}

// ==================== Fast sink evaluation ====================
// A_matrices: size N_f, each a flat row-major n*n matrix.
// shift_close: size N_f.
// n: matrix dimension per flavour.
thread_local std::vector<Real> tl_nodes;

inline Real compute_sink_value_fast(const DAGData& dag,
                                    const std::vector<std::vector<Real>>& A_matrices,
                                    int n,
                                    const std::vector<Real>& shift_close) {
    auto& node = tl_nodes;
    const int need = dag.max_node_id + 1;
    if ((int)node.size() < need) node.resize(need);
    std::fill(node.begin(), node.begin() + need, 0.0);

    // Initialize all layer-0 sources to 1.0 (Python DP has 2*N_f start states, each weight=1)
    if (!dag.layer_edges.empty()) {
        const auto& L0 = dag.layer_edges.front();
        for (const auto& e : L0) node[e.src] = 1.0;
    }

    // Walk layers in order
    for (const auto& L : dag.layer_edges) {
        int cur_src = -1; Real srcv = 0.0;
        for (const RawEdge& e : L) {
            if (e.src != cur_src) { cur_src = e.src; srcv = node[cur_src]; }
            Real mval;
            if (e.kind == EdgeKind::Mat) {
                const auto& M = A_matrices[e.label];
                mval = M[e.row * n + e.col];
            } else {
                mval = shift_close[e.label];
            }
            node[e.dst] += srcv * e.factor * mval;
        }
    }
    return node[dag.sink_node_id];
}

} // namespace typed_dag

#endif // TYPED_DAG_GENERATOR_H
