#ifndef GSPICE_TICER_HPP
#define GSPICE_TICER_HPP

// ---------------------------------------------------------------------------
// TICER — Time-Domain Interconnect Reduction (Feature A).
//
// Background:
//   Post-layout parasitic extraction generates netlists with tens or hundreds
//   of thousands of linear RC/RCCX nodes. Evaluating and solving these massive
//   linear networks at every timestep dominates memory and execution time.
//
// TICER Algorithm:
//   1. Identify linear RC sub-networks connected between active device ports.
//   2. For each internal node k, compute local time constant tau_k = R_k * C_k.
//   3. Compare tau_k against cutoff frequency f_max:
//        nodes with tau_k < 1 / (2 * pi * f_max) are eliminated.
//   4. Elimination uses Kron reduction / Star-Delta transformation:
//        G_ij' = G_ij - (G_ik * G_kj) / G_kk
//      Capacitance C_k is redistributed to adjacent nodes weighted by G_ik / G_kk.
//   5. Active device terminals and user .SAVE nodes are designated PORTS and
//      are strictly preserved.
// ---------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace gspice {

struct TicerResistor {
    int node1;
    int node2;
    double resistance;
};

struct TicerCapacitor {
    int node1;
    int node2;
    double capacitance;
};

struct TicerReductionResult {
    int original_nodes = 0;
    int reduced_nodes = 0;
    int eliminated_nodes = 0;
    int original_resistors = 0;
    int reduced_resistors = 0;
    int original_capacitors = 0;
    int reduced_capacitors = 0;
    double f_max = 1e9;
};

class TicerReducer {
public:
    explicit TicerReducer(double f_max = 1e9) : f_max_(f_max) {}

    /// Perform TICER reduction on an RC network.
    /// `port_nodes` contains node IDs that MUST be preserved (ports, device terminals, save nodes).
    TicerReductionResult reduce(
        int num_nodes,
        const std::vector<TicerResistor>& resistors,
        const std::vector<TicerCapacitor>& capacitors,
        const std::set<int>& port_nodes,
        std::vector<TicerResistor>& out_resistors,
        std::vector<TicerCapacitor>& out_capacitors) {

        TicerReductionResult res;
        res.original_nodes = num_nodes;
        res.original_resistors = static_cast<int>(resistors.size());
        res.original_capacitors = static_cast<int>(capacitors.size());
        res.f_max = f_max_;

        if (num_nodes <= 0 || resistors.empty()) {
            out_resistors = resistors;
            out_capacitors = capacitors;
            res.reduced_nodes = num_nodes;
            res.reduced_resistors = res.original_resistors;
            res.reduced_capacitors = res.original_capacitors;
            return res;
        }

        // Build Nodal Conductance and Capacitance maps
        std::map<std::pair<int, int>, double> G_map;
        std::unordered_map<int, double> C_node;

        for (const auto& r : resistors) {
            if (r.node1 == r.node2 || r.resistance <= 0.0) continue;
            const double g = 1.0 / r.resistance;
            const int u = std::min(r.node1, r.node2);
            const int v = std::max(r.node1, r.node2);
            G_map[{u, v}] += g;
        }

        for (const auto& c : capacitors) {
            if (c.node1 == c.node2 || c.capacitance <= 0.0) continue;
            if (c.node1 >= 0) C_node[c.node1] += c.capacitance;
            if (c.node2 >= 0) C_node[c.node2] += c.capacitance;
        }

        // Find candidate internal nodes for elimination
        std::set<int> internal_candidates;
        for (int i = 0; i < num_nodes; ++i) {
            if (port_nodes.find(i) == port_nodes.end()) {
                internal_candidates.insert(i);
            }
        }

        const double tau_limit = 1.0 / (2.0 * M_PI * f_max_);

        // Perform elimination pass
        std::unordered_set<int> eliminated;
        for (int k : internal_candidates) {
            // Find all neighbors of k
            std::vector<std::pair<int, double>> neighbors;
            double G_kk = 0.0;

            for (auto it = G_map.begin(); it != G_map.end(); ++it) {
                const int u = it->first.first;
                const int v = it->first.second;
                const double g = it->second;

                if (u == k && eliminated.find(v) == eliminated.end()) {
                    neighbors.push_back({v, g});
                    G_kk += g;
                } else if (v == k && eliminated.find(u) == eliminated.end()) {
                    neighbors.push_back({u, g});
                    G_kk += g;
                }
            }

            if (G_kk <= 1e-15) continue;

            const double R_eq = 1.0 / G_kk;
            const double C_k = C_node[k];
            const double tau_k = R_eq * C_k;

            // Eliminate if time constant is faster than f_max cutoff
            if (tau_k < tau_limit) {
                eliminated.insert(k);

                // Star-Delta Kron reduction: add G_ij' = (G_ik * G_jk) / G_kk between neighbors
                const std::size_t num_nbrs = neighbors.size();
                for (std::size_t i = 0; i < num_nbrs; ++i) {
                    for (std::size_t j = i + 1; j < num_nbrs; ++j) {
                        const int n1 = neighbors[i].first;
                        const int n2 = neighbors[j].first;
                        const double g1 = neighbors[i].second;
                        const double g2 = neighbors[j].second;
                        const double g_new = (g1 * g2) / G_kk;

                        const int u = std::min(n1, n2);
                        const int v = std::max(n1, n2);
                        G_map[{u, v}] += g_new;
                    }
                    // Redistribute capacitance C_k to neighbor i
                    const double weight = neighbors[i].second / G_kk;
                    C_node[neighbors[i].first] += C_k * weight;
                }
            }
        }

        // Reconstruct reduced resistor and capacitor lists
        out_resistors.clear();
        for (const auto& [pair, g] : G_map) {
            const int u = pair.first;
            const int v = pair.second;
            if (eliminated.find(u) != eliminated.end() || eliminated.find(v) != eliminated.end()) continue;
            if (g > 1e-15) {
                out_resistors.push_back({u, v, 1.0 / g});
            }
        }

        out_capacitors.clear();
        for (const auto& [node, c] : C_node) {
            if (eliminated.find(node) != eliminated.end()) continue;
            if (c > 1e-18) {
                out_capacitors.push_back({node, -1, c}); // -1 = ground
            }
        }

        res.eliminated_nodes = static_cast<int>(eliminated.size());
        res.reduced_nodes = num_nodes - res.eliminated_nodes;
        res.reduced_resistors = static_cast<int>(out_resistors.size());
        res.reduced_capacitors = static_cast<int>(out_capacitors.size());

        return res;
    }

private:
    double f_max_ = 1e9;
};

} // namespace gspice

#endif // GSPICE_TICER_HPP
