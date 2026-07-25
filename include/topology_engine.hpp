#ifndef GSPICE_TOPOLOGY_ENGINE_HPP
#define GSPICE_TOPOLOGY_ENGINE_HPP

#include <vector>
#include <numeric>

namespace gspice {

class TopologyEngine {
public:
    explicit TopologyEngine(int numNodes = 0) {
        reset(numNodes);
    }

    void reset(int numNodes) {
        parent_.resize(numNodes);
        std::iota(parent_.begin(), parent_.end(), 0);
        rank_.assign(numNodes, 0);
    }

    int find(int node) {
        if (node < 0 || node >= static_cast<int>(parent_.size())) return node;
        if (parent_[node] != node) {
            parent_[node] = find(parent_[node]); // Path compression
        }
        return parent_[node];
    }

    bool unite(int nodeA, int nodeB) {
        int rootA = find(nodeA);
        int rootB = find(nodeB);
        if (rootA == rootB) return false;

        // Ground node (0) always stays root
        if (rootA == 0) {
            parent_[rootB] = rootA;
        } else if (rootB == 0) {
            parent_[rootA] = rootB;
        } else if (rank_[rootA] < rank_[rootB]) {
            parent_[rootA] = rootB;
        } else if (rank_[rootA] > rank_[rootB]) {
            parent_[rootB] = rootA;
        } else {
            parent_[rootB] = rootA;
            rank_[rootA]++;
        }
        return true;
    }

    bool isCollapsed(int nodeA, int nodeB) {
        return find(nodeA) == find(nodeB);
    }

private:
    std::vector<int> parent_;
    std::vector<int> rank_;
};

} // namespace gspice

#endif // GSPICE_TOPOLOGY_ENGINE_HPP
