#ifndef GSPICE_CIRCUIT_TOPOLOGY_HPP
#define GSPICE_CIRCUIT_TOPOLOGY_HPP

#include "matrix.hpp"
#include "sparse_matrix.hpp"
#include <algorithm>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace gspice {

struct NodeCollapse {
    int first = -1;
    int second = -1; // -1 denotes ground.
};

// Fixed-dimension circuit-wide node equivalence. Device stamps remain in their
// original coordinates. At an elaboration boundary the assembled equations are
// projected with P^T A P and P^T b; unused alias rows enforce x_alias=x_root.
// This keeps branch indices and transient vectors stable.
class CircuitTopology {
public:
    CircuitTopology() = default;
    explicit CircuitTopology(int nodeCount) { reset(nodeCount); }

    void reset(int nodeCount) {
        if (nodeCount < 0) throw std::invalid_argument("negative circuit node count");
        node_count_ = nodeCount;
        canonical_.resize(static_cast<std::size_t>(nodeCount));
        std::iota(canonical_.begin(), canonical_.end(), 0);
        ++revision_;
    }

    bool rebuild(int nodeCount, const std::vector<NodeCollapse>& collapses) {
        if (nodeCount < 0) throw std::invalid_argument("negative circuit node count");
        const int ground = nodeCount;
        std::vector<int> parent(static_cast<std::size_t>(nodeCount + 1));
        std::iota(parent.begin(), parent.end(), 0);
        const auto find = [&](int value, auto&& self) -> int {
            if (parent[static_cast<std::size_t>(value)] != value) {
                parent[static_cast<std::size_t>(value)] =
                    self(parent[static_cast<std::size_t>(value)], self);
            }
            return parent[static_cast<std::size_t>(value)];
        };
        const auto unite = [&](int lhs, int rhs) {
            int a = find(lhs, find);
            int b = find(rhs, find);
            if (a == b) return;
            if (a == ground || b == ground) {
                parent[static_cast<std::size_t>(a == ground ? b : a)] = ground;
            } else if (a < b) {
                parent[static_cast<std::size_t>(b)] = a;
            } else {
                parent[static_cast<std::size_t>(a)] = b;
            }
        };

        for (const auto& collapse : collapses) {
            if (collapse.first < 0 || collapse.first >= nodeCount) {
                throw std::runtime_error("node collapse references an invalid first circuit node");
            }
            const int second = collapse.second < 0 ? ground : collapse.second;
            if (second < 0 || second > nodeCount) {
                throw std::runtime_error("node collapse references an invalid second circuit node");
            }
            unite(collapse.first, second);
        }

        std::vector<int> updated(static_cast<std::size_t>(nodeCount));
        for (int node = 0; node < nodeCount; ++node) {
            const int root = find(node, find);
            updated[static_cast<std::size_t>(node)] = root == ground ? -1 : root;
        }
        const bool changed = node_count_ != nodeCount || updated != canonical_;
        node_count_ = nodeCount;
        canonical_ = std::move(updated);
        if (changed) ++revision_;
        return changed;
    }

    int canonicalNode(int node) const {
        if (node < 0) return -1;
        if (node >= node_count_) return node;
        return canonical_[static_cast<std::size_t>(node)];
    }

    bool hasAliases() const {
        for (int node = 0; node < node_count_; ++node) {
            if (canonical_[static_cast<std::size_t>(node)] != node) return true;
        }
        return false;
    }

    int aliasCount() const {
        int count = 0;
        for (int node = 0; node < node_count_; ++node) {
            if (canonical_[static_cast<std::size_t>(node)] != node) ++count;
        }
        return count;
    }

    std::uint64_t revision() const { return revision_; }

    template <typename T>
    void projectSolution(Vector<T>& values) const {
        for (int node = 0; node < node_count_ && node < values.getSize(); ++node) {
            const int root = canonical_[static_cast<std::size_t>(node)];
            if (root == node) continue;
            values[node] = root < 0 ? T(0) : values[root];
        }
    }

    template <typename T>
    void projectRhs(Vector<T>& rhs) const {
        const auto original = rhs.snapshotData();
        rhs.clear();
        const int count = std::min(node_count_, static_cast<int>(original.size()));
        for (int row = 0; row < static_cast<int>(original.size()); ++row) {
            const int mapped = row < count
                ? canonical_[static_cast<std::size_t>(row)]
                : row;
            if (mapped >= 0) rhs.add(mapped, original[static_cast<std::size_t>(row)]);
        }
    }

    template <typename T>
    void apply(SparseMatrix<T>& matrix, Vector<T>& rhs) const {
        if (!hasAliases()) return;
        const auto entries = matrix.getEntries();
        const auto originalRhs = rhs.snapshotData();
        matrix.clear();
        rhs.clear();

        for (const auto& entry : entries) {
            const int row = canonicalUnknown(entry.row);
            const int col = canonicalUnknown(entry.col);
            if (row >= 0 && col >= 0) matrix.add(row, col, entry.value);
        }
        for (int row = 0; row < static_cast<int>(originalRhs.size()); ++row) {
            const int mapped = canonicalUnknown(row);
            if (mapped >= 0) rhs.add(mapped, originalRhs[static_cast<std::size_t>(row)]);
        }

        for (int alias = 0; alias < node_count_ && alias < matrix.getSize(); ++alias) {
            const int root = canonical_[static_cast<std::size_t>(alias)];
            if (root == alias) continue;
            matrix.add(alias, alias, T(1));
            if (root >= 0) matrix.add(alias, root, T(-1));
        }
    }

    // Harmonic and periodic-small-signal systems use
    // flat_index=base_unknown*slots+harmonic. Apply the same circuit-node
    // equivalence independently to every harmonic slot.
    template <typename T>
    void applyInterleaved(
        SparseMatrix<T>& matrix,
        Vector<T>& rhs,
        int slots) const {
        if (!hasAliases()) return;
        if (slots <= 0) throw std::invalid_argument("topology slot count must be positive");
        const auto entries = matrix.getEntries();
        const auto originalRhs = rhs.snapshotData();
        matrix.clear();
        rhs.clear();
        const auto mapIndex = [&](int index) {
            if (index < 0) return -1;
            const int base = index / slots;
            const int slot = index % slots;
            if (base >= node_count_) return index;
            const int root = canonical_[static_cast<std::size_t>(base)];
            return root < 0 ? -1 : root * slots + slot;
        };
        for (const auto& entry : entries) {
            const int row = mapIndex(entry.row);
            const int col = mapIndex(entry.col);
            if (row >= 0 && col >= 0) matrix.add(row, col, entry.value);
        }
        for (int row = 0; row < static_cast<int>(originalRhs.size()); ++row) {
            const int mapped = mapIndex(row);
            if (mapped >= 0) rhs.add(mapped, originalRhs[static_cast<std::size_t>(row)]);
        }
        for (int alias = 0; alias < node_count_; ++alias) {
            const int root = canonical_[static_cast<std::size_t>(alias)];
            if (root == alias) continue;
            for (int slot = 0; slot < slots; ++slot) {
                const int row = alias * slots + slot;
                if (row >= matrix.getSize()) continue;
                matrix.add(row, row, T(1));
                if (root >= 0) matrix.add(row, root * slots + slot, T(-1));
            }
        }
    }

private:
    int canonicalUnknown(int unknown) const {
        return unknown >= 0 && unknown < node_count_
            ? canonical_[static_cast<std::size_t>(unknown)]
            : unknown;
    }

    int node_count_ = 0;
    std::vector<int> canonical_;
    std::uint64_t revision_ = 0;
};

} // namespace gspice

#endif
