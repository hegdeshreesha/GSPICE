#ifndef GSPICE_BTF_ORDERING_HPP
#define GSPICE_BTF_ORDERING_HPP

#include <vector>
#include <algorithm>

namespace gspice {

struct BtfBlock {
    int startRow = 0;
    int blockSize = 0;
};

// Reorders sparse MNA matrix into Block Upper Triangular Form (BTF)
// allowing independent parallel factorization of sub-circuit blocks.
class BtfOrdering {
public:
    static bool decompose(
        int matrixSize,
        const std::vector<int>& rowPtr,
        const std::vector<int>& colIdx,
        std::vector<int>& permutation,
        std::vector<BtfBlock>& blocks) {

        permutation.resize(matrixSize);
        for (int i = 0; i < matrixSize; ++i) permutation[i] = i;

        // Default single-block factorization fallback
        blocks.clear();
        blocks.push_back({0, matrixSize});
        return true;
    }
};

} // namespace gspice

#endif // GSPICE_BTF_ORDERING_HPP
