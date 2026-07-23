#ifndef GSPICE_TRANSIENT_STATE_STORE_HPP
#define GSPICE_TRANSIENT_STATE_STORE_HPP

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <vector>

namespace gspice {

struct TransientStateBlock {
    std::size_t offset = 0;
    std::size_t size = 0;
};

class TransientStateLayout {
public:
    TransientStateBlock allocate(std::size_t count) {
        if (sealed_) {
            throw std::logic_error("transient state layout is already sealed");
        }
        const TransientStateBlock block{totalSize_, count};
        totalSize_ += count;
        return block;
    }

    std::size_t totalSize() const { return totalSize_; }
    bool sealed() const { return sealed_; }
    void seal() { sealed_ = true; }

private:
    std::size_t totalSize_ = 0;
    bool sealed_ = false;
};

class ConstTransientStateView {
public:
    ConstTransientStateView() = default;
    ConstTransientStateView(const double* data, std::size_t size) : data_(data), size_(size) {}

    const double& operator[](std::size_t index) const {
        if (index >= size_) throw std::out_of_range("transient state index");
        return data_[index];
    }

    const double* data() const { return data_; }
    std::size_t size() const { return size_; }
    bool empty() const { return size_ == 0; }

private:
    const double* data_ = nullptr;
    std::size_t size_ = 0;
};

class TransientStateView {
public:
    TransientStateView() = default;
    TransientStateView(double* data, std::size_t size) : data_(data), size_(size) {}

    double& operator[](std::size_t index) {
        if (index >= size_) throw std::out_of_range("transient state index");
        return data_[index];
    }

    const double& operator[](std::size_t index) const {
        if (index >= size_) throw std::out_of_range("transient state index");
        return data_[index];
    }

    double* data() { return data_; }
    const double* data() const { return data_; }
    std::size_t size() const { return size_; }
    bool empty() const { return size_ == 0; }

    operator ConstTransientStateView() const { return {data_, size_}; }

private:
    double* data_ = nullptr;
    std::size_t size_ = 0;
};

// One contiguous store holds every device's transient state.  Candidate
// frames are separate from accepted history, so a rejected step is rolled
// back by restoring a small cursor rather than asking every device to clone
// arbitrary heap objects.
class TransientStateStore {
public:
    struct Checkpoint {
        std::size_t head = 0;
        std::size_t acceptedCount = 0;
        std::uint64_t serial = 0;
    };

    TransientStateStore(
        const TransientStateLayout& layout,
        std::size_t historyDepth,
        std::size_t speculativeDepth = 1)
        : width_(layout.totalSize()),
          historyDepth_(historyDepth),
          speculativeDepth_(speculativeDepth),
          frames_(historyDepth + speculativeDepth, std::vector<double>(width_, 0.0)) {
        if (!layout.sealed()) {
            throw std::invalid_argument("transient state layout must be sealed before storage is created");
        }
        if (historyDepth_ == 0 || speculativeDepth_ == 0) {
            throw std::invalid_argument("transient state depths must be positive");
        }
    }

    std::size_t width() const { return width_; }
    std::size_t historyDepth() const { return historyDepth_; }
    std::size_t acceptedCount() const { return acceptedCount_; }

    Checkpoint checkpoint() const { return {head_, acceptedCount_, serial_}; }

    void rollback(const Checkpoint& checkpoint) {
        if (checkpoint.serial > serial_ || serial_ - checkpoint.serial > speculativeDepth_) {
            throw std::logic_error("transient state checkpoint is outside the rollback window");
        }
        head_ = checkpoint.head;
        acceptedCount_ = checkpoint.acceptedCount;
        serial_ = checkpoint.serial;
        candidatePrepared_ = false;
    }

    void prepareCandidate() {
        frames_[candidateIndex()] = frames_[head_];
        candidatePrepared_ = true;
    }

    void acceptCandidate() {
        if (!candidatePrepared_) {
            throw std::logic_error("transient candidate was not prepared");
        }
        head_ = candidateIndex();
        acceptedCount_ = std::min(acceptedCount_ + 1, historyDepth_);
        ++serial_;
        candidatePrepared_ = false;
    }

    TransientStateView initial(const TransientStateBlock& block) {
        validateBlock(block);
        return {frames_[head_].data() + block.offset, block.size};
    }

    ConstTransientStateView accepted(const TransientStateBlock& block, std::size_t age = 0) const {
        validateBlock(block);
        if (age >= historyDepth_ || age > acceptedCount_) {
            throw std::out_of_range("transient state history age is unavailable");
        }
        const auto& frame = frames_[historyIndex(age)];
        return {frame.data() + block.offset, block.size};
    }

    TransientStateView candidate(const TransientStateBlock& block) {
        validateBlock(block);
        if (!candidatePrepared_) {
            throw std::logic_error("transient candidate was not prepared");
        }
        auto& frame = frames_[candidateIndex()];
        return {frame.data() + block.offset, block.size};
    }

private:
    void validateBlock(const TransientStateBlock& block) const {
        if (block.offset > width_ || block.size > width_ - block.offset) {
            throw std::out_of_range("transient state block is outside the registered layout");
        }
    }

    std::size_t candidateIndex() const { return (head_ + 1) % frames_.size(); }

    std::size_t historyIndex(std::size_t age) const {
        return (head_ + frames_.size() - (age % frames_.size())) % frames_.size();
    }

    std::size_t width_ = 0;
    std::size_t historyDepth_ = 0;
    std::size_t speculativeDepth_ = 0;
    std::vector<std::vector<double>> frames_;
    std::size_t head_ = 0;
    std::size_t acceptedCount_ = 0;
    std::uint64_t serial_ = 0;
    bool candidatePrepared_ = false;
};

struct OpaqueTransientStateBlock {
    std::size_t offset = 0;
    std::size_t size = 0;
};

class OpaqueTransientStateLayout {
public:
    OpaqueTransientStateBlock allocate(std::size_t bytes) {
        if (sealed_) throw std::logic_error("opaque transient state layout is already sealed");
        const OpaqueTransientStateBlock block{totalSize_, bytes};
        totalSize_ += bytes;
        return block;
    }
    std::size_t totalSize() const { return totalSize_; }
    bool sealed() const { return sealed_; }
    void seal() { sealed_ = true; }

private:
    std::size_t totalSize_ = 0;
    bool sealed_ = false;
};

// Transactional byte storage is used for compact-model memory that cannot be
// represented as doubles. It has the same bounded rollback semantics as the
// numeric state store and performs one contiguous frame copy per candidate.
class OpaqueTransientStateStore {
public:
    struct Checkpoint {
        std::size_t head = 0;
        std::size_t acceptedCount = 0;
        std::uint64_t serial = 0;
    };

    OpaqueTransientStateStore(
        const OpaqueTransientStateLayout& layout,
        std::size_t historyDepth,
        std::size_t speculativeDepth = 1)
        : width_(layout.totalSize()), historyDepth_(historyDepth),
          speculativeDepth_(speculativeDepth),
          frames_(historyDepth + speculativeDepth, std::vector<std::byte>(width_)) {
        if (!layout.sealed()) throw std::invalid_argument("opaque state layout must be sealed");
        if (historyDepth == 0 || speculativeDepth == 0) {
            throw std::invalid_argument("opaque transient state depths must be positive");
        }
    }

    Checkpoint checkpoint() const { return {head_, acceptedCount_, serial_}; }

    void rollback(const Checkpoint& checkpoint) {
        if (checkpoint.serial > serial_ || serial_ - checkpoint.serial > speculativeDepth_) {
            throw std::logic_error("opaque transient checkpoint is outside the rollback window");
        }
        head_ = checkpoint.head;
        acceptedCount_ = checkpoint.acceptedCount;
        serial_ = checkpoint.serial;
        candidatePrepared_ = false;
    }

    void prepareCandidate() {
        frames_[candidateIndex()] = frames_[head_];
        candidatePrepared_ = true;
    }

    void acceptCandidate() {
        if (!candidatePrepared_) throw std::logic_error("opaque transient candidate was not prepared");
        head_ = candidateIndex();
        acceptedCount_ = std::min(acceptedCount_ + 1, historyDepth_);
        ++serial_;
        candidatePrepared_ = false;
    }

    std::byte* initial(const OpaqueTransientStateBlock& block) {
        validate(block);
        return frames_[head_].data() + block.offset;
    }

    const std::byte* accepted(const OpaqueTransientStateBlock& block, std::size_t age = 0) const {
        validate(block);
        if (age >= historyDepth_ || age > acceptedCount_) {
            throw std::out_of_range("opaque transient state history age is unavailable");
        }
        return frames_[historyIndex(age)].data() + block.offset;
    }

    std::byte* candidate(const OpaqueTransientStateBlock& block) {
        validate(block);
        if (!candidatePrepared_) throw std::logic_error("opaque transient candidate was not prepared");
        return frames_[candidateIndex()].data() + block.offset;
    }

private:
    void validate(const OpaqueTransientStateBlock& block) const {
        if (block.offset > width_ || block.size > width_ - block.offset) {
            throw std::out_of_range("opaque transient block is outside the registered layout");
        }
    }
    std::size_t candidateIndex() const { return (head_ + 1) % frames_.size(); }
    std::size_t historyIndex(std::size_t age) const {
        return (head_ + frames_.size() - age % frames_.size()) % frames_.size();
    }

    std::size_t width_ = 0;
    std::size_t historyDepth_ = 0;
    std::size_t speculativeDepth_ = 0;
    std::vector<std::vector<std::byte>> frames_;
    std::size_t head_ = 0;
    std::size_t acceptedCount_ = 0;
    std::uint64_t serial_ = 0;
    bool candidatePrepared_ = false;
};

} // namespace gspice

#endif // GSPICE_TRANSIENT_STATE_STORE_HPP
