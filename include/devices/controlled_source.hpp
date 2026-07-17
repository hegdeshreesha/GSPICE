#ifndef GSPICE_CONTROLLED_SOURCE_HPP
#define GSPICE_CONTROLLED_SOURCE_HPP

#include "device.hpp"
#include <complex>
#include <string>
#include <utility>

namespace gspice {

class VoltageControlledCurrentSource : public Device {
public:
    VoltageControlledCurrentSource(
        const std::string& name,
        int nodePos,
        int nodeNeg,
        int ctrlPos,
        int ctrlNeg,
        double transconductance)
        : Device(name),
          nodePos_(nodePos),
          nodeNeg_(nodeNeg),
          ctrlPos_(ctrlPos),
          ctrlNeg_(ctrlNeg),
          transconductance_(transconductance) {}

    void dcStamp(
        SparseMatrixReal& J,
        VectorReal& b,
        const VectorReal& x,
        double timeStep,
        double currentTime,
        const std::vector<VectorReal>& x_hist) override {
        (void)b;
        (void)x;
        (void)timeStep;
        (void)currentTime;
        (void)x_hist;
        stamp(J, transconductance_);
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        (void)b;
        (void)omega;
        (void)x_dc;
        stamp(J, std::complex<double>{transconductance_, 0.0});
    }

private:
    template <typename MatrixT, typename ValueT>
    void stamp(MatrixT& J, ValueT gain) const {
        J.add(nodePos_, ctrlPos_, gain);
        J.add(nodePos_, ctrlNeg_, -gain);
        J.add(nodeNeg_, ctrlPos_, -gain);
        J.add(nodeNeg_, ctrlNeg_, gain);
    }

    int nodePos_;
    int nodeNeg_;
    int ctrlPos_;
    int ctrlNeg_;
    double transconductance_;
};

class VoltageControlledVoltageSource : public Device {
public:
    VoltageControlledVoltageSource(
        const std::string& name,
        int nodePos,
        int nodeNeg,
        int ctrlPos,
        int ctrlNeg,
        double gain,
        int branchIndex = -1)
        : Device(name),
          nodePos_(nodePos),
          nodeNeg_(nodeNeg),
          ctrlPos_(ctrlPos),
          ctrlNeg_(ctrlNeg),
          gain_(gain),
          branchIndex_(branchIndex) {}

    void setBranchIndex(int index) { branchIndex_ = index; }
    int getBranchIndex() const { return branchIndex_; }

    void dcStamp(
        SparseMatrixReal& J,
        VectorReal& b,
        const VectorReal& x,
        double timeStep,
        double currentTime,
        const std::vector<VectorReal>& x_hist) override {
        (void)b;
        (void)x;
        (void)timeStep;
        (void)currentTime;
        (void)x_hist;
        if (branchIndex_ < 0) return;
        stamp(J, gain_);
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        (void)b;
        (void)omega;
        (void)x_dc;
        if (branchIndex_ < 0) return;
        stamp(J, std::complex<double>{gain_, 0.0});
    }

private:
    template <typename MatrixT, typename ValueT>
    void stamp(MatrixT& J, ValueT gain) const {
        J.add(nodePos_, branchIndex_, ValueT{1.0});
        J.add(nodeNeg_, branchIndex_, ValueT{-1.0});
        J.add(branchIndex_, nodePos_, ValueT{1.0});
        J.add(branchIndex_, nodeNeg_, ValueT{-1.0});
        J.add(branchIndex_, ctrlPos_, -gain);
        J.add(branchIndex_, ctrlNeg_, gain);
    }

    int nodePos_;
    int nodeNeg_;
    int ctrlPos_;
    int ctrlNeg_;
    double gain_;
    int branchIndex_;
};

class CurrentControlledCurrentSource : public Device {
public:
    CurrentControlledCurrentSource(
        const std::string& name,
        int nodePos,
        int nodeNeg,
        std::string controlSource,
        double gain)
        : Device(name),
          nodePos_(nodePos),
          nodeNeg_(nodeNeg),
          controlSource_(std::move(controlSource)),
          gain_(gain) {}

    const std::string& getControlSource() const { return controlSource_; }
    void setControlBranchIndex(int index) { controlBranchIndex_ = index; }
    int getControlBranchIndex() const { return controlBranchIndex_; }

    void dcStamp(
        SparseMatrixReal& J,
        VectorReal& b,
        const VectorReal& x,
        double timeStep,
        double currentTime,
        const std::vector<VectorReal>& x_hist) override {
        (void)b;
        (void)x;
        (void)timeStep;
        (void)currentTime;
        (void)x_hist;
        if (controlBranchIndex_ < 0) return;
        stamp(J, gain_);
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        (void)b;
        (void)omega;
        (void)x_dc;
        if (controlBranchIndex_ < 0) return;
        stamp(J, std::complex<double>{gain_, 0.0});
    }

private:
    template <typename MatrixT, typename ValueT>
    void stamp(MatrixT& J, ValueT gain) const {
        J.add(nodePos_, controlBranchIndex_, gain);
        J.add(nodeNeg_, controlBranchIndex_, -gain);
    }

    int nodePos_;
    int nodeNeg_;
    std::string controlSource_;
    double gain_;
    int controlBranchIndex_ = -1;
};

class CurrentControlledVoltageSource : public Device {
public:
    CurrentControlledVoltageSource(
        const std::string& name,
        int nodePos,
        int nodeNeg,
        std::string controlSource,
        double transresistance,
        int branchIndex = -1)
        : Device(name),
          nodePos_(nodePos),
          nodeNeg_(nodeNeg),
          controlSource_(std::move(controlSource)),
          transresistance_(transresistance),
          branchIndex_(branchIndex) {}

    const std::string& getControlSource() const { return controlSource_; }
    void setControlBranchIndex(int index) { controlBranchIndex_ = index; }
    int getControlBranchIndex() const { return controlBranchIndex_; }
    void setBranchIndex(int index) { branchIndex_ = index; }
    int getBranchIndex() const { return branchIndex_; }

    void dcStamp(
        SparseMatrixReal& J,
        VectorReal& b,
        const VectorReal& x,
        double timeStep,
        double currentTime,
        const std::vector<VectorReal>& x_hist) override {
        (void)b;
        (void)x;
        (void)timeStep;
        (void)currentTime;
        (void)x_hist;
        if (branchIndex_ < 0 || controlBranchIndex_ < 0) return;
        stamp(J, transresistance_);
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        (void)b;
        (void)omega;
        (void)x_dc;
        if (branchIndex_ < 0 || controlBranchIndex_ < 0) return;
        stamp(J, std::complex<double>{transresistance_, 0.0});
    }

private:
    template <typename MatrixT, typename ValueT>
    void stamp(MatrixT& J, ValueT transresistance) const {
        J.add(nodePos_, branchIndex_, ValueT{1.0});
        J.add(nodeNeg_, branchIndex_, ValueT{-1.0});
        J.add(branchIndex_, nodePos_, ValueT{1.0});
        J.add(branchIndex_, nodeNeg_, ValueT{-1.0});
        J.add(branchIndex_, controlBranchIndex_, -transresistance);
    }

    int nodePos_;
    int nodeNeg_;
    std::string controlSource_;
    double transresistance_;
    int branchIndex_;
    int controlBranchIndex_ = -1;
};

} // namespace gspice

#endif // GSPICE_CONTROLLED_SOURCE_HPP
