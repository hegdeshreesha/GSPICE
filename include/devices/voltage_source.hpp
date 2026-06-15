#ifndef GSPICE_VOLTAGE_SOURCE_HPP
#define GSPICE_VOLTAGE_SOURCE_HPP

#include "device.hpp"
#include <string>
#include <vector>
#include <cmath>

namespace gspice {

class VoltageSource : public Device {
public:
    enum class WaveformType {
        DC,
        PULSE,
        SIN,
        PWL
    };

    struct PulseParams {
        double v1 = 0.0;
        double v2 = 0.0;
        double td = 0.0;
        double tr = 1e-12;
        double tf = 1e-12;
        double pw = 0.0;
        double per = 0.0;
    };

    struct SinParams {
        double vo = 0.0;
        double va = 0.0;
        double freq = 0.0;
        double td = 0.0;
        double theta = 0.0;
        double phase_deg = 0.0;
    };

    VoltageSource(const std::string& name, int nodePos, int nodeNeg, double voltage, int branchIndex = -1)
        : Device(name), nodePos_(nodePos), nodeNeg_(nodeNeg), dcValue_(voltage), branchIndex_(branchIndex) {}

    void setPulse(const PulseParams& pulse) {
        waveformType_ = WaveformType::PULSE;
        pulse_ = pulse;
    }

    void setSin(const SinParams& sin) {
        waveformType_ = WaveformType::SIN;
        sin_ = sin;
    }

    void setPwl(const std::vector<double>& times, const std::vector<double>& values) {
        if (times.empty() || values.empty() || times.size() != values.size()) return;
        waveformType_ = WaveformType::PWL;
        pwlTimes_ = times;
        pwlValues_ = values;
    }

    void setAcMagnitude(double value) {
        acMagnitude_ = value;
    }

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, double timeStep, const std::vector<VectorReal>& x_hist) override {
        if (branchIndex_ < 0) return;
        J.add(nodePos_, branchIndex_, 1.0);
        J.add(nodeNeg_, branchIndex_, -1.0);
        J.add(branchIndex_, nodePos_, 1.0);
        J.add(branchIndex_, nodeNeg_, -1.0);
        b.add(branchIndex_, evaluateSourceValue(timeStep, x_hist));
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        if (branchIndex_ < 0) return;
        J.add(nodePos_, branchIndex_, {1.0, 0.0});
        J.add(nodeNeg_, branchIndex_, {-1.0, 0.0});
        J.add(branchIndex_, nodePos_, {1.0, 0.0});
        J.add(branchIndex_, nodeNeg_, {-1.0, 0.0});
        b.add(branchIndex_, {acMagnitude_, 0.0});
    }

    void pacStamp(SparseMatrixReal& J, VectorReal& b, double f_in, double f_fund, int n_harms, const VectorReal& x_periodic) override {
        hbStamp(J, b, f_fund, n_harms, x_periodic); 
        int K = 2 * n_harms + 1;
        b.add(branchIndex_ * K, 1.0);
    }

    void hbStamp(SparseMatrixReal& J, VectorReal& b, double f_fund, int n_harms, const VectorReal& x_hb) override {
        if (branchIndex_ < 0) return;
        int K = 2 * n_harms + 1;
        for (int k = 0; k < K; ++k) {
            J.add(nodePos_ * K + k, branchIndex_ * K + k, 1.0);
            J.add(nodeNeg_ * K + k, branchIndex_ * K + k, -1.0);
            J.add(branchIndex_ * K + k, nodePos_ * K + k, 1.0);
            J.add(branchIndex_ * K + k, nodeNeg_ * K + k, -1.0);
        }
        b.add(branchIndex_ * K, dcValue_);
    }

    void setBranchIndex(int index) { branchIndex_ = index; }
    int getBranchIndex() const { return branchIndex_; }

private:
    static constexpr double PI_ = 3.14159265358979323846;

    double evaluateSourceValue(double timeStep, const std::vector<VectorReal>& x_hist) const {
        if (timeStep <= 0.0 || x_hist.empty()) {
            return dcValue_;
        }
        const double t = timeStep * static_cast<double>(x_hist.size());

        if (waveformType_ == WaveformType::PULSE) return evalPulse(t);
        if (waveformType_ == WaveformType::SIN) return evalSin(t);
        if (waveformType_ == WaveformType::PWL) return evalPwl(t);
        return dcValue_;
    }

    double evalPulse(double t) const {
        if (t < pulse_.td) return pulse_.v1;
        double tr = (pulse_.tr > 0.0) ? pulse_.tr : 1e-18;
        double tf = (pulse_.tf > 0.0) ? pulse_.tf : 1e-18;

        double local = t - pulse_.td;
        if (pulse_.per > 0.0) {
            local = std::fmod(local, pulse_.per);
            if (local < 0.0) local += pulse_.per;
        }

        if (local < tr) {
            return pulse_.v1 + (pulse_.v2 - pulse_.v1) * (local / tr);
        }
        if (local < tr + pulse_.pw) {
            return pulse_.v2;
        }
        if (local < tr + pulse_.pw + tf) {
            return pulse_.v2 + (pulse_.v1 - pulse_.v2) * ((local - tr - pulse_.pw) / tf);
        }
        return pulse_.v1;
    }

    double evalSin(double t) const {
        if (t < sin_.td) return sin_.vo;
        double tau = t - sin_.td;
        double phase = (sin_.phase_deg * PI_ / 180.0) + (2.0 * PI_ * sin_.freq * tau);
        double env = (sin_.theta == 0.0) ? 1.0 : std::exp(-sin_.theta * tau);
        return sin_.vo + sin_.va * env * std::sin(phase);
    }

    double evalPwl(double t) const {
        if (pwlTimes_.empty() || pwlValues_.empty()) return dcValue_;
        if (t <= pwlTimes_.front()) return pwlValues_.front();
        for (size_t i = 1; i < pwlTimes_.size(); ++i) {
            if (t <= pwlTimes_[i]) {
                double dt = pwlTimes_[i] - pwlTimes_[i - 1];
                if (dt <= 0.0) return pwlValues_[i];
                double alpha = (t - pwlTimes_[i - 1]) / dt;
                return pwlValues_[i - 1] + alpha * (pwlValues_[i] - pwlValues_[i - 1]);
            }
        }
        return pwlValues_.back();
    }

    int nodePos_, nodeNeg_;
    double dcValue_ = 0.0;
    double acMagnitude_ = 1.0;
    int branchIndex_;
    WaveformType waveformType_ = WaveformType::DC;
    PulseParams pulse_;
    SinParams sin_;
    std::vector<double> pwlTimes_;
    std::vector<double> pwlValues_;
};

} // namespace gspice

#endif // GSPICE_VOLTAGE_SOURCE_HPP
