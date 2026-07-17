#ifndef GSPICE_CURRENT_SOURCE_HPP
#define GSPICE_CURRENT_SOURCE_HPP

#include "device.hpp"
#include "devices/voltage_source.hpp"
#include <string>
#include <vector>
#include <cmath>

namespace gspice {

class CurrentSource : public Device {
public:
    /**
     * @param name Name of source
     * @param nodePos Positive node (current flows FROM this node)
     * @param nodeNeg Negative node (current flows TO this node)
     * @param dcValue DC current value
     */
    CurrentSource(const std::string& name, int nodePos, int nodeNeg, double dcValue)
        : Device(name), nodePos_(nodePos), nodeNeg_(nodeNeg), dcValue_(dcValue) {}

    void setPulse(const VoltageSource::PulseParams& pulse) {
        waveformType_ = VoltageSource::WaveformType::PULSE;
        pulse_ = pulse;
    }

    void setSin(const VoltageSource::SinParams& sin) {
        waveformType_ = VoltageSource::WaveformType::SIN;
        sin_ = sin;
    }

    void setPwl(const std::vector<double>& times, const std::vector<double>& values) {
        if (times.empty() || values.empty() || times.size() != values.size()) return;
        waveformType_ = VoltageSource::WaveformType::PWL;
        pwlTimes_ = times;
        pwlValues_ = values;
    }

    void setAcMagnitude(double value) {
        acMagnitude_ = value;
    }

    void setDcValue(double value) {
        dcValue_ = value;
    }

    double getDcValue() const {
        return dcValue_;
    }

    void collectBreakpoints(double t_stop, std::vector<double>& points) const override {
        if (t_stop <= 0.0) return;
        auto add = [&](double t) {
            if (t >= 0.0 && t <= t_stop) points.push_back(t);
        };
        if (waveformType_ == VoltageSource::WaveformType::PULSE) {
            const double tr = (pulse_.tr > 0.0) ? pulse_.tr : 1e-18;
            const double tf = (pulse_.tf > 0.0) ? pulse_.tf : 1e-18;
            const double period = pulse_.per > 0.0 ? pulse_.per : t_stop + 1.0;
            for (double base = pulse_.td; base <= t_stop; base += period) {
                add(base);
                add(base + tr);
                add(base + tr + pulse_.pw);
                add(base + tr + pulse_.pw + tf);
                if (pulse_.per <= 0.0) break;
            }
        } else if (waveformType_ == VoltageSource::WaveformType::PWL) {
            for (double t : pwlTimes_) add(t);
        } else if (waveformType_ == VoltageSource::WaveformType::SIN && sin_.td > 0.0) {
            add(sin_.td);
        }
    }

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, double timeStep, double currentTime, const std::vector<VectorReal>& x_hist) override {
        // Current Source RHS Stamp:
        // Ax = b. A current source is a constant 'b'.
        // Flowing out of nodePos: RHS[nodePos] -= dcValue
        // Flowing into nodeNeg:  RHS[nodeNeg] += dcValue
        const double current = evaluateAt(currentTime);
        b.add(nodePos_, -current);
        b.add(nodeNeg_, current);
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        b.add(nodePos_, {-acMagnitude_, 0.0});
        b.add(nodeNeg_, {acMagnitude_, 0.0});
    }

    void hbStamp(SparseMatrixReal& J, VectorReal& b, double f_fund, int n_harms, const VectorReal& x_hb) override {
        int K = 2 * n_harms + 1;
        // RHS: DC component at index 0
        b.add(nodePos_ * K, -dcValue_);
        b.add(nodeNeg_ * K, dcValue_);
    }

private:
    static constexpr double PI_ = 3.14159265358979323846;

    double evaluateAt(double time) const {
        if (time <= 0.0) return dcValue_;
        if (waveformType_ == VoltageSource::WaveformType::PULSE) return evalPulse(time);
        if (waveformType_ == VoltageSource::WaveformType::SIN) return evalSin(time);
        if (waveformType_ == VoltageSource::WaveformType::PWL) return evalPwl(time);
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
        if (local < tr) return pulse_.v1 + (pulse_.v2 - pulse_.v1) * (local / tr);
        if (local < tr + pulse_.pw) return pulse_.v2;
        if (local < tr + pulse_.pw + tf) return pulse_.v2 + (pulse_.v1 - pulse_.v2) * ((local - tr - pulse_.pw) / tf);
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

    int nodePos_;
    int nodeNeg_;
    double dcValue_;
    double acMagnitude_ = 1.0;
    VoltageSource::WaveformType waveformType_ = VoltageSource::WaveformType::DC;
    VoltageSource::PulseParams pulse_;
    VoltageSource::SinParams sin_;
    std::vector<double> pwlTimes_;
    std::vector<double> pwlValues_;
};

} // namespace gspice

#endif // GSPICE_CURRENT_SOURCE_HPP
