#ifndef GSPICE_OSDI_DEVICE_HPP
#define GSPICE_OSDI_DEVICE_HPP

#include "device.hpp"
#include "osdi.h"
#include <vector>

namespace gspice {

class OSDIDevice : public Device {
public:
    OSDIDevice(const std::string& name, const OsdiDescriptor& desc, const std::vector<int>& nodes)
        : Device(name), desc_(desc), nodes_(nodes) {
        // Create an instance of the model logic
        instance_data_ = desc_.create_instance(nullptr);
    }

    void dcStamp(SparseMatrixReal& J, VectorReal& b, const VectorReal& x, double timeStep, const std::vector<VectorReal>& x_hist) override {
        int n = static_cast<int>(nodes_.size());
        std::vector<double> voltages(n);
        for (int i = 0; i < n; ++i) {
            voltages[i] = (nodes_[i] >= 0) ? x[nodes_[i]] : 0.0;
        }

        std::vector<double> currents(n, 0.0);
        std::vector<double> charges(n, 0.0);
        std::vector<double> jacobian(n * n, 0.0);

        OsdiEvaluationData eval_data = {
            voltages.data(),
            currents.data(),
            charges.data(),
            jacobian.data()
        };

        // Call the compiled OSDI model logic
        desc_.evaluate(instance_data_, &eval_data);

        // 1. Stamp Jacobian (conductances)
        for (int i = 0; i < n; ++i) {
            if (nodes_[i] < 0) continue;
            for (int j = 0; j < n; ++j) {
                if (nodes_[j] < 0) continue;
                J.add(nodes_[i], nodes_[j], jacobian[i * n + j]);
            }
        }

        // 2. Stamp RHS (currents)
        // Ieq_i = I_i - sum_j(g_ij * V_j)
        for (int i = 0; i < n; ++i) {
            if (nodes_[i] < 0) continue;
            double Ieq = currents[i];
            for (int j = 0; j < n; ++j) {
                Ieq -= jacobian[i * n + j] * voltages[j];
            }
            b.add(nodes_[i], -Ieq);
        }
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {
        // Linearize at DC OP and stamp small-signal parameters
        // (Similar to dcStamp but with omega and complex matrix)
    }

private:
    OsdiDescriptor desc_;
    std::vector<int> nodes_;
    void* instance_data_;
};

} // namespace gspice

#endif // GSPICE_OSDI_DEVICE_HPP
