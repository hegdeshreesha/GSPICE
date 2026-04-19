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

        // Standard Newton-Raphson update: J * V_next = J * V_k - I_total
        // RHS_entry = sum_j(g_ij * V_j) - I_i_total
        
        for (int i = 0; i < n; ++i) {
            if (nodes_[i] < 0) continue;
            
            // 1. Stamp Jacobian (conductances)
            for (int j = 0; j < n; ++j) {
                if (nodes_[j] < 0) continue;
                J.add(nodes_[i], nodes_[j], jacobian[i * n + j]);
            }

            // 2. Stamp RHS (current equivalents)
            double rhs_val = -currents[i];
            for (int j = 0; j < n; ++j) {
                rhs_val += jacobian[i * n + j] * voltages[j];
            }
            b.add(nodes_[i], rhs_val);
        }
    }

    void acStamp(SparseMatrixComplex& J, VectorComplex& b, double omega, const VectorReal& x_dc) override {}

private:
    OsdiDescriptor desc_;
    std::vector<int> nodes_;
    void* instance_data_;
};

} // namespace gspice

#endif // GSPICE_OSDI_DEVICE_HPP
