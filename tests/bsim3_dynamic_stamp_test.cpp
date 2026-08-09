#include "devices/bsim3_model.hpp"

#include <cassert>
#include <cmath>

int main() {
    gspice::Bsim3Mosfet device("M1", 0, 1, 2, 3);
    gspice::VectorReal x(4);
    x[0] = 0.2;
    x[1] = 0.9;

    gspice::SparseMatrixComplex ac(4);
    gspice::VectorComplex ac_rhs(4);
    device.acStamp(ac, ac_rhs, 2.0e9, x);
    const auto ac_dense = ac.toDense();
    bool has_cap_susceptance = false;
    for (int row = 0; row < 4; ++row)
        for (int col = 0; col < 4; ++col)
            has_cap_susceptance = has_cap_susceptance ||
                std::abs(ac_dense(row, col).imag()) > 0.0;
    assert(has_cap_susceptance);

    gspice::VectorReal previous = x;
    previous[1] = 0.8;
    std::vector<gspice::VectorReal> history{previous};
    gspice::TransientContext context;
    context.timeStep = 1.0e-12;
    context.currentTime = context.timeStep;
    context.method = gspice::TransientIntegrationMethod::BackwardEuler;
    context.a0 = 1.0 / context.timeStep;
    context.a1 = -context.a0;
    context.xHistory = &history;

    gspice::SparseMatrixReal tran(4);
    gspice::VectorReal rhs(4);
    device.tranStamp(tran, rhs, x, context);
    const auto tran_dense = tran.toDense();
    bool has_dynamic_stamp = false;
    for (int row = 0; row < 4; ++row)
        for (int col = 0; col < 4; ++col)
            has_dynamic_stamp = has_dynamic_stamp || std::abs(tran_dense(row, col)) > 0.0;
    assert(has_dynamic_stamp);
}
