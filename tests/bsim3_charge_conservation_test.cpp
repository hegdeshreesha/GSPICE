#include "devices/bsim3_model.hpp"

#include <cassert>
#include <cmath>

class TestBsim3Charge final : public gspice::Bsim3Mosfet {
public:
    using gspice::Bsim3Mosfet::Bsim3Mosfet;
    void acStamp(gspice::SparseMatrixComplex&, gspice::VectorComplex&, double,
                 const gspice::VectorReal&) override {}
};

int main() {
    TestBsim3Charge device("M1", 0, 1, 2, 3);
    gspice::VectorReal x(4);
    x[0] = 0.2;
    x[1] = 0.9;
    gspice::DaeRequest request;
    request.staticResidual = false;
    request.staticJacobian = false;
    request.dynamicResidual = true;
    request.dynamicJacobian = true;
    gspice::DaeEvaluation evaluation;
    assert(device.evaluateDae(x, request, evaluation));
    assert(evaluation.finite());
    double charge_sum = 0.0;
    for (const auto& term : evaluation.dynamicResidual) charge_sum += term.value;
    assert(std::abs(charge_sum) < 1.0e-24);
}
