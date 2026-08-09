#include "devices/bsim3_model.hpp"

#include <cassert>
#include <cmath>

class TestBsim3 final : public gspice::Bsim3Mosfet {
public:
    using gspice::Bsim3Mosfet::Bsim3Mosfet;
    void acStamp(gspice::SparseMatrixComplex&, gspice::VectorComplex&, double,
                 const gspice::VectorReal&) override {}
};

static gspice::DaeEvaluation evaluate(gspice::Bsim3Mosfet& device, const gspice::VectorReal& x) {
    gspice::DaeRequest request;
    gspice::DaeEvaluation result;
    assert(device.evaluateDae(x, request, result));
    assert(result.finite());
    return result;
}

int main() {
    TestBsim3 device("M1", 0, 1, 2, 3);
    gspice::VectorReal x(4);
    x[0] = 0.3;
    x[1] = 1.5;
    x[2] = 0.1;
    x[3] = 0.0;

    const auto base = evaluate(device, x);
    double jac[4][4]{};
    for (const auto& term : base.staticJacobian) jac[term.equation][term.unknown] = term.value;

    for (int column = 0; column < 4; ++column) {
        const double h = 1.0e-6;
        gspice::VectorReal plus = x;
        gspice::VectorReal minus = x;
        plus[column] += h;
        minus[column] -= h;
        const auto p = evaluate(device, plus);
        const auto m = evaluate(device, minus);
        for (int row = 0; row < 4; ++row) {
            const double fp = p.staticResidual[row].value;
            const double fm = m.staticResidual[row].value;
            const double numerical = (fp - fm) / (2.0 * h);
            assert(std::abs(jac[row][column] - numerical) < 1.0e-8 * std::max(1.0, std::abs(numerical)));
        }
    }
}
