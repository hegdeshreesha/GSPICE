#include "devices/bsim4_model.hpp"

#include <cassert>
#include <cmath>

int main() {
    const auto model = gspice::Bsim4ParameterSet::from({
        {"LEVEL", 54.0}, {"U0", 500.0}, {"CGSO", 1.0e-10},
        {"CGDO", 1.0e-10}, {"CGBO", 5.0e-11}, {"CJ", 1.0e-3},
        {"CJSW", 2.0e-10}, {"PB", 0.8}, {"PBSW", 0.8}
    })
                           .prepare(1.0e-6, 1.0e-6);
    gspice::Bsim4Mosfet device("M1", 0, 1, 2, 3, model);
    gspice::VectorReal x(4);
    x[0] = 0.8;
    x[1] = 0.9;
    gspice::DaeRequest request;
    request.staticResidual = true;
    request.staticJacobian = true;
    request.dynamicResidual = true;
    request.dynamicJacobian = true;
    gspice::DaeEvaluation evaluation;
    assert(device.evaluateDae(x, request, evaluation));
    assert(evaluation.finite());
    double charge_sum = 0.0;
    for (const auto& term : evaluation.dynamicResidual) charge_sum += term.value;
    assert(std::abs(charge_sum) < 1.0e-24);
    for (int column = 0; column < 4; ++column) {
        double derivative_sum = 0.0;
        for (const auto& term : evaluation.dynamicJacobian)
            if (term.unknown == column) derivative_sum += term.value;
        assert(std::abs(derivative_sum) < 1.0e-18);
    }

    gspice::SparseMatrixComplex ac(4);
    gspice::VectorComplex rhs(4);
    device.acStamp(ac, rhs, 1.0e9, x);
    const auto dense = ac.toDense();
    bool has_imaginary = false;
    for (int row = 0; row < 4; ++row)
        for (int column = 0; column < 4; ++column)
            has_imaginary = has_imaginary || std::abs(dense(row, column).imag()) > 0.0;
    assert(has_imaginary);
}
