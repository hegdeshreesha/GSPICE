#include "devices/gdi_device.hpp"
#include "devices/bsim4_model.hpp"
#include "devices/gmc_juncap_express.hpp"
#include "gsdi_device_adapter.hpp"

#include <cassert>
#include <cmath>

int main() {
    gspice::GmcJuncapExpressModule module;
    auto instance = module.createInstance(
        {{"PHITD", 0.026}, {"ISATFOR1", 1e-14}, {"MFOR1", 1.0}},
        {{"VJ", 0.2}, {"DVJ", 1.0}});

    gspice::GdiDevice device("D1", std::move(instance), {1, 0});
    gspice::VectorReal x(2);
    x[0] = 0.0;
    x[1] = 0.6;

    gspice::DaeRequest request;
    request.staticResidual = true;
    request.staticJacobian = true;
    request.dynamicResidual = true;
    request.dynamicJacobian = true;

    gspice::DaeEvaluation evaluation;
    assert(device.evaluateDae(x, request, evaluation));
    assert(evaluation.staticResidual.size() == 2);
    assert(evaluation.staticJacobian.size() == 4);
    assert(evaluation.dynamicResidual.size() == 2);
    assert(evaluation.dynamicJacobian.size() == 4);
    assert(evaluation.staticResidual[0].equation == 1);
    assert(evaluation.staticResidual[1].equation == 0);
    assert(std::abs(evaluation.staticResidual[0].value + evaluation.staticResidual[1].value) < 1e-30);
    assert(evaluation.finite());
    assert(device.getNoisePSD(1.0, x) > 0.0);

    const auto model = gspice::Bsim4ParameterSet::from({{"LEVEL", 54.0}, {"U0", 500.0}})
                           .prepare(1.0e-6, 1.0e-6);
    auto bsim = std::make_unique<gspice::Bsim4Mosfet>("M1", 0, 1, 2, 3, model);
    gspice::GdiDevice bsim_device(
        "M1",
        std::make_unique<gspice::GsdiDaeDeviceAdapter>(std::move(bsim), 4),
        std::vector<int>{3, 2, 1, 0});
    gspice::VectorReal xb(4);
    xb[3] = 0.8;
    xb[2] = 0.9;
    gspice::DaeEvaluation bsim_eval;
    assert(bsim_device.evaluateDae(xb, request, bsim_eval));
    assert(bsim_eval.finite());
    assert(bsim_eval.staticResidual[0].equation == 3);
    gspice::SparseMatrixComplex ac(4);
    gspice::VectorComplex rhs(4);
    bsim_device.acStamp(ac, rhs, 1.0e9, xb);
    const auto dense = ac.toDense();
    bool has_imaginary = false;
    for (int row = 0; row < 4; ++row)
        for (int column = 0; column < 4; ++column)
            has_imaginary = has_imaginary || std::abs(dense(row, column).imag()) > 0.0;
    assert(has_imaginary);
    return 0;
}
