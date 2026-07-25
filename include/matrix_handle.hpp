#ifndef GSPICE_MATRIX_HANDLE_HPP
#define GSPICE_MATRIX_HANDLE_HPP

#include <cstddef>

namespace gspice {

struct MatrixEntryHandle {
    int row = -1;
    int col = -1;
    double* valuePtr = nullptr;
};

struct DeviceMatrixHandles {
    MatrixEntryHandle g_rr;
    MatrixEntryHandle g_rc;
    MatrixEntryHandle g_cr;
    MatrixEntryHandle g_cc;

    void stamp(double val_rr, double val_rc, double val_cr, double val_cc) const {
        if (g_rr.valuePtr) *g_rr.valuePtr += val_rr;
        if (g_rc.valuePtr) *g_rc.valuePtr += val_rc;
        if (g_cr.valuePtr) *g_cr.valuePtr += val_cr;
        if (g_cc.valuePtr) *g_cc.valuePtr += val_cc;
    }
};

} // namespace gspice

#endif // GSPICE_MATRIX_HANDLE_HPP
