#ifndef GSPICE_OSDI_H
#define GSPICE_OSDI_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

// OSDI Parameter Types
typedef enum {
    OSDI_PARA_TY_REAL = 0,
    OSDI_PARA_TY_INT = 1,
    OSDI_PARA_TY_STR = 2
} OsdiParameterType;

// Model Metadata
typedef struct {
    const char* name;
    uint32_t num_terminals;
    const char** terminal_names;
} OsdiModelDescriptor;

// Evaluation Input/Output
typedef struct {
    const double* voltages;     // Input: node voltages
    double* currents;           // Output: terminal currents
    double* charges;            // Output: terminal charges
    double* jacobian;           // Output: dI/dV and dQ/dV (conductances/capacitances)
} OsdiEvaluationData;

// Function pointers for model logic
typedef void (*OsdiEvaluateFunc)(void* instance_data, OsdiEvaluationData* data);
typedef void* (*OsdiCreateInstanceFunc)(void* model_data);

// Descriptor for a compiled model
typedef struct {
    const char* model_name;
    uint32_t version_major;
    uint32_t version_minor;
    OsdiModelDescriptor metadata;
    OsdiEvaluateFunc evaluate;
    OsdiCreateInstanceFunc create_instance;
} OsdiDescriptor;

#ifdef __cplusplus
}
#endif

#endif // GSPICE_OSDI_H
