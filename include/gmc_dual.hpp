#ifndef GSPICE_GMC_DUAL_HPP
#define GSPICE_GMC_DUAL_HPP

#include <array>
#include <cmath>
#include <cstddef>
#include <type_traits>

namespace gspice {

template <std::size_t N>
struct GmcDual;

// Detection trait for GmcDual<N> used by GMC-generated and hand-written model
// code to dispatch between dual-number and plain-double evaluation.
template <typename T>
struct is_gmc_dual : std::false_type {};

template <std::size_t N>
struct is_gmc_dual<GmcDual<N>> : std::true_type {};

template <typename T>
inline constexpr bool is_gmc_dual_v = is_gmc_dual<T>::value;

// Forward-mode automatic differentiation over N independent voltage
// differences. GMC-generated compact-model code and hand-ported model cores
// (PSP, BSIM, JUNCAP) evaluate currents, charges, and their analytic
// Jacobians from a single expression. N is the number of local unknowns the
// device participates in (terminals plus kept internal/collapsible nodes).
template <std::size_t N>
struct GmcDual {
    double value = 0.0;
    std::array<double, N> derivative{};

    GmcDual() = default;
    GmcDual(double input) : value(input) {}

    static constexpr std::size_t size() { return N; }

    static GmcDual constant(double input) {
        GmcDual result;
        result.value = input;
        return result;
    }

    static GmcDual variable(double input, std::size_t index) {
        GmcDual result = constant(input);
        if (index < N) result.derivative[index] = 1.0;
        return result;
    }
};

template <std::size_t N>
inline GmcDual<N>& operator+=(GmcDual<N>& lhs, const GmcDual<N>& rhs) {
    lhs.value += rhs.value;
    for (std::size_t i = 0; i < N; ++i) lhs.derivative[i] += rhs.derivative[i];
    return lhs;
}

template <std::size_t N>
inline GmcDual<N>& operator-=(GmcDual<N>& lhs, const GmcDual<N>& rhs) {
    lhs.value -= rhs.value;
    for (std::size_t i = 0; i < N; ++i) lhs.derivative[i] -= rhs.derivative[i];
    return lhs;
}

template <std::size_t N>
inline GmcDual<N>& operator*=(GmcDual<N>& lhs, const GmcDual<N>& rhs) {
    lhs = lhs * rhs;
    return lhs;
}

template <std::size_t N>
inline GmcDual<N>& operator/=(GmcDual<N>& lhs, const GmcDual<N>& rhs) {
    lhs = lhs / rhs;
    return lhs;
}

template <std::size_t N>
inline GmcDual<N>& operator+=(GmcDual<N>& lhs, double rhs) {
    lhs.value += rhs;
    return lhs;
}

template <std::size_t N>
inline GmcDual<N>& operator-=(GmcDual<N>& lhs, double rhs) {
    lhs.value -= rhs;
    return lhs;
}

template <std::size_t N>
inline GmcDual<N>& operator*=(GmcDual<N>& lhs, double rhs) {
    lhs.value *= rhs;
    for (double& derivative : lhs.derivative) derivative *= rhs;
    return lhs;
}

template <std::size_t N>
inline GmcDual<N>& operator/=(GmcDual<N>& lhs, double rhs) {
    lhs.value /= rhs;
    for (double& derivative : lhs.derivative) derivative /= rhs;
    return lhs;
}

template <std::size_t N>
inline GmcDual<N> operator+(GmcDual<N> lhs, const GmcDual<N>& rhs) {
    lhs.value += rhs.value;
    for (std::size_t i = 0; i < N; ++i) lhs.derivative[i] += rhs.derivative[i];
    return lhs;
}

template <std::size_t N>
inline GmcDual<N> operator-(GmcDual<N> lhs, const GmcDual<N>& rhs) {
    lhs.value -= rhs.value;
    for (std::size_t i = 0; i < N; ++i) lhs.derivative[i] -= rhs.derivative[i];
    return lhs;
}

template <std::size_t N>
inline GmcDual<N> operator-(GmcDual<N> value) {
    value.value = -value.value;
    for (double& derivative : value.derivative) derivative = -derivative;
    return value;
}

template <std::size_t N>
inline GmcDual<N> operator*(const GmcDual<N>& lhs, const GmcDual<N>& rhs) {
    GmcDual<N> result;
    result.value = lhs.value * rhs.value;
    for (std::size_t i = 0; i < N; ++i) {
        result.derivative[i] = lhs.derivative[i] * rhs.value + lhs.value * rhs.derivative[i];
    }
    return result;
}

template <std::size_t N>
inline GmcDual<N> operator/(const GmcDual<N>& lhs, const GmcDual<N>& rhs) {
    GmcDual<N> result;
    result.value = lhs.value / rhs.value;
    const double scale = rhs.value * rhs.value;
    for (std::size_t i = 0; i < N; ++i) {
        result.derivative[i] = (lhs.derivative[i] * rhs.value - lhs.value * rhs.derivative[i]) / scale;
    }
    return result;
}

template <std::size_t N>
inline GmcDual<N> operator+(GmcDual<N> lhs, double rhs) { lhs.value += rhs; return lhs; }
template <std::size_t N>
inline GmcDual<N> operator+(double lhs, GmcDual<N> rhs) { rhs.value += lhs; return rhs; }
template <std::size_t N>
inline GmcDual<N> operator-(GmcDual<N> lhs, double rhs) { lhs.value -= rhs; return lhs; }
template <std::size_t N>
inline GmcDual<N> operator-(double lhs, const GmcDual<N>& rhs) { return GmcDual<N>::constant(lhs) - rhs; }
template <std::size_t N>
inline GmcDual<N> operator*(GmcDual<N> lhs, double rhs) {
    lhs.value *= rhs;
    for (double& derivative : lhs.derivative) derivative *= rhs;
    return lhs;
}
template <std::size_t N>
inline GmcDual<N> operator*(double lhs, GmcDual<N> rhs) { return rhs * lhs; }
template <std::size_t N>
inline GmcDual<N> operator/(GmcDual<N> lhs, double rhs) { return lhs * (1.0 / rhs); }
template <std::size_t N>
inline GmcDual<N> operator/(double lhs, const GmcDual<N>& rhs) {
    return GmcDual<N>::constant(lhs) / rhs;
}

template <std::size_t N>
inline bool operator<(const GmcDual<N>& lhs, const GmcDual<N>& rhs) { return lhs.value < rhs.value; }
template <std::size_t N>
inline bool operator>(const GmcDual<N>& lhs, const GmcDual<N>& rhs) { return lhs.value > rhs.value; }
template <std::size_t N>
inline bool operator<=(const GmcDual<N>& lhs, const GmcDual<N>& rhs) { return lhs.value <= rhs.value; }
template <std::size_t N>
inline bool operator>=(const GmcDual<N>& lhs, const GmcDual<N>& rhs) { return lhs.value >= rhs.value; }
template <std::size_t N>
inline bool operator==(const GmcDual<N>& lhs, const GmcDual<N>& rhs) { return lhs.value == rhs.value; }
template <std::size_t N>
inline bool operator!=(const GmcDual<N>& lhs, const GmcDual<N>& rhs) { return lhs.value != rhs.value; }
template <std::size_t N>
inline bool operator<(const GmcDual<N>& lhs, double rhs) { return lhs.value < rhs; }
template <std::size_t N>
inline bool operator>(const GmcDual<N>& lhs, double rhs) { return lhs.value > rhs; }
template <std::size_t N>
inline bool operator<=(const GmcDual<N>& lhs, double rhs) { return lhs.value <= rhs; }
template <std::size_t N>
inline bool operator>=(const GmcDual<N>& lhs, double rhs) { return lhs.value >= rhs; }
template <std::size_t N>
inline bool operator==(const GmcDual<N>& lhs, double rhs) { return lhs.value == rhs; }
template <std::size_t N>
inline bool operator!=(const GmcDual<N>& lhs, double rhs) { return lhs.value != rhs; }
template <std::size_t N>
inline bool operator<(double lhs, const GmcDual<N>& rhs) { return lhs < rhs.value; }
template <std::size_t N>
inline bool operator>(double lhs, const GmcDual<N>& rhs) { return lhs > rhs.value; }
template <std::size_t N>
inline bool operator<=(double lhs, const GmcDual<N>& rhs) { return lhs <= rhs.value; }
template <std::size_t N>
inline bool operator>=(double lhs, const GmcDual<N>& rhs) { return lhs >= rhs.value; }
template <std::size_t N>
inline bool operator==(double lhs, const GmcDual<N>& rhs) { return lhs == rhs.value; }
template <std::size_t N>
inline bool operator!=(double lhs, const GmcDual<N>& rhs) { return lhs != rhs.value; }

template <std::size_t N>
inline GmcDual<N> gmcExp(const GmcDual<N>& input) {
    GmcDual<N> result;
    result.value = std::exp(input.value);
    for (std::size_t i = 0; i < N; ++i) result.derivative[i] = result.value * input.derivative[i];
    return result;
}

template <std::size_t N>
inline GmcDual<N> gmcLog(const GmcDual<N>& input) {
    GmcDual<N> result;
    result.value = std::log(input.value);
    for (std::size_t i = 0; i < N; ++i) result.derivative[i] = input.derivative[i] / input.value;
    return result;
}

template <std::size_t N>
inline GmcDual<N> gmcSqrt(const GmcDual<N>& input) {
    GmcDual<N> result;
    result.value = std::sqrt(input.value);
    const double scale = 0.5 / result.value;
    for (std::size_t i = 0; i < N; ++i) result.derivative[i] = scale * input.derivative[i];
    return result;
}

template <std::size_t N>
inline GmcDual<N> gmcPow(const GmcDual<N>& input, double exponent) {
    GmcDual<N> result;
    result.value = std::pow(input.value, exponent);
    const double scale = exponent * std::pow(input.value, exponent - 1.0);
    for (std::size_t i = 0; i < N; ++i) {
        result.derivative[i] = scale * input.derivative[i];
    }
    return result;
}

template <std::size_t N>
inline GmcDual<N> gmcTanh(const GmcDual<N>& input) {
    GmcDual<N> result;
    result.value = std::tanh(input.value);
    const double scale = 1.0 - result.value * result.value;
    for (std::size_t i = 0; i < N; ++i) result.derivative[i] = scale * input.derivative[i];
    return result;
}

template <std::size_t N>
inline GmcDual<N> gmcSin(const GmcDual<N>& input) {
    GmcDual<N> result;
    result.value = std::sin(input.value);
    const double scale = std::cos(input.value);
    for (std::size_t i = 0; i < N; ++i) result.derivative[i] = scale * input.derivative[i];
    return result;
}

template <std::size_t N>
inline GmcDual<N> gmcCos(const GmcDual<N>& input) {
    GmcDual<N> result;
    result.value = std::cos(input.value);
    const double scale = -std::sin(input.value);
    for (std::size_t i = 0; i < N; ++i) result.derivative[i] = scale * input.derivative[i];
    return result;
}

template <std::size_t N>
inline GmcDual<N> gmcTan(const GmcDual<N>& input) {
    GmcDual<N> result;
    result.value = std::tan(input.value);
    const double scale = 1.0 + result.value * result.value;
    for (std::size_t i = 0; i < N; ++i) result.derivative[i] = scale * input.derivative[i];
    return result;
}

template <std::size_t N>
inline GmcDual<N> gmcAtan(const GmcDual<N>& input) {
    GmcDual<N> result;
    result.value = std::atan(input.value);
    const double scale = 1.0 / (1.0 + input.value * input.value);
    for (std::size_t i = 0; i < N; ++i) result.derivative[i] = scale * input.derivative[i];
    return result;
}

template <std::size_t N>
inline GmcDual<N> gmcFmod(const GmcDual<N>& a, const GmcDual<N>& b) {
    // Between discontinuities d/dx fmod(x, b) == 1 w.r.t. x; the b-derivative
    // is dropped (b is a constant in the emitted uses, e.g. nf % 2 parity).
    GmcDual<N> result;
    result.value = std::fmod(a.value, b.value);
    for (std::size_t i = 0; i < N; ++i) result.derivative[i] = a.derivative[i];
    return result;
}

template <std::size_t N>
inline GmcDual<N> gmcFmod(const GmcDual<N>& a, double b) {
    GmcDual<N> result;
    result.value = std::fmod(a.value, b);
    for (std::size_t i = 0; i < N; ++i) result.derivative[i] = a.derivative[i];
    return result;
}

template <std::size_t N>
inline GmcDual<N> gmcSinh(const GmcDual<N>& input) {
    GmcDual<N> result;
    result.value = std::sinh(input.value);
    const double scale = std::cosh(input.value);
    for (std::size_t i = 0; i < N; ++i) result.derivative[i] = scale * input.derivative[i];
    return result;
}

template <std::size_t N>
inline GmcDual<N> gmcCosh(const GmcDual<N>& input) {
    GmcDual<N> result;
    result.value = std::cosh(input.value);
    const double scale = std::sinh(input.value);
    for (std::size_t i = 0; i < N; ++i) result.derivative[i] = scale * input.derivative[i];
    return result;
}

template <std::size_t N>
inline GmcDual<N> gmcAbs(const GmcDual<N>& input) {
    return input.value < 0.0 ? -input : input;
}

template <std::size_t N>
inline GmcDual<N> gmcMin(const GmcDual<N>& lhs, const GmcDual<N>& rhs) {
    return lhs.value <= rhs.value ? lhs : rhs;
}

template <std::size_t N>
inline GmcDual<N> gmcMax(const GmcDual<N>& lhs, const GmcDual<N>& rhs) {
    return lhs.value >= rhs.value ? lhs : rhs;
}

template <std::size_t N>
inline GmcDual<N> gmcMin(const GmcDual<N>& lhs, double rhs) {
    return lhs.value <= rhs ? lhs : GmcDual<N>::constant(rhs);
}

template <std::size_t N>
inline GmcDual<N> gmcMax(const GmcDual<N>& lhs, double rhs) {
    return lhs.value >= rhs ? lhs : GmcDual<N>::constant(rhs);
}

// Piecewise-linear limiter with zero derivative in the clamped region, used
// by compact-model voltage/current limiting.
template <std::size_t N>
inline GmcDual<N> gmcLim(const GmcDual<N>& input, double low, double high) {
    if (input.value <= low) return GmcDual<N>::constant(low);
    if (input.value >= high) return GmcDual<N>::constant(high);
    return input;
}

// Compatibility alias used by the PSP103 core and the GMC C++ emitter.
using GmcDual4 = GmcDual<4>;

// Sizing used by BSIM-class four-terminal models with internal nodes.
using GmcDual8 = GmcDual<8>;

} // namespace gspice

#endif // GSPICE_GMC_DUAL_HPP
