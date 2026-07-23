#ifndef GSPICE_FOURIER_HPP
#define GSPICE_FOURIER_HPP

// Ensure M_PI is available on MSVC (must be before any math includes).
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace gspice {

class Fourier {
public:
    // -----------------------------------------------------------------------
    // forward()
    //   Converts N real time-domain samples to frequency-domain coefficients.
    //   Returns (N/2 + 1) complex bins: bin 0 = DC, bin k = amplitude at
    //   k·f₀, where f₀ = 1/T. Bins are scaled so that bin k carries the
    //   two-sided amplitude of the k-th harmonic (same convention as the
    //   original GSPICE DFT implementation).
    // -----------------------------------------------------------------------
    static std::vector<std::complex<double>>
    forward(const std::vector<double>& time_samples) {
        const int N_in = static_cast<int>(time_samples.size());
        if (N_in <= 0) return {};
        const int N = nextPow2(N_in);        // pad to power of 2
        const int K = N / 2;

        // Promote real input to complex, zero-pad.
        std::vector<std::complex<double>> buf(static_cast<std::size_t>(N));
        for (int i = 0; i < N_in; ++i) {
            buf[static_cast<std::size_t>(i)] = {time_samples[static_cast<std::size_t>(i)], 0.0};
        }

        fftInPlace(buf, /*inverse=*/false);

        // Extract positive-frequency bins and apply the same normalisation
        // as the original DFT: divide by N, multiply harmonics (k>0) by 2.
        const int bins = K + 1;
        std::vector<std::complex<double>> harmonics(static_cast<std::size_t>(bins));
        for (int k = 0; k < bins; ++k) {
            harmonics[static_cast<std::size_t>(k)] =
                buf[static_cast<std::size_t>(k)] / static_cast<double>(N);
            if (k > 0) harmonics[static_cast<std::size_t>(k)] *= 2.0;
        }
        return harmonics;
    }

    // -----------------------------------------------------------------------
    // inverse()
    //   Converts K+1 complex harmonics back to N real time-domain samples.
    //   The input is the same format as returned by forward().
    // -----------------------------------------------------------------------
    static std::vector<double>
    inverse(const std::vector<std::complex<double>>& harmonics, int N) {
        if (N <= 0 || harmonics.empty()) return std::vector<double>(static_cast<std::size_t>(N), 0.0);

        const int N_fft = nextPow2(N);
        const int K = static_cast<int>(harmonics.size()) - 1;

        // Reconstruct the full two-sided spectrum.
        std::vector<std::complex<double>> buf(static_cast<std::size_t>(N_fft), {0.0, 0.0});
        // DC
        buf[0] = harmonics[0] * static_cast<double>(N_fft);
        // Positive frequencies (k > 0): each bin carries 2x amplitude so
        // restore the one-sided value (divide by 2), then scale by N_fft.
        for (int k = 1; k <= std::min(K, N_fft / 2 - 1); ++k) {
            const std::complex<double> val =
                harmonics[static_cast<std::size_t>(k)] * (static_cast<double>(N_fft) / 2.0);
            buf[static_cast<std::size_t>(k)]        = val;
            buf[static_cast<std::size_t>(N_fft - k)] = std::conj(val);
        }
        if (N_fft / 2 <= K) {
            buf[static_cast<std::size_t>(N_fft / 2)] =
                harmonics[static_cast<std::size_t>(N_fft / 2)] *
                (static_cast<double>(N_fft) / 2.0);
        }

        fftInPlace(buf, /*inverse=*/true);

        // Return only the first N (original) samples; scale by 1/N_fft is
        // already baked into fftInPlace when inverse=true.
        std::vector<double> out(static_cast<std::size_t>(N));
        for (int i = 0; i < N; ++i) {
            out[static_cast<std::size_t>(i)] = buf[static_cast<std::size_t>(i)].real();
        }
        return out;
    }

    // -----------------------------------------------------------------------
    // forwardFull()
    //   Transforms a complex vector in-place to frequency domain.
    //   Exposed for the Harmonic Balance engine which works entirely in
    //   the complex-valued frequency domain.
    //   N must be a power of 2.
    // -----------------------------------------------------------------------
    static void forwardFull(std::vector<std::complex<double>>& buf) {
        fftInPlace(buf, /*inverse=*/false);
    }

    // -----------------------------------------------------------------------
    // inverseFull()
    //   Transforms a complex vector in-place to time domain (IFFT).
    //   Values are divided by N (normalised IDFT convention).
    //   N must be a power of 2.
    // -----------------------------------------------------------------------
    static void inverseFull(std::vector<std::complex<double>>& buf) {
        fftInPlace(buf, /*inverse=*/true);
    }

    // -----------------------------------------------------------------------
    // Utility: smallest power of 2 >= n.
    // -----------------------------------------------------------------------
    static int nextPow2(int n) {
        if (n <= 1) return 1;
        int p = 1;
        while (p < n) p <<= 1;
        return p;
    }

private:
    // -----------------------------------------------------------------------
    // Iterative Cooley-Tukey radix-2 decimation-in-time FFT.
    //
    // Operates in-place on buf[0..N-1] where N = buf.size() must be a power
    // of two. When inverse=false the transform is the unnormalised forward DFT.
    // When inverse=true the output is divided by N (normalised IDFT).
    // -----------------------------------------------------------------------
    static void fftInPlace(std::vector<std::complex<double>>& buf, bool inverse) {
        const int N = static_cast<int>(buf.size());
        if (N <= 1) return;

        // ---- Bit-reversal permutation -----------------------------------------
        // Swap element i with bit_reverse(i, log2(N)) to permute into
        // decimation-in-time order without recursion.
        const int bits = log2Int(N);
        for (int i = 0; i < N; ++i) {
            const int j = bitReverse(i, bits);
            if (j > i) {
                std::swap(buf[static_cast<std::size_t>(i)],
                          buf[static_cast<std::size_t>(j)]);
            }
        }

        // ---- Butterfly stages -----------------------------------------------
        // Stage s processes pairs of butterfly groups of size 2^s.
        // The twiddle factor for group of size m is W_m = exp(-j·2π/m)
        // (or +j for the inverse). Each stage is O(N).
        const double sign = inverse ? 1.0 : -1.0;
        for (int s = 1; s <= bits; ++s) {
            const int m   = 1 << s;          // group size
            const int mh  = m >> 1;          // half-group size
            // Principal twiddle factor for this stage.
            const std::complex<double> Wm{
                std::cos(2.0 * M_PI / static_cast<double>(m)),
                sign * std::sin(2.0 * M_PI / static_cast<double>(m))
            };
            for (int k = 0; k < N; k += m) {
                std::complex<double> W = {1.0, 0.0};
                for (int j = 0; j < mh; ++j) {
                    const auto u = buf[static_cast<std::size_t>(k + j)];
                    const auto v = buf[static_cast<std::size_t>(k + j + mh)] * W;
                    buf[static_cast<std::size_t>(k + j)]       = u + v;
                    buf[static_cast<std::size_t>(k + j + mh)]  = u - v;
                    W *= Wm;
                }
            }
        }

        // Normalise for inverse transform.
        if (inverse) {
            const double inv_N = 1.0 / static_cast<double>(N);
            for (auto& c : buf) c *= inv_N;
        }
    }

    // ---- Integer log₂ for powers of two ------------------------------------
    static int log2Int(int n) {
        int k = 0;
        while ((1 << k) < n) ++k;
        return k;
    }

    // ---- Bit-reversal of an integer with `bits` significant bits -----------
    static int bitReverse(int x, int bits) {
        int r = 0;
        for (int i = 0; i < bits; ++i) {
            r = (r << 1) | (x & 1);
            x >>= 1;
        }
        return r;
    }
};

} // namespace gspice

#endif // GSPICE_FOURIER_HPP
