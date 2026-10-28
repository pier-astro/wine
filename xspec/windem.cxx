/*
 * WINE Emission Convolution Model for XSPEC
 *
 * Model: windem
 * Convolves emission spectra with relativistic outflow profile using FFT.
 * Uses 2x oversampling and works in flux per d(log E) for correct Jacobian.
 *
 * Public coordinates are converted to physical angles once, before the
 * shared convolution. Cf is the occupied fraction of one hemisphere.
 */

#include <xsTypes.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <set>
#include <vector>

#include "wine_lineshape.h"
#include "wine_geometry.h"
#include "wine_xspec.h"

#if __has_include("fftw3.h")
#  include "fftw3.h"
#else
#  include "fftw/fftw3.h"
#endif


static void rebinConservative(const std::vector<double>& oldEdge,
                              const std::vector<double>& oldFlux,
                              const std::vector<double>& newEdge,
                              std::vector<double>& newFlux) {
    newFlux.assign(newEdge.size() - 1, 0.0);
    size_t i = 0;
    size_t j = 0;
    while (i < oldFlux.size() && j < newFlux.size()) {
        const double lo = std::max(oldEdge[i], newEdge[j]);
        const double hi = std::min(oldEdge[i + 1], newEdge[j + 1]);
        if (hi > lo) {
            newFlux[j] += oldFlux[i] * (hi - lo) / (oldEdge[i + 1] - oldEdge[i]);
        }
        if (oldEdge[i + 1] <= newEdge[j + 1]) {
            ++i;
        } else {
            ++j;
        }
    }
}

static int nextPowerOfTwo(int n) {
    int result = 1;
    while (result < n) result *= 2;
    return result;
}

static void convolveEmission(const RealArray& energy, double beta,
                             const wine::Geometry& geometry, int spectrum,
                             RealArray& flux, RealArray& fluxErr) {
    const int n_bins = static_cast<int>(energy.size()) - 1;
    fluxErr.resize(0);
    if (n_bins <= 0) {
        return;
    }
    if (geometry.Cf == 0.0) {
        flux = 0.0;
        return;
    }
    for (int i = 0; i < n_bins; ++i) {
        if (!(energy[i] > 0.0 && energy[i + 1] > energy[i])) {
            throw FunctionUtility::FunctionException(
                "windem: energy edges must be positive and strictly increasing.");
        }
    }
    if (std::abs(beta) < 1e-10) {
        for (int i = 0; i < n_bins; ++i) {
            flux[i] *= geometry.Cf / 2.0;
        }
        return;
    }

    static std::set<int> warnedSpectra;
    if (warnedSpectra.insert(spectrum).second) {
        FunctionUtility::xsWrite(
            "\n\t*** WINE windem ENERGY-GRID NOTICE ***\n"
            "\tConvolution uses only the input spectrum evaluated by XSPEC. "
            "Missing input energies can suppress emission near a grid boundary. "
            "Extend the model grid beyond the fitted band, within the input "
            "table's energy coverage, and check stability in the fitted band.\n",
            10);
    }

    std::vector<double> inputEdge(energy.size());
    for (int i = 0; i < n_bins; ++i) {
        inputEdge[i] = energy[i];
    }
    inputEdge[n_bins] = energy[n_bins];

    const int oversample = 2;
    const int n_os = n_bins * oversample;
    if (n_os < 2) {
        return;
    }

    const double logE_min = std::log(inputEdge.front());
    const double dlog = (std::log(inputEdge.back()) - logE_min) / n_os;
    std::vector<double> uniformEdge(n_os + 1);
    for (int i = 0; i <= n_os; ++i) {
        uniformEdge[i] = std::exp(logE_min + dlog * i);
    }
    std::vector<double> inputFlux(n_bins);
    for (int i = 0; i < n_bins; ++i) inputFlux[i] = flux[i];
    std::vector<double> F_uniform;
    rebinConservative(inputEdge, inputFlux, uniformEdge, F_uniform);

    std::vector<double> E_kernel(n_os);
    std::vector<double> kernel_raw(n_os);
    std::vector<double> kernel_logE(n_os);
    const int mid = n_os / 2;
    for (int i = 0; i < n_os; ++i) {
        int k = i - mid;
        E_kernel[i] = std::exp(static_cast<double>(k) * dlog);
    }

    wine::lineshape(E_kernel.data(), n_os, 1.0, beta, geometry.incl,
                    geometry.tin, geometry.tout, kernel_raw.data());
    for (int i = 0; i < n_os; ++i) {
        kernel_logE[i] = kernel_raw[i] * E_kernel[i];
    }

    const double fraction = geometry.Cf / 2.0;
    double kernelSum = 0.0;
    for (int i = 0; i < n_os; ++i) {
        kernel_logE[i] *= dlog;
        kernelSum += kernel_logE[i];
    }
    if (kernelSum > 0.0) {
        for (double& value : kernel_logE) value *= fraction / kernelSum;
    }

    const int n_fft = nextPowerOfTwo(3 * n_os - 2);
    fftw_complex* F_fft = fftw_alloc_complex(n_fft);
    fftw_complex* K_fft = fftw_alloc_complex(n_fft);
    fftw_complex* G_fft = fftw_alloc_complex(n_fft);

    for (int i = 0; i < n_fft; ++i) {
        F_fft[i][0] = i < n_os ? F_uniform[i] : 0.0;
        F_fft[i][1] = 0.0;
        K_fft[i][0] = 0.0;
        K_fft[i][1] = 0.0;
    }
    for (int i = 0; i < n_os; ++i) {
        const int lag = i - mid;
        K_fft[lag >= 0 ? lag : n_fft + lag][0] = kernel_logE[i];
    }

    fftw_plan plan_F = fftw_plan_dft_1d(n_fft, F_fft, F_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_K = fftw_plan_dft_1d(n_fft, K_fft, K_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_G = fftw_plan_dft_1d(n_fft, G_fft, G_fft, FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(plan_F);
    fftw_execute(plan_K);

    for (int i = 0; i < n_fft; ++i) {
        double a = F_fft[i][0];
        double b = F_fft[i][1];
        double c = K_fft[i][0];
        double d = K_fft[i][1];
        G_fft[i][0] = a * c - b * d;
        G_fft[i][1] = a * d + b * c;
    }

    fftw_execute(plan_G);

    std::vector<double> result_uniform(n_os);
    const double inv_n = 1.0 / static_cast<double>(n_fft);
    for (int i = 0; i < n_os; ++i) {
        result_uniform[i] = G_fft[i][0] * inv_n;
    }

    fftw_destroy_plan(plan_F);
    fftw_destroy_plan(plan_K);
    fftw_destroy_plan(plan_G);
    fftw_free(F_fft);
    fftw_free(K_fft);
    fftw_free(G_fft);

    std::vector<double> result;
    rebinConservative(uniformEdge, result_uniform, inputEdge, result);
    for (int i = 0; i < n_bins; ++i) flux[i] = result[i];
}

extern "C" void windem(const RealArray& energy, const RealArray& params,
                        int spectrum, RealArray& flux, RealArray& fluxErr,
                        const string& init) {
    const double beta = wine::boundedParameter(params[0], 0.0, 0.99, "vout", spectrum);
    const double Cf = wine::boundedParameter(params[1], 0.0, 1.0, "Cf", spectrum);
    const double hole = wine::boundedParameter(params[2], 0.0, 1.0, "hole", spectrum);
    const double incl = wine::boundedParameter(params[3], 0.0, 90.0, "incl", spectrum);
    const wine::Geometry geometry = wine::freeGeometry(Cf, hole, incl * wine::kDeg2Rad);
    convolveEmission(energy, beta, geometry, spectrum, flux, fluxErr);
}

extern "C" void windemlos(const RealArray& energy, const RealArray& params,
                           int spectrum, RealArray& flux, RealArray& fluxErr,
                           const string& init) {
    const double beta = wine::boundedParameter(params[0], 0.0, 0.99, "vout", spectrum);
    const double Cf = wine::boundedParameter(params[1], 0.0, 1.0, "Cf", spectrum);
    const double skew = wine::boundedParameter(params[2], 0.0, 1.0, "skew", spectrum);
    const double incl = wine::boundedParameter(params[3], 0.0, 90.0, "incl", spectrum);
    const wine::Geometry geometry = wine::losGeometry(Cf, skew, incl * wine::kDeg2Rad);
    convolveEmission(energy, beta, geometry, spectrum, flux, fluxErr);
}

extern "C" void windemoff(const RealArray& energy, const RealArray& params,
                           int spectrum, RealArray& flux, RealArray& fluxErr,
                           const string& init) {
    const double beta = wine::boundedParameter(params[0], 0.0, 0.99, "vout", spectrum);
    const double Cf = wine::boundedParameter(params[1], 0.0, 1.0, "Cf", spectrum);
    const double hole = wine::boundedParameter(params[2], 0.0, 1.0, "hole", spectrum);
    const double view = wine::boundedParameter(params[3], 0.0, 1.0, "view", spectrum);
    if (params[4] != 0.0 && params[4] != 1.0) {
        throw FunctionUtility::FunctionException(
            "windemoff: region must be 0 (polar) or 1 (equatorial).");
    }
    const wine::Geometry geometry = wine::offGeometry(Cf, hole, view, static_cast<int>(params[4]));
    convolveEmission(energy, beta, geometry, spectrum, flux, fluxErr);
}
