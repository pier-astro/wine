/*
 * RED-WINE Emission Convolution Model for XSPEC
 *
 * Model: wem_Cf
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

#include "redwine_lineshape.h"
#include "redwine_geometry.h"
#include "redwine_xspec.h"

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
                             const redwine::Geometry& geometry, int spectrum,
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
                "wem_Cf: energy edges must be positive and strictly increasing.");
        }
    }
    const int oversample = 2;
    const int n_os = n_bins * oversample;
    const double logE_min = std::log(energy[0]);
    const double dlog = (std::log(energy[n_bins]) - logE_min) / n_os;
    // When even the largest Doppler shift lies inside half a kernel bin,
    // the resolved profile is a delta at zero shift.
    if (std::atanh(std::abs(beta)) < 0.5 * dlog) {
        for (int i = 0; i < n_bins; ++i) {
            flux[i] *= geometry.Cf / 2.0;
        }
        return;
    }

    static std::set<int> warnedSpectra;
    if (warnedSpectra.insert(spectrum).second) {
        FunctionUtility::xsWrite(
            "\n\t*** RED-WINE wem ENERGY-GRID NOTICE ***\n"
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

    redwine::lineshape(E_kernel.data(), n_os, 1.0, beta, geometry.incl,
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

extern "C" void wemCf(const RealArray& energy, const RealArray& params,
                         int spectrum, RealArray& flux, RealArray& fluxErr,
                         const string& init) {
    const double beta = redwine::boundedParameter(params[0], -0.99, 0.99, "vout", spectrum);
    const auto geometry = redwine::freeCf(params, spectrum);
    convolveEmission(energy, beta, geometry, spectrum, flux, fluxErr);
}

extern "C" void wemlosCf(const RealArray& energy, const RealArray& params,
                         int spectrum, RealArray& flux, RealArray& fluxErr,
                         const string& init) {
    const double beta = redwine::boundedParameter(params[0], -0.99, 0.99, "vout", spectrum);
    const auto geometry = redwine::losCf(params, spectrum);
    convolveEmission(energy, beta, geometry, spectrum, flux, fluxErr);
}

extern "C" void wemoffCf(const RealArray& energy, const RealArray& params,
                         int spectrum, RealArray& flux, RealArray& fluxErr,
                         const string& init) {
    const double beta = redwine::boundedParameter(params[0], -0.99, 0.99, "vout", spectrum);
    const auto geometry = redwine::offCf(params, spectrum);
    convolveEmission(energy, beta, geometry, spectrum, flux, fluxErr);
}

extern "C" void wem(const RealArray& energy, const RealArray& params,
                         int spectrum, RealArray& flux, RealArray& fluxErr,
                         const string& init) {
    const double beta = redwine::boundedParameter(params[0], -0.99, 0.99, "vout", spectrum);
    const auto geometry = redwine::freeangles(params, spectrum);
    convolveEmission(energy, beta, geometry, spectrum, flux, fluxErr);
}

extern "C" void wemlos(const RealArray& energy, const RealArray& params,
                         int spectrum, RealArray& flux, RealArray& fluxErr,
                         const string& init) {
    const double beta = redwine::boundedParameter(params[0], -0.99, 0.99, "vout", spectrum);
    const auto geometry = redwine::losangles(params, spectrum);
    convolveEmission(energy, beta, geometry, spectrum, flux, fluxErr);
}

extern "C" void wemoff(const RealArray& energy, const RealArray& params,
                         int spectrum, RealArray& flux, RealArray& fluxErr,
                         const string& init) {
    const double beta = redwine::boundedParameter(params[0], -0.99, 0.99, "vout", spectrum);
    const auto geometry = redwine::offangles(params, spectrum);
    convolveEmission(energy, beta, geometry, spectrum, flux, fluxErr);
}
