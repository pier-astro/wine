/*
 * WINE Emission Convolution Model for XSPEC
 *
 * Model: windem
 * Convolves emission spectra with relativistic outflow profile using FFT.
 * Uses 2x oversampling and works in flux per d(log E) for correct Jacobian.
 *
 * Parameters:
 *   vout  - Outflow velocity in units of c (beta = v/c)
 *   incl  - Inclination angle in degrees
 *   tout  - Outer opening angle in degrees
 *   tin   - Inner opening angle in degrees
 *
 * windemlos replaces incl with uincl in [0, 1] and evaluates
 * incl = tin + uincl * (tout - tin), forcing the line of sight through the wind.
 *
 * windemoff replaces incl with uview in [0, 1] and a fixed region switch:
 *   region = 0: incl = (1 - uview) * tin       (polar cavity)
 *   region = 1: incl = tout + uview * (90 - tout) (outside the wind cone)
 */

#include <xsTypes.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <set>
#include <vector>

#include "wine_lineshape.h"

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

extern "C" void windem(const RealArray& energy, const RealArray& params,
                       int spectrum, RealArray& flux, RealArray& fluxErr,
                       const string& init) {
    double beta = params[0];
    double incl_deg = params[1];
    double tout_deg = params[2];
    double tin_deg = params[3];

    if (!std::isfinite(beta) || beta < 0.0 || beta >= 1.0 ||
        tout_deg <= tin_deg) {
        throw FunctionUtility::FunctionException(
            "windem: require 0 <= vout < 1 and tout > tin.");
    }

    double incl = incl_deg * wine::kDeg2Rad;
    double theta_out = tout_deg * wine::kDeg2Rad;
    double theta_in = tin_deg * wine::kDeg2Rad;

    const int n_bins = static_cast<int>(energy.size()) - 1;
    if (n_bins <= 0) {
        return;
    }
    for (int i = 0; i < n_bins; ++i) {
        if (!(energy[i] > 0.0 && energy[i + 1] > energy[i])) {
            throw FunctionUtility::FunctionException(
                "windem: energy edges must be positive and strictly increasing.");
        }
    }
    fluxErr.resize(0);

    if (std::abs(beta) < 1e-10) {
        double solid_angle = 2.0 * wine::kPi * (std::cos(theta_in) - std::cos(theta_out));
        double fraction = solid_angle / (4.0 * wine::kPi);
        for (int i = 0; i < n_bins; ++i) {
            flux[i] *= fraction;
        }
        return;
    }

    static std::set<int> warnedSpectra;
    if (warnedSpectra.insert(spectrum).second) {
        FunctionUtility::xsWrite(
            "\n\t*** WINE windem ENERGY-GRID NOTICE ***\n"
            "\tEmission shifted beyond a finite model grid is physically lost. "
            "Extend both energy boundaries beyond the fitted detector band.\n",
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

    wine::lineshape(E_kernel.data(), n_os, 1.0, beta, incl, theta_in, theta_out, kernel_raw.data());
    for (int i = 0; i < n_os; ++i) {
        kernel_logE[i] = kernel_raw[i] * E_kernel[i];
    }

    const double fraction =
        (std::cos(theta_in) - std::cos(theta_out)) / 2.0;
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

extern "C" void windemlos(const RealArray& energy, const RealArray& params,
                           int spectrum, RealArray& flux, RealArray& fluxErr,
                           const string& init) {
    const double uincl = params[1];
    const double tout_deg = params[2];
    const double tin_deg = params[3];
    // XSPEC keeps fitted uincl in [0, 1], but Levenberg evaluates finite
    // differences beyond hard limits. Evaluate those probes with windem.
    if (!std::isfinite(uincl) ||
        !std::isfinite(tout_deg) || !std::isfinite(tin_deg) ||
        tout_deg <= tin_deg) {
        throw FunctionUtility::FunctionException(
            "windemlos: require finite uincl and tout > tin.");
    }

    RealArray physicalParams(params);
    physicalParams[1] = tin_deg + uincl * (tout_deg - tin_deg);
    windem(energy, physicalParams, spectrum, flux, fluxErr, init);
}

extern "C" void windemoff(const RealArray& energy, const RealArray& params,
                           int spectrum, RealArray& flux, RealArray& fluxErr,
                           const string& init) {
    const double uview = params[1];
    const double tout_deg = params[2];
    const double tin_deg = params[3];
    const double region = params[4];
    // Keep the same numerical extension for finite-difference probes.
    if (!std::isfinite(uview) ||
        !std::isfinite(tout_deg) || !std::isfinite(tin_deg) ||
        tin_deg < 0.0 || tout_deg > 90.0 || tout_deg <= tin_deg ||
        (region != 0.0 && region != 1.0)) {
        throw FunctionUtility::FunctionException(
            "windemoff: require finite uview, "
            "0 <= tin < tout <= 90 deg, and region = 0 or 1.");
    }

    RealArray physicalParams(params);
    physicalParams[1] = region == 0.0
        ? (1.0 - uview) * tin_deg
        : tout_deg + uview * (90.0 - tout_deg);
    windem(energy, physicalParams, spectrum, flux, fluxErr, init);
}
