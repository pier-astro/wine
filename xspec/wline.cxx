/* Single relativistic line: public coordinates share the same physical
 * geometry conversion as the emission convolution. Norm is applied by XSPEC.
 */
#include <xsTypes.h>
#include <cmath>
#include <vector>

#include "redwine_lineshape.h"
#include "redwine_geometry.h"
#include "redwine_xspec.h"

static void lineEmission(const RealArray& energy, double beta,
                         const redwine::Geometry& geometry, double lineE,
                         RealArray& flux, RealArray& fluxErr) {
    const int n_bins = static_cast<int>(energy.size()) - 1;
    flux.resize(std::max(0, n_bins));
    flux = 0.0;
    fluxErr.resize(0);

    if (n_bins <= 0 || lineE <= 0.0 || geometry.Cf == 0.0) {
        return;
    }

    if (std::abs(beta) < 1e-10) {
        const double lineflux = geometry.Cf / 2.0;

        int idx = -1;
        for (int i = 0; i < n_bins; ++i) {
            if (lineE >= energy[i] && lineE < energy[i + 1]) {
                idx = i;
                break;
            }
        }
        if (idx >= 0) {
            flux[idx] = lineflux;
        }
        return;
    }

    std::vector<double> E_center(n_bins);
    std::vector<double> dE(n_bins);
    for (int i = 0; i < n_bins; ++i) {
        E_center[i] = std::sqrt(energy[i] * energy[i + 1]);
        dE[i] = energy[i + 1] - energy[i];
    }

    std::vector<double> profile(n_bins);
    redwine::lineshape(E_center.data(), n_bins, lineE, beta, geometry.incl, geometry.tin, geometry.tout, profile.data());

    for (int i = 0; i < n_bins; ++i) {
        // XSPEC applies the additive normalization, keep internal norm at 1.
        flux[i] = profile[i] * dE[i];
    }
}

extern "C" void wlineCf(const RealArray& energy, const RealArray& params,
                         int spectrum, RealArray& flux, RealArray& fluxErr,
                         const string& init) {
    const double beta = redwine::boundedParameter(params[0], -0.99, 0.99, "vout", spectrum);
    const auto geometry = redwine::freeCf(params, spectrum);
    lineEmission(energy, beta, geometry, params[4], flux, fluxErr);
}

extern "C" void wlinelosCf(const RealArray& energy, const RealArray& params,
                         int spectrum, RealArray& flux, RealArray& fluxErr,
                         const string& init) {
    const double beta = redwine::boundedParameter(params[0], -0.99, 0.99, "vout", spectrum);
    const auto geometry = redwine::losCf(params, spectrum);
    lineEmission(energy, beta, geometry, params[4], flux, fluxErr);
}

extern "C" void wlineoffCf(const RealArray& energy, const RealArray& params,
                         int spectrum, RealArray& flux, RealArray& fluxErr,
                         const string& init) {
    const double beta = redwine::boundedParameter(params[0], -0.99, 0.99, "vout", spectrum);
    const auto geometry = redwine::offCf(params, spectrum);
    lineEmission(energy, beta, geometry, params[5], flux, fluxErr);
}

extern "C" void wline(const RealArray& energy, const RealArray& params,
                         int spectrum, RealArray& flux, RealArray& fluxErr,
                         const string& init) {
    const double beta = redwine::boundedParameter(params[0], -0.99, 0.99, "vout", spectrum);
    const auto geometry = redwine::freeangles(params, spectrum);
    lineEmission(energy, beta, geometry, params[4], flux, fluxErr);
}

extern "C" void wlinelos(const RealArray& energy, const RealArray& params,
                         int spectrum, RealArray& flux, RealArray& fluxErr,
                         const string& init) {
    const double beta = redwine::boundedParameter(params[0], -0.99, 0.99, "vout", spectrum);
    const auto geometry = redwine::losangles(params, spectrum);
    lineEmission(energy, beta, geometry, params[4], flux, fluxErr);
}

extern "C" void wlineoff(const RealArray& energy, const RealArray& params,
                         int spectrum, RealArray& flux, RealArray& fluxErr,
                         const string& init) {
    const double beta = redwine::boundedParameter(params[0], -0.99, 0.99, "vout", spectrum);
    const auto geometry = redwine::offangles(params, spectrum);
    lineEmission(energy, beta, geometry, params[5], flux, fluxErr);
}
