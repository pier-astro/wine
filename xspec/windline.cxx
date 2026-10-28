/* Single relativistic line: public coordinates share the same physical
 * geometry conversion as the emission convolution. Norm is applied by XSPEC.
 */
#include <xsTypes.h>
#include <cmath>
#include <vector>

#include "wine_lineshape.h"
#include "wine_geometry.h"
#include "wine_xspec.h"

static void lineEmission(const RealArray& energy, double beta,
                         const wine::Geometry& geometry, double lineE,
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
    wine::lineshape(E_center.data(), n_bins, lineE, beta, geometry.incl, geometry.tin, geometry.tout, profile.data());

    for (int i = 0; i < n_bins; ++i) {
        // XSPEC applies the additive normalization, keep internal norm at 1.
        flux[i] = profile[i] * dE[i];
    }
}

extern "C" void windline(const RealArray& energy, const RealArray& params,
                         int spectrum, RealArray& flux, RealArray& fluxErr,
                         const string& init) {
    const double beta = wine::boundedParameter(params[0], 0.0, 0.99, "vout", spectrum);
    const double Cf = wine::boundedParameter(params[1], 0.0, 1.0, "Cf", spectrum);
    const double hole = wine::boundedParameter(params[2], 0.0, 1.0, "hole", spectrum);
    const double incl = wine::boundedParameter(params[3], 0.0, 90.0, "incl", spectrum);
    const wine::Geometry geometry = wine::freeGeometry(Cf, hole, incl * wine::kDeg2Rad);
    lineEmission(energy, beta, geometry, params[4], flux, fluxErr);
}

extern "C" void windlinelos(const RealArray& energy, const RealArray& params,
                         int spectrum, RealArray& flux, RealArray& fluxErr,
                         const string& init) {
    const double beta = wine::boundedParameter(params[0], 0.0, 0.99, "vout", spectrum);
    const double Cf = wine::boundedParameter(params[1], 0.0, 1.0, "Cf", spectrum);
    const double skew = wine::boundedParameter(params[2], 0.0, 1.0, "skew", spectrum);
    const double incl = wine::boundedParameter(params[3], 0.0, 90.0, "incl", spectrum);
    const wine::Geometry geometry = wine::losGeometry(Cf, skew, incl * wine::kDeg2Rad);
    lineEmission(energy, beta, geometry, params[4], flux, fluxErr);
}

extern "C" void windlineoff(const RealArray& energy, const RealArray& params,
                         int spectrum, RealArray& flux, RealArray& fluxErr,
                         const string& init) {
    const double beta = wine::boundedParameter(params[0], 0.0, 0.99, "vout", spectrum);
    const double Cf = wine::boundedParameter(params[1], 0.0, 1.0, "Cf", spectrum);
    const double hole = wine::boundedParameter(params[2], 0.0, 1.0, "hole", spectrum);
    const double view = wine::boundedParameter(params[3], 0.0, 1.0, "view", spectrum);
    if (params[4] != 0.0 && params[4] != 1.0) {
        throw FunctionUtility::FunctionException(
            "windlineoff: region must be 0 (polar) or 1 (equatorial).");
    }
    const wine::Geometry geometry = wine::offGeometry(Cf, hole, view, static_cast<int>(params[4]));
    lineEmission(energy, beta, geometry, params[5], flux, fluxErr);
}
