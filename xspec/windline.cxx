/*
 * WINE Bare Emission Line Model for XSPEC
 *
 * Model: windline
 * Computes the relativistic line profile for a single emission line.
 *
 * Parameters:
 *   vout   - Outflow velocity in units of c (beta = v/c)
 *   incl   - Inclination angle in degrees
 *   tout   - Outer opening angle in degrees
 *   tin    - Inner opening angle in degrees
 *   LineE  - Rest-frame line energy in keV
 *   norm   - XSPEC additive normalization (applied externally)
 *
 * windlinelos replaces incl with uincl in [0, 1] and evaluates
 * incl = tin + uincl * (tout - tin), forcing the line of sight through the wind.
 *
 * windlineoff replaces incl with uview in [0, 1] and a fixed region switch:
 *   region = 0: incl = (1 - uview) * tin       (polar cavity)
 *   region = 1: incl = tout + uview * (90 - tout) (outside the wind cone)
 */

#include <xsTypes.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <cmath>
#include <vector>

#include "wine_lineshape.h"


extern "C" void windline(const RealArray& energy, const RealArray& params,
                         int spectrum, RealArray& flux, RealArray& fluxErr,
                         const string& init) {
    const double beta = params[0];
    const double incl_deg = params[1];
    const double tout_deg = params[2];
    const double tin_deg = params[3];
    const double lineE = params[4];

    if (tout_deg <= tin_deg) {
        return;
    }

    const int n_bins = static_cast<int>(energy.size()) - 1;
    if (n_bins <= 0 || lineE <= 0.0) {
        return;
    }

    flux.resize(n_bins);
    flux = 0.0;
    fluxErr.resize(0);

    const double incl = incl_deg * wine::kDeg2Rad;
    const double theta_out = tout_deg * wine::kDeg2Rad;
    const double theta_in = tin_deg * wine::kDeg2Rad;

    if (std::abs(beta) < 1e-10) {
        const double solid_angle = 2.0 * wine::kPi * (std::cos(theta_in) - std::cos(theta_out));
        const double fraction = solid_angle / (4.0 * wine::kPi);
        const double lineflux = fraction;

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
    wine::lineshape(E_center.data(), n_bins, lineE, beta, incl, theta_in, theta_out, profile.data());

    for (int i = 0; i < n_bins; ++i) {
        // XSPEC applies the additive normalization, keep internal norm at 1.
        flux[i] = profile[i] * dE[i];
    }
}

extern "C" void windlinelos(const RealArray& energy, const RealArray& params,
                             int spectrum, RealArray& flux, RealArray& fluxErr,
                             const string& init) {
    const double uincl = params[1];
    const double tout_deg = params[2];
    const double tin_deg = params[3];
    // The XSPEC parameter bounds constrain fits; derivative probes can leave
    // that range, so evaluate them with the underlying windline model.
    if (!std::isfinite(uincl) ||
        !std::isfinite(tout_deg) || !std::isfinite(tin_deg) ||
        tout_deg <= tin_deg) {
        throw FunctionUtility::FunctionException(
            "windlinelos: require finite uincl and tout > tin.");
    }

    RealArray physicalParams(params);
    physicalParams[1] = tin_deg + uincl * (tout_deg - tin_deg);
    windline(energy, physicalParams, spectrum, flux, fluxErr, init);
}

extern "C" void windlineoff(const RealArray& energy, const RealArray& params,
                             int spectrum, RealArray& flux, RealArray& fluxErr,
                             const string& init) {
    const double uview = params[1];
    const double tout_deg = params[2];
    const double tin_deg = params[3];
    const double region = params[4];
    if (!std::isfinite(uview) ||
        !std::isfinite(tout_deg) || !std::isfinite(tin_deg) ||
        tin_deg < 0.0 || tout_deg > 90.0 || tout_deg <= tin_deg ||
        (region != 0.0 && region != 1.0)) {
        throw FunctionUtility::FunctionException(
            "windlineoff: require finite uview, "
            "0 <= tin < tout <= 90 deg, and region = 0 or 1.");
    }

    RealArray physicalParams(params);
    physicalParams[1] = region == 0.0
        ? (1.0 - uview) * tin_deg
        : tout_deg + uview * (90.0 - tout_deg);
    physicalParams[4] = params[5];
    windline(energy, physicalParams, spectrum, flux, fluxErr, init);
}
