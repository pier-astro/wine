#include <xsTypes.h>
#include <XSUtil/Numerics/LinearInterp.h>
#include <XSUtil/Numerics/Numerics.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <cmath>
#include <set>
#include <string>

// Convolution models for velocity-shifted absorption
// Usage: wind*mtable{table.fits}*powerlaw  (velocity only)
//        zwind*mtable{table.fits}*powerlaw (velocity + redshift)

static void warnFiniteGrid(const char* modelName) {
    static std::set<std::string> warnedModels;
    if (warnedModels.insert(modelName).second) {
        FunctionUtility::xsWrite("\n\t*** WINE " + std::string(modelName) +
            " ENERGY-GRID NOTICE ***\n\tFinite grids require endpoint "
            "extrapolation at their shifted edge. Use `energies extend low ...` "
            "and/or `energies extend high ...` so that edge lies outside the detector response.\n",
            10);
    }
}

static void applyVelocityShift(const RealArray& energy, RealArray& flux, Real doppler,
                               const char* modelName, RealArray& fluxErr,
                               bool additive = false, Real fluxScale = 1.0) {
    using namespace Numerics;
    using namespace Rebin;

    if (flux.size() == 0) {
        return;
    }
    if (doppler == 1.0) {
        flux *= fluxScale;
        if (fluxErr.size() > 0) fluxErr *= fluxScale;
        return;
    }
    if (energy.size() != flux.size() + 1 || doppler <= 0.0 || !std::isfinite(doppler)) {
        throw FunctionUtility::FunctionException(std::string(modelName) +
            ": received an invalid energy grid or Doppler factor.");
    }

    warnFiniteGrid(modelName);

    const RealArray shiftedEnergy(energy * doppler);
    size_t inputBin;
    size_t outputBin;
    IntegerVector startBin(flux.size());
    IntegerVector endBin(flux.size());
    RealArray startWeight(flux.size());
    RealArray endWeight(flux.size());
    const Real fuzzy = 1.0e-6;

    findFirstBins(shiftedEnergy, energy, fuzzy, inputBin, outputBin);
    initializeBins(shiftedEnergy, energy, fuzzy, inputBin, outputBin,
                   startBin, endBin, startWeight, endWeight);

    RealArray shiftedFlux(flux);
    if (additive) {
        rebin(flux, startBin, endBin, startWeight, endWeight, shiftedFlux);
    } else {
        interpolate(flux, startBin, endBin, startWeight, endWeight,
                    shiftedFlux, false);
    }
    flux = shiftedFlux * fluxScale;

    if (fluxErr.size() > 0) {
        RealArray shiftedFluxErr(fluxErr);
        if (additive) {
            rebin(fluxErr, startBin, endBin, startWeight, endWeight, shiftedFluxErr);
        } else {
            interpolate(fluxErr, startBin, endBin, startWeight, endWeight,
                        shiftedFluxErr, false);
        }
        fluxErr = shiftedFluxErr * fluxScale;
    }
}

extern "C" void wind(const RealArray& energy, const RealArray& params,
                     int spectrum, RealArray& flux, RealArray& fluxErr,
                     const std::string& init) {
    const Real vout = params[0];
    const Real beta = vout;
    // Positive outflow velocity moves rest-frame features to higher energy.
    const Real doppler = std::sqrt((1.0 + beta) / (1.0 - beta));
    applyVelocityShift(energy, flux, doppler, "wind", fluxErr);
}

extern "C" void zwind(const RealArray& energy, const RealArray& params,
                      int spectrum, RealArray& flux, RealArray& fluxErr,
                      const std::string& init) {
    const Real vout = params[0];
    const Real z = params[1];
    const Real beta = vout;
    // Cosmological redshift opposes the positive-velocity blueshift.
    const Real doppler = std::sqrt((1.0 + beta) / (1.0 - beta)) / (1.0 + z);
    applyVelocityShift(energy, flux, doppler, "zwind", fluxErr);
}

extern "C" void awind(const RealArray& energy, const RealArray& params,
                      int spectrum, RealArray& flux, RealArray& fluxErr,
                      const std::string& init) {
    const Real beta = params[0];
    const Real doppler = std::sqrt((1.0 + beta) / (1.0 - beta));
    applyVelocityShift(energy, flux, doppler, "awind", fluxErr, true);
}

extern "C" void azwind(const RealArray& energy, const RealArray& params,
                       int spectrum, RealArray& flux, RealArray& fluxErr,
                       const std::string& init) {
    const Real beta = params[0];
    const Real zfactor = 1.0 / (1.0 + params[1]);
    const Real doppler = std::sqrt((1.0 + beta) / (1.0 - beta)) * zfactor;
    applyVelocityShift(energy, flux, doppler, "azwind", fluxErr, true, zfactor);
}
