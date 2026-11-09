#pragma once

#include <xsTypes.h>
#include "redwine_geometry.h"
#include "redwine_lineshape.h"
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>
#include <string>

namespace redwine {

// XSPEC bounds constrain accepted independent parameters, but numerical
// derivatives may probe outside them. Continue with the nearest physical
// boundary (a one-sided extension), rather than an inverted shell or NaN.
// Python constructors validate user inputs strictly instead of extending them.
inline double boundedParameter(double value, double lower, double upper,
                               const char* name, int spectrum) {
    if (!std::isfinite(value)) {
        throw FunctionUtility::FunctionException(
            std::string("RED-WINE: non-finite parameter ") + name + ".");
    }
    if (value < lower || value > upper) {
        static std::set<std::pair<std::string, int>> warned;
        if (warned.insert({name, spectrum}).second) {
            std::ostringstream message;
            message << "\n*** RED-WINE: " << name << "=" << value
                    << " is outside [" << lower << ", " << upper << "]. "
                    << "Evaluating its physical boundary for this numerical probe. "
                    << "Accepted fits and parameter links must respect the bounds.\n";
            FunctionUtility::xsWrite(message.str(), 10);
        }
    }
    return std::max(lower, std::min(upper, value));
}


inline Geometry angularGeometry(double incl, double tin, double tout) {
    if (tin >= tout) {
        static const bool warned = [] {
            FunctionUtility::xsWrite(
                "RED-WINE: tin >= tout. The component returns zero emission.\n", 10);
            return true;
        }();
        (void)warned;
        return {incl, tin, tout, 0.0};
    }
    return {incl, tin, tout, std::cos(tin) - std::cos(tout)};
}

inline Geometry freeCf(const RealArray& params, int spectrum) {
    const double Cf = boundedParameter(params[1], 0.0, 1.0, "Cf", spectrum);
    const double hole = boundedParameter(params[2], 0.0, 1.0, "hole", spectrum);
    const double incl = boundedParameter(params[3], 0.0, 90.0, "incl", spectrum) * kDeg2Rad;
    return freeGeometry(Cf, hole, incl);
}

inline Geometry losCf(const RealArray& params, int spectrum) {
    const double Cf = boundedParameter(params[1], 0.0, 1.0, "Cf", spectrum);
    const double skew = boundedParameter(params[2], 0.0, 1.0, "skew", spectrum);
    const double incl = boundedParameter(params[3], 0.0, 90.0, "incl", spectrum) * kDeg2Rad;
    return losGeometry(Cf, skew, incl);
}

inline Geometry offCf(const RealArray& params, int spectrum) {
    const double Cf = boundedParameter(params[1], 0.0, 1.0, "Cf", spectrum);
    const double hole = boundedParameter(params[2], 0.0, 1.0, "hole", spectrum);
    const double view = boundedParameter(params[3], 0.0, 1.0, "view", spectrum);
    if (params[4] != 0.0 && params[4] != 1.0) {
        throw FunctionUtility::FunctionException("RED-WINE: region must be 0 or 1.");
    }
    return offGeometry(Cf, hole, view, static_cast<int>(params[4]));
}

inline Geometry freeangles(const RealArray& params, int spectrum) {
    const double tin = boundedParameter(params[1], 0.0, 90.0, "tin", spectrum) * kDeg2Rad;
    const double tout = boundedParameter(params[2], 0.0, 90.0, "tout", spectrum) * kDeg2Rad;
    const double incl = boundedParameter(params[3], 0.0, 90.0, "incl", spectrum) * kDeg2Rad;
    return angularGeometry(incl, tin, tout);
}

inline Geometry losangles(const RealArray& params, int spectrum) {
    const double tin = boundedParameter(params[1], 0.0, 90.0, "tin", spectrum) * kDeg2Rad;
    const double tout = boundedParameter(params[2], 0.0, 90.0, "tout", spectrum) * kDeg2Rad;
    const double uincl = boundedParameter(params[3], 0.0, 1.0, "uincl", spectrum);
    const double incl = tin + uincl * (tout - tin);
    return angularGeometry(incl, tin, tout);
}

inline Geometry offangles(const RealArray& params, int spectrum) {
    const double tin = boundedParameter(params[1], 0.0, 90.0, "tin", spectrum) * kDeg2Rad;
    const double tout = boundedParameter(params[2], 0.0, 90.0, "tout", spectrum) * kDeg2Rad;
    const double uview = boundedParameter(params[3], 0.0, 1.0, "uview", spectrum);
    if (params[4] != 0.0 && params[4] != 1.0) {
        throw FunctionUtility::FunctionException("RED-WINE: region must be 0 or 1.");
    }
    const double incl = params[4] == 0.0 ? (1.0 - uview) * tin
        : tout + uview * (90.0 * kDeg2Rad - tout);
    return angularGeometry(incl, tin, tout);
}

} // namespace redwine
