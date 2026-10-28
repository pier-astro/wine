#pragma once

#include <XSFunctions/Utilities/FunctionUtility.h>
#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>
#include <string>

namespace wine {

// XSPEC bounds constrain accepted independent parameters, but numerical
// derivatives may probe outside them. Continue with the nearest physical
// boundary (a one-sided extension), rather than an inverted shell or NaN.
// Python constructors validate user inputs strictly instead of extending them.
inline double boundedParameter(double value, double lower, double upper,
                               const char* name, int spectrum) {
    if (!std::isfinite(value)) {
        throw FunctionUtility::FunctionException(
            std::string("WINE: non-finite parameter ") + name + ".");
    }
    if (value < lower || value > upper) {
        static std::set<std::pair<std::string, int>> warned;
        if (warned.insert({name, spectrum}).second) {
            std::ostringstream message;
            message << "\n*** WINE: " << name << "=" << value
                    << " is outside [" << lower << ", " << upper << "]. "
                    << "Evaluating its physical boundary for this numerical probe. "
                    << "Accepted fits and parameter links must respect the bounds.\n";
            FunctionUtility::xsWrite(message.str(), 10);
        }
    }
    return std::max(lower, std::min(upper, value));
}

} // namespace wine
