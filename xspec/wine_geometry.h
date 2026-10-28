#pragma once

#include <algorithm>
#include <cmath>

namespace wine {

// Physical angles are radians. Cf = Omega/(2*pi) for one hemisphere.
struct Geometry {
    double incl;
    double tin;
    double tout;
    double Cf;
};

inline Geometry freeGeometry(double Cf, double hole, double incl) {
    const double cos_in = 1.0 - hole * (1.0 - Cf);
    const double cos_out = (1.0 - hole) * (1.0 - Cf);
    return {incl, std::acos(cos_in), std::acos(cos_out), Cf};
}

inline Geometry losGeometry(double Cf, double skew, double incl) {
    const double half_pi = std::acos(0.0);
    const double c = incl == half_pi ? 0.0 : std::cos(incl);
    const double a = 1.0 - 2.0 * skew;
    const double mix = (1.0 - skew) * (1.0 - Cf) + skew * Cf;
    const double b = -(mix + a * c);
    const double constant = (1.0 - skew) * c * (1.0 - Cf);

    // For x=cos(tout), solve
    // (1-skew)*(1-x-Cf)*(c-x) = skew*x*(x+Cf-c).
    // This positive-term discriminant stays accurate at the double roots.
    const double discriminant = (mix - c) * (mix - c)
        + 4.0 * skew * (1.0 - skew) * c * (1.0 - c);
    const double root = std::sqrt(discriminant);
    double x = 0.0;
    if (b > 0.0) {
        x = (-b - root) / (2.0 * a);
    } else if (-b + root > 0.0) {
        // Rationalized physical root; also valid for skew=0.5 (a=0).
        x = 2.0 * constant / (-b + root);
    }
    const double lower = std::max(0.0, c - Cf);
    const double upper = std::min(c, 1.0 - Cf);
    x = std::max(lower, std::min(upper, x)); // Round-off at the physical bounds.
    const double cos_in = std::max(0.0, std::min(1.0, x + Cf));
    return {incl, std::acos(cos_in), std::acos(x), Cf};
}

inline Geometry offGeometry(double Cf, double hole, double view, int region) {
    Geometry geometry = freeGeometry(Cf, hole, 0.0);
    geometry.incl = region == 0
        ? view * geometry.tin
        : geometry.tout + view * (std::acos(0.0) - geometry.tout);
    return geometry;
}

} // namespace wine
