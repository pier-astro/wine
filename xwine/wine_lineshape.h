#pragma once

#include <cmath>
#include <algorithm>

namespace wine {

constexpr double kPi = 3.141592653589793;
constexpr double kDeg2Rad = kPi / 180.0;

inline double w_phi_hollowcap(double theta, double incl, double theta_in, double theta_out) {
    double costheta = std::cos(theta);
    double sintheta = std::sin(theta);
    double cosinc = std::cos(incl);
    double sininc = std::sin(incl);

    if (std::abs(sintheta) < 1e-10 || std::abs(sininc) < 1e-10) {
        double cosang = costheta * cosinc;
        if (cosang >= std::cos(theta_out) && cosang <= std::cos(theta_in)) {
            return 2.0 * kPi;
        }
        return 0.0;
    }

    double C_in = (cosinc * costheta - std::cos(theta_in)) / (sininc * sintheta);
    double C_out = (cosinc * costheta - std::cos(theta_out)) / (sininc * sintheta);

    double w_in, w_out;

    if (C_in >= 1.0) {
        w_in = 2.0 * kPi;
    } else if (C_in <= -1.0) {
        w_in = 0.0;
    } else {
        w_in = 2.0 * (kPi - std::acos(C_in));
    }

    if (C_out >= 1.0) {
        w_out = 2.0 * kPi;
    } else if (C_out <= -1.0) {
        w_out = 0.0;
    } else {
        w_out = 2.0 * (kPi - std::acos(C_out));
    }

    return std::max(0.0, w_out - w_in);
}

inline void lineshape(const double* E, int n, double E0, double beta, double incl,
                      double theta_in, double theta_out, double* result) {
    double sqrt_1_beta2 = std::sqrt(1.0 - beta * beta);
    const double fourpi = 4.0 * kPi;

    for (int i = 0; i < n; ++i) {
        double E_E0 = (E[i] > 0.0) ? (E[i] / E0) : 0.0;
        double costheta = (1.0 - sqrt_1_beta2 / E_E0) / beta;

        if (costheta < -1.0 || costheta > 1.0) {
            result[i] = 0.0;
            continue;
        }

        double theta = std::acos(costheta);
        double w = w_phi_hollowcap(theta, incl, theta_in, theta_out);

        if (w <= 0.0) {
            result[i] = 0.0;
            continue;
        }

        result[i] = sqrt_1_beta2 / (beta * E0) * w / fourpi;
    }
}

}  // namespace wine
