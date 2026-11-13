#include "trajectory.hpp"
#include <cmath>

namespace mag {

namespace {
constexpr double kRate = 0.15;   // mean path-parameter rate [rad/s]
constexpr double kMod  = 0.5;    // speed modulation depth (0..1)
constexpr double kNu   = 0.4;    // speed modulation frequency [rad/s]
}  // namespace

GiTrajectory::GiTrajectory(double momentMagnitude)
    : momentMagnitude_(momentMagnitude) {}

CapsuleState GiTrajectory::at(double t) const {
    // Path parameter s(t) and its rate; ds/dt > 0 so the capsule never reverses.
    const double s  = kRate * (t + kMod * std::sin(kNu * t) / kNu);
    const double ds = kRate * (1.0 + kMod * std::cos(kNu * t));

    // Curve c(s) and its tangent dc/ds. The tangent is never zero.
    const Vec3 c(0.05 * std::sin(s),
                 0.04 * std::sin(2.0 * s + 0.5),
                 0.03 * std::cos(0.5 * s));
    const Vec3 dc(0.05 * std::cos(s),
                  0.08 * std::cos(2.0 * s + 0.5),
                  -0.015 * std::sin(0.5 * s));

    CapsuleState state;
    state.position = c;
    state.velocity = dc * ds;
    state.moment   = momentMagnitude_ * dc.normalized();
    return state;
}

}  // namespace mag
