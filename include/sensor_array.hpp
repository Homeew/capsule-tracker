#pragma once
#include "dipole_model.hpp"
#include <vector>

namespace mag {

constexpr double kPi = 3.14159265358979323846;

/// Sensor positions [m]; each sensor is a 3-axis magnetometer.
using SensorArray = std::vector<Vec3>;

/// Two rings of sensors around the z axis (a wearable belt layout).
/// The upper ring is rotated by half a step so the rings interleave.
inline SensorArray twoLayerRing(int perRing = 8, double radius = 0.15,
                                double zLower = -0.05, double zUpper = 0.05) {
    SensorArray sensors;
    sensors.reserve(2 * perRing);
    const double step = 2.0 * kPi / perRing;
    for (int layer = 0; layer < 2; ++layer) {
        const double z      = (layer == 0) ? zLower : zUpper;
        const double offset = 0.5 * step * layer;
        for (int i = 0; i < perRing; ++i) {
            const double angle = i * step + offset;
            sensors.emplace_back(radius * std::cos(angle),
                                 radius * std::sin(angle), z);
        }
    }
    return sensors;
}

/// Noise-free stacked measurement vector [µT], three entries per sensor.
inline Eigen::VectorXd measure(const SensorArray& sensors,
                               const Vec3& p, const Vec3& m) {
    Eigen::VectorXd z(3 * sensors.size());
    for (std::size_t i = 0; i < sensors.size(); ++i)
        z.segment<3>(3 * i) = dipoleField(sensors[i], p, m);
    return z;
}

}  // namespace mag
