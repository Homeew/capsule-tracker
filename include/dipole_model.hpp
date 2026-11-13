#pragma once
#include <Eigen/Dense>
#include <cmath>

/// Magnetic dipole model. Lengths in metres, moments in A·m², fields in µT.
namespace mag {

using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;

/// mu0 / (4*pi) expressed in µT·m/A (1e-7 T·m/A).
constexpr double kMu0Over4Pi = 0.1;

/// Field [µT] at `sensor` from a dipole at `p` with moment `m`.
inline Vec3 dipoleField(const Vec3& sensor, const Vec3& p, const Vec3& m) {
    const Vec3   r  = sensor - p;
    const double R2 = r.squaredNorm();
    const double R  = std::sqrt(R2);
    const double R3 = R2 * R;
    const double R5 = R3 * R2;
    return kMu0Over4Pi * (3.0 * m.dot(r) * r / R5 - m / R3);
}

/// Analytical derivatives of the field with respect to dipole position
/// (dBdp, µT/m) and moment (dBdm, µT per A·m²).
inline void dipoleJacobian(const Vec3& sensor, const Vec3& p, const Vec3& m,
                           Mat3& dBdp, Mat3& dBdm) {
    const Vec3   r  = sensor - p;
    const double R2 = r.squaredNorm();
    const double R  = std::sqrt(R2);
    const double R3 = R2 * R;
    const double R5 = R3 * R2;
    const double R7 = R5 * R2;
    const double mr = m.dot(r);
    const Mat3   I  = Mat3::Identity();

    // The field is linear in m.
    dBdm = kMu0Over4Pi * (3.0 * r * r.transpose() / R5 - I / R3);

    // dB/dr, then dB/dp = -dB/dr because r = sensor - p.
    const Mat3 dBdr = kMu0Over4Pi *
        (3.0 * (r * m.transpose() + m * r.transpose() + mr * I) / R5
         - 15.0 * mr * r * r.transpose() / R7);
    dBdp = -dBdr;
}

}  // namespace mag
