#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "dipole_model.hpp"
#include "sensor_array.hpp"
#include "trajectory.hpp"
#include "lm_baseline.hpp"

using mag::Vec3;
using mag::Mat3;

TEST_CASE("On-axis field matches the textbook value", "[dipole]") {
    // On the dipole axis: B = (mu0/4pi) * 2m / R^3
    const Vec3 B = mag::dipoleField({0, 0, 0.1}, {0, 0, 0}, {0, 0, 0.2});
    REQUIRE(B.x() == Approx(0.0).margin(1e-12));
    REQUIRE(B.y() == Approx(0.0).margin(1e-12));
    REQUIRE(B.z() == Approx(0.1 * 2.0 * 0.2 / 1e-3));   // 40 µT
}

TEST_CASE("Analytical Jacobian matches finite differences", "[dipole]") {
    const Vec3 sensor(0.15, 0.02, 0.05);
    const Vec3 p(0.01, -0.03, 0.02);
    const Vec3 m(0.05, -0.10, 0.16);

    Mat3 dBdp, dBdm;
    mag::dipoleJacobian(sensor, p, m, dBdp, dBdm);

    const double hp = 1e-6;   // position step [m]
    const double hm = 1e-4;   // moment step [A·m²]
    for (int j = 0; j < 3; ++j) {
        Vec3 dp = Vec3::Zero();  dp[j] = hp;
        Vec3 dm = Vec3::Zero();  dm[j] = hm;
        const Vec3 colP = (mag::dipoleField(sensor, p + dp, m) -
                           mag::dipoleField(sensor, p - dp, m)) / (2 * hp);
        const Vec3 colM = (mag::dipoleField(sensor, p, m + dm) -
                           mag::dipoleField(sensor, p, m - dm)) / (2 * hm);
        for (int i = 0; i < 3; ++i) {
            REQUIRE(dBdp(i, j) == Approx(colP[i]).epsilon(1e-6).margin(1e-6));
            REQUIRE(dBdm(i, j) == Approx(colM[i]).epsilon(1e-6).margin(1e-6));
        }
    }
}

TEST_CASE("Two-layer ring has the expected geometry", "[sensors]") {
    const mag::SensorArray sensors = mag::twoLayerRing(4, 0.15, -0.15, 0.15);
    REQUIRE(sensors.size() == 8);
    for (const Vec3& s : sensors) {
        REQUIRE(std::hypot(s.x(), s.y()) == Approx(0.15));
        REQUIRE(std::abs(s.z()) == Approx(0.15));
    }
}

TEST_CASE("Stacked measurement matches the per-sensor field", "[sensors]") {
    const mag::SensorArray sensors = mag::twoLayerRing();
    const Vec3 p(0.01, -0.02, 0.0);
    const Vec3 m(0.0, 0.1, 0.17);
    const Eigen::VectorXd z = mag::measure(sensors, p, m);
    REQUIRE(z.size() == 24);
    const Vec3 B5 = mag::dipoleField(sensors[5], p, m);
    for (int i = 0; i < 3; ++i)
        REQUIRE(z[15 + i] == Approx(B5[i]));
}

TEST_CASE("Trajectory is smooth and physically consistent", "[trajectory]") {
    const mag::GiTrajectory traj(0.2);
    const double h = 1e-5;
    for (double t = 0.0; t <= 60.0; t += 0.5) {
        const mag::CapsuleState s = traj.at(t);

        // Stays inside the sensor belt.
        REQUIRE(s.position.norm() < 0.10);

        // Velocity equals the time derivative of position.
        const Vec3 fd = (traj.at(t + h).position - traj.at(t - h).position) / (2 * h);
        REQUIRE((s.velocity - fd).norm() < 1e-7);

        // Realistic capsule speed: between 1 and 30 mm/s.
        REQUIRE(s.velocity.norm() > 0.001);
        REQUIRE(s.velocity.norm() < 0.030);

        // Moment has fixed magnitude and points along the direction of travel.
        REQUIRE(s.moment.norm() == Approx(0.2));
        REQUIRE(s.moment.cross(s.velocity).norm() < 1e-12);
        REQUIRE(s.moment.dot(s.velocity) > 0.0);
    }
}

TEST_CASE("LM recovers a static dipole from noiseless data", "[lm]") {
    const mag::SensorArray sensors = mag::twoLayerRing();
    const Vec3 pTrue(0.02, -0.01, 0.03);
    const Vec3 mTrue = 0.2 * Vec3(0.3, 0.4, 0.866).normalized();
    const Eigen::VectorXd z = mag::measure(sensors, pTrue, mTrue);

    const mag::LmBaseline lm(sensors);
    const mag::LmResult res =
        lm.solve(z, pTrue + Vec3(0.01, -0.01, 0.01), Vec3(0.0, 0.0, 0.2));

    REQUIRE(res.converged);
    REQUIRE((res.position - pTrue).norm() < 1e-6);   // under 1 µm
    REQUIRE((res.moment - mTrue).norm() < 1e-6);
}

TEST_CASE("Warm-started LM tracks the trajectory without noise", "[lm]") {
    const mag::SensorArray sensors = mag::twoLayerRing();
    const mag::GiTrajectory traj(0.2);
    const mag::LmBaseline lm(sensors);

    mag::CapsuleState s0 = traj.at(0.0);
    Vec3 p = s0.position + Vec3(0.005, 0.005, -0.005);
    Vec3 m = s0.moment;
    for (double t = 0.0; t <= 30.0; t += 0.5) {
        const mag::CapsuleState s = traj.at(t);
        const mag::LmResult res =
            lm.solve(mag::measure(sensors, s.position, s.moment), p, m);
        REQUIRE(res.converged);
        REQUIRE((res.position - s.position).norm() < 1e-6);
        p = res.position;    // warm start the next frame
        m = res.moment;
    }
}
