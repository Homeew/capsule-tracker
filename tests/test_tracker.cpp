#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "dipole_model.hpp"

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
