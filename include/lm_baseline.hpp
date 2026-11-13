#pragma once
#include "sensor_array.hpp"

namespace mag {

/// Output of one Levenberg-Marquardt solve.
struct LmResult {
    Vec3   position;     ///< estimated position [m]
    Vec3   moment;       ///< estimated moment [A·m²]
    int    iterations;   ///< iterations used
    double cost;         ///< final 0.5 * ||residual||² [µT²]
    bool   converged;
};

/// Frame-by-frame estimator: fits position and moment to one measurement
/// vector using Levenberg-Marquardt with the analytical Jacobian.
class LmBaseline {
public:
    explicit LmBaseline(SensorArray sensors, int maxIterations = 100,
                        double stepTolerance = 1e-10);

    /// Solve one frame, starting from (p0, m0).
    LmResult solve(const Eigen::VectorXd& z,
                   const Vec3& p0, const Vec3& m0) const;

private:
    SensorArray sensors_;
    int         maxIterations_;
    double      stepTolerance_;
};

}  // namespace mag
