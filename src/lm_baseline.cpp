#include "lm_baseline.hpp"
#include <algorithm>

namespace mag {

using Vec6 = Eigen::Matrix<double, 6, 1>;
using Mat6 = Eigen::Matrix<double, 6, 6>;

LmBaseline::LmBaseline(SensorArray sensors, int maxIterations,
                       double stepTolerance)
    : sensors_(std::move(sensors)),
      maxIterations_(maxIterations),
      stepTolerance_(stepTolerance) {}

LmResult LmBaseline::solve(const Eigen::VectorXd& z,
                           const Vec3& p0, const Vec3& m0) const {
    const int n = 3 * static_cast<int>(sensors_.size());

    auto residual = [&](const Vec6& x) -> Eigen::VectorXd {
        return measure(sensors_, x.head<3>(), x.tail<3>()) - z;
    };

    auto jacobian = [&](const Vec6& x) -> Eigen::MatrixXd {
        Eigen::MatrixXd J(n, 6);
        Mat3 dBdp, dBdm;
        for (int i = 0; i < static_cast<int>(sensors_.size()); ++i) {
            dipoleJacobian(sensors_[i], x.head<3>(), x.tail<3>(), dBdp, dBdm);
            J.block<3, 3>(3 * i, 0) = dBdp;
            J.block<3, 3>(3 * i, 3) = dBdm;
        }
        return J;
    };

    Vec6 x;
    x << p0, m0;
    Eigen::VectorXd r = residual(x);
    double cost   = 0.5 * r.squaredNorm();
    double lambda = 1e-3;

    LmResult result{};
    int it = 0;
    for (; it < maxIterations_; ++it) {
        const Eigen::MatrixXd J = jacobian(x);
        const Mat6 A = J.transpose() * J;
        const Vec6 g = J.transpose() * r;

        // Marquardt damping: scale by the diagonal so position and moment
        // columns (different units) are damped fairly.
        Mat6 D = Mat6::Zero();
        D.diagonal() = A.diagonal().cwiseMax(1e-12);
        const Vec6 dx = (A + lambda * D).ldlt().solve(-g);

        if (dx.norm() < stepTolerance_ * (1.0 + x.norm())) {
            result.converged = true;
            break;
        }

        const Vec6 xNew = x + dx;
        const Eigen::VectorXd rNew = residual(xNew);
        const double costNew = 0.5 * rNew.squaredNorm();

        if (costNew < cost) {          // accept the step
            x = xNew;
            r = rNew;
            cost = costNew;
            lambda = std::max(lambda / 10.0, 1e-12);
        } else {                       // reject and damp harder
            lambda *= 10.0;
            if (lambda > 1e12) break;
        }
    }

    result.position   = x.head<3>();
    result.moment     = x.tail<3>();
    result.iterations = it;
    result.cost       = cost;
    return result;
}

}  // namespace mag
