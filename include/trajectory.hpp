#pragma once
#include "dipole_model.hpp"

namespace mag {

/// Ground-truth capsule state at one instant.
struct CapsuleState {
    Vec3 position;   ///< [m]
    Vec3 velocity;   ///< [m/s]
    Vec3 moment;     ///< [A·m²], aligned with the direction of travel
};

/// Smooth, looping 3D path that mimics a capsule moving through the
/// intestine with a pulsating (peristalsis-like) speed.
class GiTrajectory {
public:
    explicit GiTrajectory(double momentMagnitude = 0.2);

    /// Ground truth at time t [s].
    CapsuleState at(double t) const;

    double momentMagnitude() const { return momentMagnitude_; }

private:
    double momentMagnitude_;
};

}  // namespace mag
