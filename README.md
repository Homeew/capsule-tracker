# Capsule Tracker

![CI](https://github.com/Homeew/capsule-tracker/actions/workflows/ci.yml/badge.svg)

Real-time magnetic tracking of a moving capsule endoscope in C++17:
Extended and Unscented Kalman filters versus frame-by-frame
Levenberg-Marquardt on a simulated GI trajectory.

## Motivation

In wireless capsule endoscopy, a capsule carrying a small permanent magnet
travels through the GI tract while an external magnetometer array measures its
field. Most localization pipelines solve each frame independently with
nonlinear least squares. This project treats the capsule as a **moving
target**: a Kalman filter fuses the measurements over time with a motion
model, which reduces noise and also estimates velocity.

It is the dynamic follow-up to
[magnetic-dipole-simulator](https://github.com/Homeew/magnetic-dipole-simulator),
which covers the static, single-frame problem.

## Status

- [x] Dipole field model with analytical Jacobian (position and moment)
- [x] Two-layer ring magnetometer array (16 three-axis sensors)
- [x] Simulated GI trajectory with pulsating, peristalsis-like speed
- [ ] Frame-by-frame Levenberg-Marquardt baseline
- [ ] Extended Kalman filter
- [ ] Unscented Kalman filter
- [ ] Noise-sweep benchmark and plots

## Model

**State (9):** position `p` [m], velocity `v` [m/s], magnetic moment `m` [A·m²].

**Motion:** constant velocity with white-acceleration process noise; the
moment follows a random walk.

**Measurement:** each three-axis sensor at `s` observes the dipole field

    B(s) = (mu0 / 4 pi) * [ 3 (m . r) r / |r|^5 - m / |r|^3 ],   r = s - p

Fields are expressed in microtesla throughout, which keeps the filter
covariances well conditioned.

## Simulation setup

| Parameter | Value |
|-----------|-------|
| Sensors | 16, two interleaved rings, radius 0.15 m, z = ±0.05 m |
| Magnet moment | 0.2 A·m², aligned with the direction of travel |
| Capsule speed | about 2 to 21 mm/s |
| Workspace | within 0.10 m of the array centre |

## Build

```bash
cmake -S . -B build
cmake --build build --parallel 4
```

Requires g++ with C++17 and CMake 3.18 or newer. Eigen 3.4 and Catch2 are
downloaded automatically.

## Test

```bash
./build/tests
```

The tests check the analytical Jacobian against finite differences, the
sensor geometry, and the physical consistency of the trajectory.

## Project structure

    include/
      dipole_model.hpp   - dipole field and analytical Jacobian
      sensor_array.hpp   - two-layer ring array and stacked measurements
      trajectory.hpp     - simulated GI trajectory
    src/
      trajectory.cpp
    tests/
      test_tracker.cpp   - Catch2 unit tests

## Research context

This project builds on the author's published work on magnetic localization
for wireless capsule endoscopy:

- O. Yaghoobian, M. Einlou, K.A. Wahid. *Two-stage cylindrical magnetic
  localization with neural network calibration for wearable capsule
  endoscopy.* Scientific Reports, 2026.
- O. Yaghoobian, H. Keshmiri, M. Einlou, K.A. Wahid. *Layer-Conditioned
  Magnetic Localization with Residual-Weighted Fusion for Wireless Capsule
  Endoscopy Under Asymmetric Sensor Visibility and GI Motion.* IEEE Access, 2026.
- O. Yaghoobian, K.A. Wahid. *SNR-Weighted Layer-Conditioned Magnetometer
  Array with Ellipsoid Calibration and Neural Network Initialization for
  Capsule Endoscopy Localization Under Asymmetric Sensor Visibility.*
  Sensors, 2026.
- O. Yaghoobian, K.A. Wahid. *A New Localization Algorithm for a Two-Layer
  Magnetic Sensor Array in Capsule Endoscopy.* ICCIT 2025.
