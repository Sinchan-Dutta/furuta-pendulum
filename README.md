# Furuta Pendulum: From Swing-Up + LQR to Constrained Nonlinear MPC with EKF State Estimation

This repository documents the development of a control stack for a Furuta Pendulum (rotary inverted pendulum). It contains two tracks:

- **MATLAB**: a classical switched controller, with an energy-based swing-up followed by an LQR stabilizer, on a plant that includes the full DC motor electrical dynamics.
- **Python**: a unified Nonlinear Model Predictive Controller (NMPC) built with CasADi and IPOPT, extended to a closed loop with noisy sensors in which the NMPC acts on the state estimate of an Extended Kalman Filter (EKF).

The project was developed as part of the Controls and Systems course (some developments were made later), with an emphasis on physically realistic modeling (a DC motor actuator model, voltage limits, explicit state constraints, sensor noise and actuator latency) rather than a purely textbook simulation.

**Note on the two implementations**: the MATLAB and Python controllers were developed somewhat independently rather than as a single linear progression. The Python line (`python/`) is the primary, original line of development, including the parameter sourcing, the constrained optimal control formulation, and the noise and estimation extensions. The MATLAB controller (`matlab/`) is a separate reference exploration of a more classical architecture, with the actuator modeled in more electrical detail (armature inductance and motor current as an explicit fifth state, rather than the quasi-static voltage-to-torque relationship used in Python). It is included as a point of comparison, not as an earlier version of the Python implementation. The two use different physical and actuator parameters (see "Differences between the MATLAB and Python tracks").

## Project structure

```
furuta-pendulum/
├── docs/
│   └── pendulum_dynamics.pdf             # Handwritten derivation of the dynamics, linearization and motor model
├── matlab/
│   ├── src/
│   │   └── swingup_lqr_controller.m
│   └── results/
│       ├── furuta_motor.mp4
│       └── state_variation_with_time.png
├── python/
│   ├── src/
│   │   ├── dynamics.py                   # CasADi symbolic nonlinear EOM + RK4 integrator
│   │   ├── motor_model.py                # DC motor voltage-to-torque model
│   │   ├── noise_models.py               # Gaussian sensor noise model
│   │   ├── ekf.py                        # Extended Kalman Filter
│   │   ├── unified_nmpc_without_noise.py # NMPC with perfect state feedback
│   │   └── unified_nmpc_with_noise.py    # NMPC + noisy sensors + EKF
│   └── results/
│       ├── nmpc_closed_loop.mp4
│       ├── nmpc_results_without_noise.png
│       ├── nmpc_closed_loop_ekf.mp4
│       ├── nmpc_results_ekf.png
│       └── nmpc_results_with_noise.png
├── requirements.txt
└── README.md
```

## Background

A Furuta Pendulum consists of a horizontal arm, driven by a single motor at its base, with an unactuated pendulum rod attached to the far end of the arm. The control objective is to swing the pendulum up from its hanging-down rest position and stabilize it in the inverted (upright) position, using only the single arm-side actuator. Since the pendulum has no actuator of its own, it can only be driven indirectly, through the dynamic coupling between the arm's rotation and the pendulum's swing.

This system is a standard benchmark in nonlinear and underactuated control because it captures, in a small and analyzable package, many challenges found in more complex underactuated systems: a large-angle, genuinely nonlinear swing-up maneuver, followed by precision stabilization at an unstable equilibrium, all while respecting real actuator limits.

## Physical model

The arm angle is denoted θ and the pendulum angle φ, with φ = 0 the hanging-down position and φ = π the upright (inverted) equilibrium. `docs/pendulum_dynamics.pdf` contains the derivation: the Euler–Lagrange method applied to the kinetic and potential energy of the centre rod, arm and pendulum, the resulting coupled nonlinear ODEs, their linearization about the equilibria (hanging and inverted), and the corresponding transfer functions. It also derives the actuator model in both forms used in this repository: the full electrical model with armature inductance (used in MATLAB), and the quasi-static simplification obtained by assuming the electrical time constant is much faster than the control sample time, i.e. `Lm/Rm << sample time` (used in Python). The equations of motion are:

```
Den      = α·β + β²·sin²(φ) − (γ·cos(φ))²
U        = −β·sin(2φ)·θ̇·φ̇ + γ·sin(φ)·φ̇² − c·θ̇ + τ
V        = 0.5·β·sin(2φ)·θ̇² − δ·sin(φ) − b·φ̇
θ̈        = (β·U − γ·cos(φ)·V) / Den
φ̈        = (−γ·cos(φ)·U + (α + β·sin²(φ))·V) / Den
```

where c and b are the arm and pendulum viscous damping coefficients, τ is the torque applied at the arm shaft, and

```
α = Jr + mp·la²      β = Jp      γ = ½·mp·la·lp      δ = ½·mp·g·lp
```

with `Jr` the arm (base) inertia and `Jp` the pendulum inertia about its pivot. The two tracks use the same equations with different numerical parameters:

- **Python** takes the arm and pendulum dimensions, masses, inertias and damping from published Quanser QUBE-Servo 2 specifications (`la = 0.085 m`, `lp = 0.129 m`, `mp = 0.024 kg`, measured `Jr = 2.306e-4`, `Jp = 1.3313e-4 kg·m²`).
- **MATLAB** uses illustrative values that are not tied to a specific platform (`la = 0.25 m`, `ma = 0.3 kg`, `lp = 0.15 m`, `mp = 0.05 kg`, base inertia `J = 1e-3 kg·m²`), with `α = J + ma·la²/3 + mp·la²` and `β = mp·lp²/3`, i.e. the arm and pendulum modeled as uniform thin rods.

## Actuator model

Torque is not a free input. The arm is driven by a DC gearmotor with full electrical dynamics:

```
Lm·(di/dt) = Vin − Rm·i − Kb·N·θ̇
τ          = N·Kt·i
```

where `Kt` and `Kb` are the torque and back-EMF constants, `Rm` and `Lm` the armature resistance and inductance, and `N` the gear ratio. Assuming `Lm·di/dt ≈ 0`, this reduces to the quasi-static relationship

```
τ = N·Kt·(Vin − Kb·N·θ̇) / Rm
```

Either way, the achievable torque depends on the arm's speed (back-EMF reduces available torque as speed increases), and the quantity that is actually bounded is voltage, not torque. `python/src/motor_model.py` implements the quasi-static form; the MATLAB model carries the current as a fifth state.

## MATLAB: swing-up + LQR

`matlab/src/swingup_lqr_controller.m` simulates the full nonlinear plant with the electrical dynamics, with the state `[θ, θ̇, φ, φ̇, i]`, using `ode45` (`RelTol = 1e-6`, `AbsTol = 1e-7`) over 6 s from a pendulum slightly displaced from hanging (`φ = 0.1 rad`, zero current). The controller switches between two modes depending on the pendulum angle error `φ_err = wrap(φ − π)`:

- **Energy-based swing-up** (`|φ_err| ≥ 40°`): the pendulum energy `E = ½·β·φ̇² + δ·(1 − cos φ)` is driven towards the upright energy `E_d = 2δ` with `τ_swing = kE·(E_d − E)·tanh(φ̇·cos φ / ε)` (`kE = 30`, `ε = 0.1`), plus a centering term `−1.5·θ − 0.05·θ̇` that keeps the arm near the middle of its range. The desired torque is converted to a feed-forward voltage through the quasi-static motor equations, `V = Rm·τ/(N·Kt) + Ke·N·θ̇`.
- **LQR stabilization** (`|φ_err| < 40°`): a state-feedback law `V = −K·[θ, θ̇, φ_err, φ̇, i]ᵀ` designed on the linearization about the upright equilibrium, including the motor current as a fifth state. The Riccati solution is computed with the Hamiltonian eigenvector method rather than a black-box call, to keep it explicit. The weights are `Q = diag(1, 1, 10, 1, 0.01)` and `R = 0.5`, and a controllability check (rank of the controllability matrix) is run before the design.

Voltage is saturated at ±24 V in both modes.

## Python: unified Nonlinear MPC

`python/src/` contains a from-scratch CasADi/IPOPT implementation that replaces the switched architecture with a single nonlinear optimal control problem, solved at every control step over a finite receding horizon.

- **`dynamics.py`** builds the equations of motion as a CasADi symbolic function, and a fourth-order Runge–Kutta integrator on top of it. Working symbolically allows the same model to be reused, unmodified, inside the NMPC's prediction, as the simulated plant, and (through automatic differentiation) as the source of the EKF's Jacobian.
- **`motor_model.py`** converts the NMPC's voltage decision variable into the torque the dynamics model expects, so the optimizer solves for a realistic, voltage-bounded input.
- **`unified_nmpc_without_noise.py`** is the NMPC with perfect state feedback.
- **`unified_nmpc_with_noise.py`** is the same controller in a closed loop with noisy sensors and an EKF.

### NMPC formulation

At every step the controller solves, over a horizon of N = 20 steps at dt = 0.025 s (0.5 s lookahead), for the state trajectory `X` and the voltage sequence `U`.

**Cost**, summed over the horizon:

```
Σ_k  Qθ·θ_k² + Qθ̇·θ̇_k² + Qφ·(1 + cos φ_k) + (k/N)²·Qφ̇·φ̇_k²   +   Σ_k R_u·u_k²   +   Q_T·(1 + cos φ_N)
```

- The upright cost `1 + cos φ` is zero at φ = π and handles angle wraparound correctly.
- The pendulum-velocity penalty is ramped by `(k/N)²`, so fast motion is cheap early in the horizon (needed to swing up) and expensive near its end.
- A heavy terminal cost `Q_T = 1e5` on `1 + cos φ_N` pushes the optimizer to finish the horizon upright.
- Weights: Qθ = 40, Qθ̇ = 20, Qφ = 1, Qφ̇ = 5, R_u = 0.1 (without noise); Qθ = 50, Qθ̇ = 15, Qφ = 1, Qφ̇ = 5, R_u = 0.5 (with noise).

**Constraints:**

- The true nonlinear dynamics (RK4, with the motor model in the loop) as an equality constraint at every step, and the initial state equal to the current state (or estimate).
- Arm angle |θ| ≤ π, pendulum angle |φ| ≤ 2π, and actuator voltage |V| ≤ 4 V.

**Solver:** IPOPT through CasADi, built once and reused, with the initial state as a parameter. Warm starting is enabled, with at most 50 iterations and a tolerance of 1e-4 for speed. Only the first voltage of each solution is applied (receding horizon). The solver is warm-started from the previous solution; the noiseless run uses a swing-up trajectory as the first guess, the noisy run starts from zeros.

Unlike the MATLAB controller, this one does not switch between separate swing-up and stabilization laws: a single optimization discovers both the energy-pumping swing-up and the final stabilizing approach as one continuous trajectory under the same constraints.

## Sensor noise model

`python/src/noise_models.py` provides `SensorModel`, which turns the true state `[θ, θ̇, φ, φ̇]` into a measurement by adding zero-mean Gaussian noise to all four states and then wrapping the two angles to [−π, π], as an encoder would report them.

| State | Noise standard deviation |
| ----- | ------------------------ |
| θ (arm angle)         | 0.002 rad   |
| φ (pendulum angle)    | 0.003 rad   |
| θ̇ (arm velocity)      | 0.05 rad/s  |
| φ̇ (pendulum velocity) | 0.08 rad/s  |

## Extended Kalman Filter

`python/src/ekf.py` implements an EKF over the four states `[θ, θ̇, φ, φ̇]`.

- **Predict:** the estimate is propagated with the same nonlinear RK4 model, using the previously applied torque. The covariance is propagated as `P ← F·P·Fᵀ + Q`, where `F ≈ I + A·dt` is a first-order discretization of the Jacobian `A = ∂f/∂x`, obtained symbolically with CasADi.
- **Update:** all four states are measured (`H = I`). The innovation is computed with the angular error wrapped to the shortest path in [−π, π], followed by the standard Kalman gain update.
- **Continuous state:** the estimated state itself is not wrapped, so the NMPC always sees a continuous angle.
- **Tuning:** `R` equals the sensor variances; `Q = diag(1e-6, 1e-4, 1e-6, 1e-4)` (larger on the velocities, which are harder to predict); `P0 = 1e-4·I`.

## Closed loop with noise and estimation

`unified_nmpc_with_noise.py` runs the full loop. At every step:

1. The sensor model produces a noisy measurement of the true plant state.
2. The EKF predicts with the previous torque and fuses the measurement.
3. The NMPC solves using the **EKF estimate**, never the true state.
4. The first voltage is applied to the true plant with **one step of actuator latency** (the previous command is the one that acts).

The simulation runs 200 steps (5 s). The plant starts hanging at rest, with the EKF initialized at the true state, and a disturbance of +0.3 rad is injected on φ at step 100 (t = 2.5 s) to test recovery. The plant and the controller's internal model use the same parameters, so there is no model mismatch in this experiment.

## Differences between the MATLAB and Python tracks

The two implementations are not directly interchangeable and should not be compared as before/after versions of the same system.

|                         | MATLAB (`matlab/src/swingup_lqr_controller.m`)                                                                                    | Python (`python/src/`)                                                                                                                                  |
| ----------------------- | --------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Control architecture    | Switched: energy-based swing-up, then LQR once the pendulum is within 40° of upright                                              | Unified: a single NMPC solves swing-up and stabilization together, no switching                                                                         |
| Mechanical parameters   | Illustrative values (`la = 0.25 m`, `lp = 0.15 m`, `mp = 0.05 kg`), rods modeled as uniform, not tied to a platform                | Arm/pendulum dimensions and measured inertias taken from published Quanser QUBE-Servo 2 specifications                                                  |
| Actuator model          | Full electrical model: armature resistance *and* inductance, with motor current as an explicit fifth state                         | Quasi-static voltage-to-torque relationship (instantaneous current, no inductance)                                                                      |
| Motor constants         | `Rm = 2.0 Ω`, `Lm = 5 mH`, `Kt = Ke = 0.05 N·m/A`, `N = 20`                                                                       | `R = 8.4 Ω`, `Kt = Kb = 0.042 N·m/A`, `N = 8`                                                                                                           |
| Voltage limit           | ±24 V (saturation)                                                                                                                | ±4 V (hard constraint inside the optimization)                                                                                                          |
| State constraints       | None (no workspace bound on arm or pendulum angle)                                                                                | Arm angle \|θ\| ≤ π and pendulum angle \|φ\| ≤ 2π, as hard constraints                                                                                   |
| Sensing                 | Full state, including motor current, assumed known                                                                                | Perfect state feedback, or noisy measurements with EKF state estimation                                                                                 |
| Simulation              | Continuous-time (`ode45`), 6 s, controller evaluated inside the ODE                                                               | Discrete-time (RK4, dt = 0.025 s), 5 s, one control update per step                                                                                     |

Because of these differences, the results in `matlab/results/` and `python/results/` come from two different physical systems, not the same plant under two control strategies. They are included together to compare control *architectures*, not to give a strict apples-to-apples numerical comparison.

## Results

### MATLAB (`matlab/results/`)

- **`state_variation_with_time.png`**: the pendulum swings up and settles at the upright position in roughly 2.5 to 3 s, after which the arm returns slowly towards zero. The motor current peaks at about 6 A during the final swing, and the voltage panel shows the commanded voltage in both modes, limited to ±24 V. The pendulum angle settles near −π, which is the same physical position as the dashed π reference line (the state angle is not wrapped).
- **`furuta_motor.mp4`**: 3D animation of the closed-loop response.

### Python (`python/results/`)

- **`nmpc_results_without_noise.png`, `nmpc_closed_loop.mp4`**: NMPC with perfect state feedback. The pendulum swings up and settles at the upright position (φ ≈ π) in roughly one second. The voltage sits at its ±4 V limit during the swing-up and then drops to small corrective commands.
- **`nmpc_results_with_noise.png`**: EKF validation, showing the true state, the noisy measurement and the EKF estimate for all four states. The estimate follows the true state closely, with a small lag at the peaks of the velocity signals. The measured angles appear near +π while the true angle is near −π because measurements are wrapped to [−π, π] and the plant state is continuous; they represent the same physical angle.
- **`nmpc_results_ekf.png`, `nmpc_closed_loop_ekf.mp4`**: the closed-loop trajectory and animation when the NMPC runs on the EKF estimate. The pendulum is held upright and recovers from the +0.3 rad disturbance at t = 2.5 s.

## Known issues

- **MATLAB toolbox dependency**: `wrapToPi` and `deg2rad` come from the Mapping Toolbox in older MATLAB releases.
- **MATLAB switching**: the controller switches hard between swing-up and LQR at 40° from upright, inside the ODE function, which makes the right-hand side discontinuous at the switching instant.

## Current limitations

- All four states, including the velocities, are measured directly in the Python EKF (`H = I`). A real encoder-based setup measures only the angles, so estimating the velocities from angle measurements alone is the natural next step.
- The Python plant and the controller's model share the same parameters, so robustness to parameter mismatch is not tested.
- The only disturbance is a single step on φ, and actuator latency is fixed at one sample.
- The MATLAB controller has not been tested with noise or constraints on the arm angle.
- Everything is simulation; nothing has been run on hardware.

## Requirements

**Python**: see `requirements.txt` (CasADi with the bundled IPOPT solver, NumPy, Matplotlib, and `imageio_ffmpeg`, which Matplotlib uses to write the `.mp4` animations).

```
pip install -r requirements.txt
```

**MATLAB**: base MATLAB with `yline` support (R2018b or later); `wrapToPi` and `deg2rad` require the Mapping Toolbox in older releases.

## Running

**Python**, from `python/src/`:

```
python unified_nmpc_without_noise.py    # NMPC with perfect state feedback
python unified_nmpc_with_noise.py       # NMPC + noisy sensors + EKF
```

Each script runs the closed loop, writes its animation (`.mp4`) and plots (`.png`) to the current directory, and shows the figures. The files in `python/results/` were produced this way.

**MATLAB**: open `matlab/src/swingup_lqr_controller.m` and run it. It prints the controllability check and the LQR gains, plots the state trajectories, and animates the pendulum while writing `furuta_motor.mp4` to the current folder.
