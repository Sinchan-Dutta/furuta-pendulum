# Furuta Pendulum: From Linear LQR to Constrained Nonlinear MPC

This repository documents the development of a control stack for a Furuta Pendulum (rotary inverted pendulum), progressing from a classical linearized LQR/energy-based swing-up baseline to a unified Nonlinear Model Predictive Controller (NMPC), built with CasADi and IPOPT.

The project was developed as part of a master's internship application, with an emphasis on physically realistic modeling (literature-sourced hardware parameters, a DC motor actuator model, and explicit state/input constraints) rather than a purely textbook simulation.

**Note on the two implementations**: the MATLAB and Python controllers were developed somewhat independently rather than as a single linear progression. The Python NMPC line (`python/`) is the primary, original line of development for this project, including the dynamics derivation, parameter sourcing, and constrained optimal control formulation described below. The MATLAB controller (`matlab/`) was developed as a separate reference exploration of a more classical architecture — a switched energy-based swing-up and LQR stabilizer, with the actuator modeled in more electrical detail (including the motor's armature inductance and current as an explicit state, rather than the quasi-static voltage-to-torque relationship used in the Python model). It is included here as a useful point of comparison against a more classical, switched control architecture, rather than as an earlier version of the Python implementation. The two use related but not identical physical and actuator parameters (see "Differences between the MATLAB and Python tracks" below).

## Project structure

```
furuta-pendulum/
├── docs/
│   └── pendulum_dynamics.pdf       # Derivation of the nonlinear equations of motion
├── matlab/
│   ├── src/
│   │   └── swingup_lqr_controller.m
│   └── results/
│       ├── furuta_motor.mp4
│       └── state_variation_with_time.png
├── python/
│   ├── src/
│   │   ├── dynamics.py             # CasADi symbolic nonlinear EOM + RK4 integrator
│   │   ├── motor_model.py          # DC motor voltage-to-torque model
│   │   └── unified_nmpc_controller.py
│   └── results/
│       ├── Figure_1.png
│       └── inverted_pendulum_stabilization.mp4
├── requirements.txt
└── README.md
```

## Background

A Furuta Pendulum consists of a horizontal arm, driven by a single motor at its base, with an unactuated pendulum rod attached to the far end of the arm. The control objective is to swing the pendulum up from its hanging-down rest position and stabilize it in the inverted (upright) position, using only the single arm-side actuator. Since the pendulum has no actuator of its own, it can only be driven indirectly, through the dynamic coupling between the arm's rotation and the pendulum's swing.

This system is a standard benchmark in nonlinear and underactuated control because it captures, in a small and analyzable package, many challenges found in more complex underactuated systems: it requires a large-angle, genuinely nonlinear swing-up maneuver, followed by precision stabilization at an unstable equilibrium, all while respecting real actuator and workspace limits.

## Physical model

The arm angle is denoted θ and the pendulum angle φ, with φ = π corresponding to the upright (inverted) equilibrium. The nonlinear equations of motion are derived from first principles via the Euler–Lagrange method in `docs/pendulum_dynamics.pdf` — starting from the kinetic and potential energy of the centre rod, arm, and pendulum, through the Lagrangian and its partial derivatives, to the final coupled nonlinear ODEs, their linearization about both the hanging-down and inverted equilibria, and the corresponding transfer functions. The same document also derives the actuator model in both forms used across this repository: the full electrical model with motor inductance (used in the MATLAB controller), and the quasi-static simplification obtained by assuming the motor's electrical time constant is much faster than the control sample rate, i.e. `Lm/Rm << sample time` (used in the Python controller). The resulting nonlinear equations of motion take the form:

```
Den      = α·β + β²·sin²(φ) − (γ·cos(φ))²
U        = −β·sin(2φ)·θ̇·φ̇ + γ·sin(φ)·φ̇² − c·θ̇ + τ
V        = 0.5·β·sin(2φ)·θ̇² − δ·sin(φ) − b·φ̇
θ̈        = (β·U − γ·cos(φ)·V) / Den
φ̈        = (−γ·cos(φ)·U + (α + β·sin²(φ))·V) / Den
```

where α, β, γ, δ are inertia/coupling coefficients built from the physical parameters (link lengths, masses, inertias), c and b are arm and pendulum viscous damping coefficients, and τ is the torque applied at the arm shaft.

Physical parameters were sourced from published specifications for the Quanser QUBE-Servo 2 rotary pendulum platform (arm length, pendulum length and mass, measured arm/pendulum inertia, and damping coefficients), rather than chosen arbitrarily, so that the resulting dynamics and actuator limits correspond to a real, physically grounded system.

## Actuator model

Torque is not a free input. The arm is driven by a DC gearmotor, and the full electrical dynamics are:

```
Lm·(di/dt) = Vin − Rm·i − Kb·N·θ̇
τ          = N·Kt·i
```

where `Kt` and `Kb` are the motor's torque and back-EMF constants, `Rm` and `Lm` are the armature resistance and inductance, and `N` is the gear ratio. Assuming the motor's electrical time constant `Lm/Rm` is much faster than the control sample rate (`Lm·di/dt ≈ 0`), this reduces to a quasi-static voltage-to-torque relationship:

```
τ = N·Kt·(Vin − Kb·N·θ̇) / Rm
```

Either way, achievable torque depends on the arm's current speed (back-EMF reduces available torque as speed increases), and the control input that is actually bounded is voltage, not torque directly. `python/src/motor_model.py` implements the quasi-static form; the MATLAB controller instead carries the full electrical dynamics as a fifth system state, capturing the motor current's transient response to a voltage command.

## MATLAB baseline: swing-up + LQR

`matlab/src/swingup_lqr_controller.m` implements a classical two-mode controller:

- **Energy-based swing-up**: while the pendulum is far from upright, torque is computed from an energy-shaping law that pumps mechanical energy into the pendulum until it approaches the target energy needed to reach the top, combined with a centering term that keeps the arm near its workspace center.
- **LQR stabilization**: once the pendulum is within a threshold of upright, control switches to a linear-quadratic regulator, with gains computed from a linearization of the dynamics about the upright equilibrium (solved via the Hamiltonian eigenvector method rather than a black-box LQR call, to make the underlying Riccati solution explicit).

This baseline includes the full motor electrical model (resistance, inductance, back-EMF) as part of the state, with the swing-up law's desired torque inverted through the motor equations to produce a feed-forward voltage command, and the LQR gain computed over the augmented five-state system (arm angle, arm velocity, pendulum angle error, pendulum velocity, motor current).

## Python: unified Nonlinear MPC

`python/src/` contains a from-scratch CasADi/IPOPT implementation that replaces the switched swing-up/LQR architecture with a single nonlinear optimal control problem, solved at every control step over a finite receding horizon.

- **`dynamics.py`**: builds the nonlinear equations of motion as a CasADi symbolic function, and a fourth-order Runge–Kutta integrator built on top of it. Working symbolically (rather than purely numerically) allows the same dynamics model to be reused, unmodified, inside the NMPC's internal prediction model.
- **`motor_model.py`**: the voltage-to-torque relationship described above, used to convert the NMPC's voltage decision variable into the torque value the dynamics model expects, so that the optimizer is solving for a realistic control input (bounded by ±12V) rather than an unconstrained or unrealistically bounded torque.
- **`unified_nmpc_controller.py`**: formulates and solves the NMPC problem. At each step, the controller optimizes a sequence of voltages over a multi-second prediction horizon to minimize a cost combining: arm and pendulum tracking error (with the pendulum's upright cost expressed as `1 + cos(φ)` to handle angle wraparound correctly), a pendulum energy-shaping term encouraging efficient energy buildup during swing-up, control effort and control rate-of-change penalties (encouraging smooth, realistic voltage commands rather than bang-bang switching), and a heavily weighted terminal cost enforcing that the pendulum is upright and at rest at the end of the horizon (without this, the optimizer has no incentive to remain balanced near the horizon's end, since the cost does not extend past it). The optimization is subject to the true nonlinear dynamics as an equality constraint at every step, together with arm workspace and pendulum rotation-range bounds and the actuator voltage limit.

Unlike the MATLAB baseline, this controller does not switch between separate swing-up and stabilization laws; a single optimization discovers both the energy-pumping swing-up maneuver and the final stabilizing approach as one continuous trajectory, subject to the same realistic constraints throughout.

## Differences between the MATLAB and Python tracks

The two implementations are not directly interchangeable and should not be compared as before/after versions of the same system. Key differences:

| | MATLAB (`matlab/src/swingup_lqr_controller.m`) | Python (`python/src/`) |
|---|---|---|
| Control architecture | Switched: energy-shaping swing-up, then LQR stabilization once near upright | Unified: a single NMPC solves swing-up and stabilization together, no switching |
| Mechanical parameters | Earlier, less rigorously sourced arm/pendulum dimensions and inertias | Arm/pendulum dimensions and inertias taken directly from published Quanser QUBE-Servo 2 specifications |
| Actuator model | Full electrical model: armature resistance *and* inductance, with motor current as an explicit fifth state, giving the motor a genuine electrical transient response to a voltage command | Quasi-static voltage-to-torque relationship (instantaneous current, no inductance), used as the link between the NMPC's voltage decision variable and the torque the mechanical dynamics see |
| Motor constants | `R_m = 2.0 Ω`, `K_t = K_e = 0.05 N·m/A`, `N = 20` | `R = 8.4 Ω`, `Kt = Kb = 0.042 N·m/A`, `N = 8`, also Quanser-derived |
| State/input constraints | Voltage saturation only (±24V); no explicit workspace bound on arm or pendulum angle | Explicit arm angle and pendulum rotation-range bounds, plus voltage saturation (±12V), all enforced as hard constraints inside the optimization |

Because of these differences, the two animations and result plots in `matlab/results/` and `python/results/` represent two different physical systems and actuator models, not the same plant under two different control strategies. They are included together to compare control *architectures* — a classical switched linear/energy-based design against a constrained nonlinear optimization-based design — rather than to provide a strict apples-to-apples numerical comparison.

## Results

`matlab/results/` and `python/results/` contain example state trajectories and rendered 3D animations of the closed-loop response for each controller, generated from the same physical model for comparison. The NMPC result shows a smooth multi-swing energy buildup, a decisive final swing to upright, and a stable hold, while respecting the actuator voltage bound throughout.

## Requirements

See `requirements.txt`. The Python implementation requires CasADi (bundled with the IPOPT solver) and standard scientific Python packages (NumPy, Matplotlib) for simulation and visualization.
