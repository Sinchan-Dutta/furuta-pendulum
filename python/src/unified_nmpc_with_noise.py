"""
Closed-loop NMPC for Furuta Pendulum with EKF State Estimation
- True plant: nonlinear dynamics + RK4 integration
- Sensor: adds Gaussian noise to true state
- EKF: fuses noisy measurements into filtered state estimate
- NMPC: uses EKF estimate (not truth) — like a real robot
"""

import casadi as ca
import numpy as np
import imageio_ffmpeg
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from dynamics import get_dynamics_function, get_integrator
from motor_model import voltage_to_torque

# ═══════════════════════════════════════════════════════════════════
# NEW IMPORTS — make sure these files exist in your project
# ═══════════════════════════════════════════════════════════════════
from noise_models import SensorModel
from ekf import EKF

matplotlib.rcParams['animation.ffmpeg_path'] = imageio_ffmpeg.get_ffmpeg_exe()

# ═══════════════════════════════════════════════════════════════════
# 1. PARAMETERS
# ═══════════════════════════════════════════════════════════════════
params = {
    'la': 0.085, 'lp': 0.129, 'mp': 0.024,
    'Jr': 2.306e-4, 'Jp': 1.3313e-4,
    'g': 9.81, 'c': 0.0005, 'b': 0.0001
}
motor_params = {'Kt': 0.042, 'Kb': 0.042, 'R': 8.4, 'N': 8}

f = get_dynamics_function(params)
dt = 0.025
F_rk4 = get_integrator(f, dt)

# ═══════════════════════════════════════════════════════════════════
# 2. NMPC SETUP — SHORT HORIZON (N=20)
# ═══════════════════════════════════════════════════════════════════
N = 20
N_SIM = 200

X = ca.SX.sym('X', 4, N + 1)
U = ca.SX.sym('U', 1, N)
x0_param = ca.SX.sym('x0_param', 4)

Q_theta, Q_theta_dot = 50.0, 15.0
Q_phi, Q_phi_dot = 1.0, 5.0
R_u = 0.5
Q_terminal = 1e5

cost = 0
for k in range(N + 1):
    theta_k     = X[0, k]
    theta_dot_k = X[1, k]
    phi_k       = X[2, k]
    phi_dot_k   = X[3, k]

    cost += Q_theta * theta_k**2
    cost += Q_theta_dot * theta_dot_k**2
    cost += Q_phi * (1 + ca.cos(phi_k))
    cost += (k / N)**2 * Q_phi_dot * phi_dot_k**2

    if k < N:
        cost += R_u * U[0, k]**2

cost += Q_terminal * (1 + ca.cos(X[2, N]))

# Dynamics constraints (same as before)
g = []
for k in range(N):
    theta_dot_k = X[1, k]
    V_k = U[:, k]
    tau_k = voltage_to_torque(V_k, theta_dot_k, motor_params)
    x_next_predicted = F_rk4(X[:, k], tau_k)
    g.append(X[:, k + 1] - x_next_predicted)

g.append(X[:, 0] - x0_param)
g = ca.vertcat(*g)

lbg = ca.DM.zeros(g.shape[0])
ubg = ca.DM.zeros(g.shape[0])

x_lb = [-ca.pi, -ca.inf, -2 * ca.pi, -ca.inf]
x_ub = [ ca.pi,  ca.inf,  2 * ca.pi,  ca.inf]
lbx_X = x_lb * (N + 1)
ubx_X = x_ub * (N + 1)

V_min, V_max = -4, 4
lbx_U = [V_min] * N
ubx_U = [V_max] * N

lbx = lbx_X + lbx_U
ubx = ubx_X + ubx_U

# ═══════════════════════════════════════════════════════════════════
# 3. BUILD SOLVER
# ═══════════════════════════════════════════════════════════════════
w = ca.vertcat(ca.reshape(X, -1, 1), ca.reshape(U, -1, 1))
nlp = {'x': w, 'f': cost, 'g': g, 'p': x0_param}

opts = {
    'ipopt.print_level': 0,
    'print_time': 0,
    'ipopt.max_iter': 50,
    'ipopt.tol': 1e-4,
    'ipopt.warm_start_init_point': 'yes',
    'ipopt.mu_init': 1e-3,
}
solver = ca.nlpsol('solver', 'ipopt', nlp, opts)

# ═══════════════════════════════════════════════════════════════════
# 4. INITIAL GUESS
# ═══════════════════════════════════════════════════════════════════
t_horizon = np.arange(N + 1) * dt

w_guess = np.zeros(w.shape[0])

# ═══════════════════════════════════════════════════════════════════
# 5. NEW: SETUP SENSOR + EKF
# ═══════════════════════════════════════════════════════════════════

# --- Sensor noise model ---
sensor = SensorModel(
    sigma_theta=0.002,       # arm angle noise (rad)
    sigma_phi=0.003,         # pendulum angle noise (rad)
    sigma_theta_dot=0.05,    # arm velocity noise (rad/s)
    sigma_phi_dot=0.08       # pendulum velocity noise (rad/s)
)

# --- EKF tuning ---
# Q: how much do you trust your model? Larger = trust measurements more.
# Process noise is larger on velocity terms (they are harder to predict).
Q_ekf = np.diag([1e-6, 1e-4, 1e-6, 1e-4])

# R: measurement noise covariance — should equal sensor variance (honest)
R_ekf = np.diag([
    sensor.sigma[0]**2,
    sensor.sigma[1]**2,
    sensor.sigma[2]**2,
    sensor.sigma[3]**2
])

# P0: initial uncertainty
P0 = np.eye(4) * 1e-4

# Initial state
x_true = np.array([0.0, 0.0, 0.0, 0.0])

# --- Create EKF instance ---
# NOTE: your EKF class must accept f, F_rk4, params, motor_params.
# If the signature differs, adapt the constructor call accordingly.
ekf = EKF(
    f_dyn=f,
    F_rk4=F_rk4,
    dt=dt,
    params=params,
    motor_params=motor_params,
    Q_ekf=Q_ekf,
    R_ekf=R_ekf,
    x0_est=x_true.copy(),
    P0=P0
)

# ═══════════════════════════════════════════════════════════════════
# 6. CLOSED-LOOP SIMULATION LOOP
# ═══════════════════════════════════════════════════════════════════

# Logs — we now log FOUR signals
X_true_log = np.zeros((N_SIM + 1, 4))   # true plant state
X_meas_log = np.zeros((N_SIM + 1, 4))   # raw noisy measurement
X_hat_log  = np.zeros((N_SIM + 1, 4))   # EKF estimate
U_log      = np.zeros(N_SIM)            # commanded voltage

X_true_log[0] = x_true
X_meas_log[0] = sensor.measure(x_true)
X_hat_log[0]  = ekf.x_hat.copy()

u_prev = 0.0    # for 1-step actuator latency

print("Starting closed-loop NMPC with EKF...")
for step in range(N_SIM):

    # ── (a) MEASURE: add noise to true state ───────────────────────
    z = sensor.measure(x_true)
    X_meas_log[step] = z

    # ── (b) EKF: predict + update using PREVIOUS applied torque ────
    tau_for_pred = voltage_to_torque(
        ca.DM([u_prev]),
        float(ekf.x_hat[1]),
        motor_params
    )
    x_hat = ekf.step(tau_for_pred, z)
    X_hat_log[step] = x_hat

    # ── (c) NMPC: solve using EKF ESTIMATE, not truth ──────────────
    sol = solver(
        x0=w_guess,
        lbx=lbx, ubx=ubx,
        lbg=lbg, ubg=ubg,
        p=x_hat                    # ← EKF estimate in place of x_true
    )
    
    w_opt = sol['x']
    w_guess = w_opt

    # Extract first control action
    n_X = 4 * (N + 1)
    U_opt = ca.reshape(w_opt[n_X:], 1, N)
    u0 = float(U_opt[0, 0])

    # ── (d) APPLY to TRUE plant (with 1-step latency) ──────────────
    tau_applied = voltage_to_torque(
        ca.DM([u_prev]),           # latency: use previous command
        float(x_true[1]),
        motor_params
    )
    x_true = np.array(F_rk4(x_true, tau_applied)).flatten()

    # Optional disturbance
    if step == 100:
        x_true[2] += 0.3
        print("  >>> DISTURBANCE <<<")

    # ── (e) Log ────────────────────────────────────────────────────
    U_log[step] = u0
    X_true_log[step + 1] = x_true
    X_meas_log[step + 1] = sensor.measure(x_true)
    X_hat_log[step + 1]  = ekf.x_hat.copy()

    u_prev = u0

    if step % 20 == 0:
        print(f"  step {step:3d}/{N_SIM} | u={u0:6.3f} V | "
              f"phi_true={x_true[2]:.3f} | phi_hat={x_hat[2]:.3f}")

print("Closed-loop simulation done.")

# ═══════════════════════════════════════════════════════════════════
# 7. ANIMATION (uses TRUE state for rendering)
# ═══════════════════════════════════════════════════════════════════
t_full = np.arange(N_SIM + 1) * dt
la, lp, h = params['la'], params['lp'], 0.3
fps = 30
t_anim = np.arange(0, t_full[-1], 1 / fps)

theta_anim = np.interp(t_anim, t_full, X_true_log[:, 0])
phi_anim   = np.interp(t_anim, t_full, X_true_log[:, 2])

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_xlim(-0.4, 0.4); ax.set_ylim(-0.4, 0.4); ax.set_zlim(0, 0.6)
ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
ax.set_title('Closed-Loop NMPC + EKF — Furuta Pendulum')

arm1_line, = ax.plot([], [], [], 'b', linewidth=3)
arm2_line, = ax.plot([], [], [], 'r', linewidth=3)

def update(i):
    theta, phi = theta_anim[i], phi_anim[i]
    xb, yb, zb = la * np.cos(theta), la * np.sin(theta), h
    xp = xb - lp * np.sin(phi) * np.sin(theta)
    yp = yb + lp * np.sin(phi) * np.cos(theta)
    zp = h - lp * np.cos(phi)
    arm1_line.set_data_3d([0, xb], [0, yb], [h, zb])
    arm2_line.set_data_3d([xb, xp], [yb, yp], [zb, zp])
    return arm1_line, arm2_line

anim = FuncAnimation(fig, update, frames=len(t_anim), interval=1000 / fps)
anim.save('nmpc_closed_loop_ekf.mp4', fps=fps)
plt.show()

# ═══════════════════════════════════════════════════════════════════
# 8. CONTROL PLOTS (true vs. commanded)
# ═══════════════════════════════════════════════════════════════════
fig, axs = plt.subplots(4, 1, figsize=(7, 10), sharex=True)
axs[0].plot(t_full, X_true_log[:, 0]); axs[0].set_ylabel(r'$\theta$ (rad)')
axs[1].plot(t_full, X_true_log[:, 2]); axs[1].set_ylabel(r'$\phi$ (rad)')
axs[2].plot(t_full, X_true_log[:, 3]); axs[2].set_ylabel(r'$\dot\phi$ (rad/s)')
axs[3].plot(t_full[:-1], U_log);       axs[3].set_ylabel('Voltage (V)')
axs[3].set_xlabel('Time (s)')
for a in axs: a.grid(True)
plt.tight_layout()
plt.savefig('nmpc_results_ekf.png', dpi=150)
plt.show()

# ═══════════════════════════════════════════════════════════════════
# 9. EKF VALIDATION PLOT — the money figure
# ═══════════════════════════════════════════════════════════════════
fig, axs = plt.subplots(4, 1, figsize=(10, 10), sharex=True)

labels = [r'$\theta$ (rad)', r'$\dot\theta$ (rad/s)',
          r'$\phi$ (rad)',   r'$\dot\phi$ (rad/s)']

for i in range(4):
    axs[i].plot(t_full, X_true_log[:, i], 'k-', lw=2, label='True')
    axs[i].plot(t_full, X_meas_log[:, i], 'r.', ms=2, alpha=0.3,
                label='Measured (noisy)')
    axs[i].plot(t_full, X_hat_log[:, i],  'b--', lw=1.5,
                label='EKF estimate')
    axs[i].set_ylabel(labels[i])
    axs[i].grid(True)
    axs[i].legend(loc='upper right', fontsize=8)

axs[-1].set_xlabel('Time (s)')
fig.suptitle('EKF State Estimation with Noisy Sensors', fontweight='bold')
plt.tight_layout()
plt.savefig('nmpc_results_with_noise.png', dpi=150)
plt.show()
