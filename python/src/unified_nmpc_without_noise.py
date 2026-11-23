"""
Closed-loop NMPC for Furuta Pendulum
- Solves NMPC at every timestep over a short horizon (N=20)
- Applies first control action
- Integrates real plant forward by dt
- Repeats for 200 steps (~5 seconds of real time at dt=0.025)
"""

import casadi as ca
import numpy as np
import imageio_ffmpeg
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from dynamics import get_dynamics_function, get_integrator
from motor_model import voltage_to_torque

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
# 2. NMPC SETUP — SHORT HORIZON (N=20, not 200)
# ═══════════════════════════════════════════════════════════════════
N = 20                    # horizon length (0.5s lookahead)
N_SIM = 200               # total closed-loop simulation steps (5s)

# Decision variables
X = ca.SX.sym('X', 4, N + 1)
U = ca.SX.sym('U', 1, N)
x0_param = ca.SX.sym('x0_param', 4)

# Cost
Q_theta, Q_theta_dot = 40.0, 20.0
Q_phi, Q_phi_dot = 1.0, 5.0
R_u = 0.1
Q_terminal = 1e5

cost = 0
for k in range(N + 1):
    theta_k     = X[0, k]
    theta_dot_k = X[1, k]
    phi_k       = X[2, k]
    phi_dot_k   = X[3, k]

    cost += Q_theta * theta_k**2
    cost += Q_theta_dot * theta_dot_k**2
    cost += Q_phi * (1 + ca.cos(phi_k))         # upright cost
    cost += (k / N)**2 * Q_phi_dot * phi_dot_k**2

    if k < N:
        cost += R_u * U[0, k]**2

# Terminal cost
cost += Q_terminal * (1 + ca.cos(X[2, N]))

# Dynamics constraints
g = []
for k in range(N):
    theta_dot_k = X[1, k]
    V_k = U[:, k]
    tau_k = voltage_to_torque(V_k, theta_dot_k, motor_params)
    x_next_predicted = F_rk4(X[:, k], tau_k)
    g.append(X[:, k + 1] - x_next_predicted)

g.append(X[:, 0] - x0_param)          # initial condition constraint
g = ca.vertcat(*g)

# Bounds
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
# 3. BUILD SOLVER (ONCE, reused every timestep)
# ═══════════════════════════════════════════════════════════════════
w = ca.vertcat(ca.reshape(X, -1, 1), ca.reshape(U, -1, 1))
nlp = {'x': w, 'f': cost, 'g': g, 'p': x0_param}

# KEY: use 'warm_start_init_point' + limited iterations for speed
opts = {
    'ipopt.print_level': 0,
    'print_time': 0,
    'ipopt.max_iter': 50,               # cap iterations for speed
    'ipopt.tol': 1e-4,                  # relax tolerance slightly
    'ipopt.warm_start_init_point': 'yes',
    'ipopt.mu_init': 1e-3,              # helps warm-starting
}
solver = ca.nlpsol('solver', 'ipopt', nlp, opts)

# ═══════════════════════════════════════════════════════════════════
# 4. INITIAL GUESS (reused and shifted across timesteps)
# ═══════════════════════════════════════════════════════════════════
t_horizon = np.arange(N + 1) * dt

# Initial guess: swing-up trajectory
phi_guess = 1.5 * (1 - np.cos(2 * np.pi * 0.4 * t_horizon)) \
            * np.minimum(t_horizon / 2.0, 1.0)
phi_dot_guess = np.gradient(phi_guess, dt)

X_guess = []
for k in range(N + 1):
    X_guess += [0.0, 0.0, phi_guess[k], phi_dot_guess[k]]
U_guess = [0.0] * N
w_guess = ca.DM(X_guess + U_guess)          # persistent warm-start vector

# ═══════════════════════════════════════════════════════════════════
# 5. CLOSED-LOOP SIMULATION LOOP  ← THIS IS THE LOOP YOU ASKED ABOUT
# ═══════════════════════════════════════════════════════════════════
x_current = np.array([0.0, 0.0, 0.01, 0.0])   # start hanging, at rest

# Logs
X_log = np.zeros((N_SIM + 1, 4))
U_log = np.zeros(N_SIM)
X_log[0] = x_current

print("Starting closed-loop NMPC...")
for step in range(N_SIM):
    # ── (a) Solve NMPC with current state ──────────────────────────
    sol = solver(
        x0=w_guess,
        lbx=lbx, ubx=ubx,
        lbg=lbg, ubg=ubg,
        p=x_current
    )
    w_opt = sol['x']
    w_guess = w_opt                              # warm-start next solve

    # ── (b) Extract first control action ───────────────────────────
    n_X = 4 * (N + 1)
    U_opt = ca.reshape(w_opt[n_X:], 1, N)
    u0 = float(U_opt[0, 0])

    # ── (c) Apply to REAL plant (use tau from motor model) ─────────
    theta_dot_now = x_current[1]
    tau_applied = voltage_to_torque(ca.DM([u0]), theta_dot_now, motor_params)
    x_next = F_rk4(x_current, tau_applied).full().flatten()

    # ── (d) Log and advance ────────────────────────────────────────
    U_log[step] = u0
    X_log[step + 1] = x_next
    x_current = x_next

    if step % 20 == 0:
        phi_now = x_current[2]
        print(f"  step {step:3d}/{N_SIM} | u={u0:6.3f} V | phi={phi_now:.3f} rad")

print("Closed-loop simulation done.")

# ═══════════════════════════════════════════════════════════════════
# 6. ANIMATION
# ═══════════════════════════════════════════════════════════════════
t_full = np.arange(N_SIM + 1) * dt
la, lp, h = params['la'], params['lp'], 0.3
fps = 30
t_anim = np.arange(0, t_full[-1], 1 / fps)

theta_anim = np.interp(t_anim, t_full, X_log[:, 0])
phi_anim   = np.interp(t_anim, t_full, X_log[:, 2])

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_xlim(-0.4, 0.4); ax.set_ylim(-0.4, 0.4); ax.set_zlim(0, 0.6)
ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
ax.set_title('Closed-Loop NMPC — Furuta Pendulum')

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
anim.save('nmpc_closed_loop.mp4', fps=fps)
plt.show()

# ═══════════════════════════════════════════════════════════════════
# 7. PLOTS
# ═══════════════════════════════════════════════════════════════════
fig, axs = plt.subplots(4, 1, figsize=(7, 10), sharex=True)
axs[0].plot(t_full, X_log[:, 0]); axs[0].set_ylabel(r'$\theta$ (rad)')
axs[1].plot(t_full, X_log[:, 2]); axs[1].set_ylabel(r'$\phi$ (rad)')
axs[2].plot(t_full, X_log[:, 3]); axs[2].set_ylabel(r'$\dot\phi$ (rad/s)')
axs[3].plot(t_full[:-1], U_log);  axs[3].set_ylabel('Voltage (V)')
axs[3].set_xlabel('Time (s)')
for a in axs: a.grid(True)
plt.tight_layout()
plt.savefig('nmpc_results_without_noise.png', dpi=150)
plt.show()