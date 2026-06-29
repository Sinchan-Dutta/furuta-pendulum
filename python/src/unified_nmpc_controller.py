import casadi as ca
import numpy as np
import imageio_ffmpeg
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from dynamics import get_dynamics_function, get_integrator
from motor_model import voltage_to_torque
matplotlib.rcParams['animation.ffmpeg_path'] = imageio_ffmpeg.get_ffmpeg_exe()


params = {
    'la': 0.085,        # arm length (m)
    'lp': 0.129,         # pendulum length (m)
    'mp': 0.024,         # pendulum mass (kg)
    'Jr': 2.306e-4,       # measured arm inertia (replaces J + (1/3)*ma*la^2)
    'Jp': 1.3313e-4,      # measured pendulum inertia (replaces (1/3)*mp*lp^2)
    'g': 9.81,
    'c': 0.0005,          # arm damping (br)
    'b': 0.0001           # pendulum damping (bp)
}

motor_params = {'Kt': 0.042, 'Kb': 0.042, 'R': 8.4, 'N': 8} 

# params = {
#     'la': 0.25,        # arm length (m)
#     'lp': 0.15,         # pendulum length (m)
#     'mp': 0.05,         # pendulum mass (kg)
#     'Jr': 7.25e-3,       # measured arm inertia (replaces J + (1/3)*ma*la^2)
#     'Jp': 3.75e-4,      # measured pendulum inertia (replaces (1/3)*mp*lp^2)
#     'g': 9.81,
#     'c': 0.0005,          # arm damping (br)
#     'b': 0.0005           # pendulum damping (bp)
# }

f = get_dynamics_function(params)
dt = 0.025
F_rk4 = get_integrator(f, dt)

# 2. NMPC setup constants
N = 200                  # horizon length, decided earlier

# 3. Decision variables
X = ca.SX.sym('X', 4, N+1)
U = ca.SX.sym('U', 1, N)

# 4. Parameter: current measured state
x0_param = ca.SX.sym('x0_param', 4)

# 5. Cost function (build up across the horizon)
x_goal = ca.DM([0, 0, ca.pi, 0])
Q_theta = 1.0
Q_theta_dot = 5
Q_phi = 1
Q_phi_dot = 5
R_u = 0.1
R_du = 1

beta = params['Jp']
delta = 0.5 * params['mp'] * params['g'] * params['lp']
Ed = 2 * delta
Q_energy = 1

cost = 0
for k in range(N+1):
  theta_k = X[0, k]
  theta_dot_k = X[1, k]
  phi_k = X[2, k]
  phi_dot_k = X[3, k]

  cost += Q_theta * theta_k**2
  cost += Q_theta_dot * theta_dot_k**2
  cost += Q_phi*(1 + ca.cos(phi_k))          # upright cost, wraparound-safe
  cost += (k/N)**2 * phi_dot_k**2
  
  if k < N:
    cost += R_u * U[0, k]**2

# extra terminal cost
theta_N = X[0, N]
theta_dot_N = X[1, N]
phi_N = X[2, N]
phi_dot_N = X[3, N]

Q_terminal = 1e5   # much larger than the running weights

cost += Q_terminal * (1 + ca.cos(phi_N))

# 6. Constraint list g (dynamics consistency + initial condition)
g = []
for k in range(N):
    theta_dot_k = X[1, k]
    V_k = U[:,k]
    tau_k = voltage_to_torque(V_k, theta_dot_k, motor_params)
    x_next_predicted = F_rk4(X[:, k], tau_k)
    g.append(X[:, k+1] - x_next_predicted)

g.append(X[:, 0] - x0_param)
g = ca.vertcat(*g)
# 7. Bounds: lbx/ubx (variable bounds: phi range, theta range, torque limits)
#    lbg/ubg (constraint bounds: equality constraints => lb=ub=0)
lbg = ca.DM.zeros(g.shape[0])
ubg = ca.DM.zeros(g.shape[0])

theta_min, theta_max = -ca.pi, ca.pi
phi_min, phi_max = -2*ca.pi, 2*ca.pi
theta_dot_min, theta_dot_max = -ca.inf, ca.inf   # no explicit velocity limit (yet)
phi_dot_min, phi_dot_max = -ca.inf, ca.inf

x_lb = [theta_min, theta_dot_min, phi_min, phi_dot_min]
x_ub = [theta_max, theta_dot_max, phi_max, phi_dot_max]

lbx_X = x_lb * (N+1)    # Python list repetition -- repeats the 4-element pattern N+1 times
ubx_X = x_ub * (N+1)

# u_min, u_max = -0.5, 0.5
# lbx_U = [u_min] * N
# ubx_U = [u_max] * N

V_min, V_max = -12, 12
lbx_U = [V_min] * N
ubx_U = [V_max] * N

lbx = lbx_X + lbx_U
ubx = ubx_X + ubx_U
# 8. Package into nlpsol, call solver
w = ca.vertcat(ca.reshape(X, -1, 1), ca.reshape(U, -1, 1))
nlp = {'x': w, 'f': cost, 'g': g, 'p': x0_param}
opts = {
    'ipopt.print_level': 0,
    'print_time': 0
}
solver = ca.nlpsol('solver', 'ipopt', nlp, opts)

x_current = [0, 0, 0.01, 0]   # hanging down, at rest

t_array = np.arange(N+1) * dt   # time at each of the N+1 horizon steps

# crude growing-amplitude oscillation for phi, peaking near pi by partway through
phi_guess_array = 1.5 * (1 - np.cos(2*np.pi*0.4*t_array)) * np.minimum(t_array / 2.0, 1.0)

theta_guess_array = np.zeros(N+1)
theta_dot_guess_array = np.zeros(N+1)
phi_dot_guess_array = np.gradient(phi_guess_array, dt)   # rough derivative of the phi guess, for consistency

X_guess = []
for k in range(N+1):
    X_guess += [theta_guess_array[k], theta_dot_guess_array[k], phi_guess_array[k], phi_dot_guess_array[k]]

U_guess = [0.0] * N

w_guess = X_guess + U_guess

sol = solver(x0=w_guess, lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg, p=x_current)

w_opt = sol['x']

n_X = 4*(N+1)
X_opt_flat = w_opt[0:n_X]
U_opt_flat = w_opt[n_X:]

X_opt = ca.reshape(X_opt_flat, 4, N+1)   # should undo the packing exactly
U_opt = ca.reshape(U_opt_flat, 1, N)

theta_opt = X_opt[0, :].full().flatten()
phi_opt = X_opt[2, :].full().flatten()
u_opt = U_opt[0, :].full().flatten()

# Plotting

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.set_xlim(-0.4, 0.4)
ax.set_ylim(-0.4, 0.4)
ax.set_zlim(0, 0.6)

arm1_line, = ax.plot([], [], [], 'b', linewidth=3)
arm2_line, = ax.plot([], [], [], 'r', linewidth=3)

la = params['la']
lp = params['lp']
h = 0.3
fps = 1/dt*2
t_anim = np.arange(0, t_array[-1], 1/fps)
theta_anim = np.interp(t_anim, t_array, theta_opt)
phi_anim = np.interp(t_anim, t_array, phi_opt)
u_anim = np.interp(t_anim, t_array[0:N], u_opt)


def update(frame_idx):
  theta = theta_anim[frame_idx]
  phi = phi_anim[frame_idx]

  xb = la*np.cos(theta)
  yb = la*np.sin(theta)
  zb = h

  xp = xb - lp*np.sin(phi)*np.sin(theta)
  yp = yb + lp*np.sin(phi)*np.cos(theta)
  zp = h - lp*np.cos(phi)

  arm1_line.set_data_3d([0, xb], [0, yb], [h, zb])
  arm2_line.set_data_3d([xb, xp], [yb, yp], [zb, zp])

  return arm1_line, arm2_line

anim = FuncAnimation(fig, update, frames=len(t_anim), interval=1000/fps)
anim.save('inverted_pendulum_stabilization.mp4', fps=fps)
plt.show()

fig, axs = plt.subplots(3, 1, figsize=(4, 10))

axs[0].plot(t_anim, phi_anim)
axs[0].set_title("Phi Opt")  # Optional: Adds clarity to the first plot

axs[1].plot(t_anim, theta_anim)
axs[1].set_title("Theta Opt")  # Optional: Adds clarity to the second plot
axs[2].plot(t_anim, u_anim)
axs[2].set_title("Voltage")

plt.tight_layout()
plt.show()