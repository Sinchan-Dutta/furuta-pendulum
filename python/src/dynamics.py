import casadi as ca

def get_dynamics_function(params):
  la = params['la']
  lp = params['lp']
  mp = params['mp']
  Jr = params['Jr']
  Jp = params['Jp']
  g = params['g']
  c = params['c']
  b = params['b']

  alpha = Jr + mp*la**2
  beta  = Jp
  gamma = (1/2)*mp*la*lp
  delta = (1/2)*mp*g*lp

  x = ca.SX.sym('x', 4)
  u = ca.SX.sym('u', 1)

  theta = x[0]
  theta_dot = x[1]
  phi = x[2]
  phi_dot = x[3]

  Den = alpha*beta + beta**2*ca.sin(phi)**2 - (gamma*ca.cos(phi))**2
  
  U = -beta*ca.sin(2*phi)*theta_dot*phi_dot + gamma*ca.sin(phi)*phi_dot**2 - c*theta_dot + u

  V = 0.5*beta*ca.sin(2*phi)*theta_dot**2 - delta*ca.sin(phi) - b*phi_dot

  theta_dd = (beta*U - gamma*ca.cos(phi)*V)/Den;
  phi_dd   = (-gamma*ca.cos(phi)*U + (alpha + beta*ca.sin(phi)**2)*V)/Den;

  dx = ca.vertcat(theta_dot, theta_dd, phi_dot, phi_dd)
  f = ca.Function('f', [x, u], [dx])
  return f


def get_integrator(f, dt):
  x = ca.SX.sym('x', 4)
  u = ca.SX.sym('u', 1)

  k1 = f(x, u)
  k2 = f(x + dt/2 * k1, u)
  k3 = f(x + dt/2 * k2, u)
  k4 = f(x + dt * k3, u)

  x_next = x + dt/6 * (k1 + 2*k2 + 2*k3 + k4)

  F_rk4 = ca.Function('F_rk4', [x, u], [x_next])
  return F_rk4


if __name__ == "__main__":
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

  f = get_dynamics_function(params)
  print(f([0.2, 0.5, 0.3, -0.4], 1.5))

  F_rk4 = get_integrator(f, 0.0001)   # dt = 10ms, a reasonable control sample rate
  print(F_rk4([0.2, 0.5, 0.3, -0.4], 1.5))

  x0 = ca.DM([0.2, 0.5, 0.3, -0.4])
  x_next_euler = x0 + f(x0, 1.5) * 0.0001
  print(x_next_euler)
