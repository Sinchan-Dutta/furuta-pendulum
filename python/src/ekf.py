import numpy as np
import casadi as ca

class EKF:
    """EKF for Furuta pendulum. State: [θ, θ̇, φ, φ̇]."""
    
    def __init__(self, f_dyn, F_rk4, dt, params, motor_params,
                 Q_ekf, R_ekf, x0_est, P0):
        """
        f_dyn       : nonlinear dynamics function
        F_rk4       : RK4 integrator (for prediction)
        dt          : timestep
        Q_ekf       : 4x4 process noise covariance
        R_ekf       : 4x4 measurement noise covariance
        x0_est      : initial state estimate
        P0          : initial covariance
        """
        self.f = f_dyn
        self.F_rk4 = F_rk4
        self.dt = dt
        self.params = params
        self.motor_params = motor_params
        
        self.x_hat = np.array(x0_est, dtype=float)
        self.P = np.array(P0, dtype=float)
        self.Q = np.array(Q_ekf, dtype=float)
        self.R = np.array(R_ekf, dtype=float)
        
        # Build Jacobian function (4x4) symbolically with CasADi
        x_sym = ca.SX.sym('x', 4)
        tau_sym = ca.SX.sym('tau', 1)
        f_sym = self.f(x_sym, tau_sym)
        A_sym = ca.jacobian(f_sym, x_sym)
        self.A_func = ca.Function('A', [x_sym, tau_sym], [A_sym])
    
    def predict(self, tau_applied):
        """Propagate estimate forward by dt using RK4."""
        # Nonlinear prediction
        self.x_hat = np.array(self.F_rk4(self.x_hat, tau_applied)).flatten()
        
        # Jacobian at current estimate
        A = np.array(self.A_func(self.x_hat, tau_applied)).astype(float)
        # Discrete-time: F ≈ I + A*dt (first-order)
        F = np.eye(4) + A * self.dt
        
        self.P = F @ self.P @ F.T + self.Q
        # NOTE: State wrapping removed here to ensure continuous state for NMPC

    def update(self, z):
        """Fuse noisy measurement z."""
        # Measurement model: h(x) = x (we measure all 4 states)
        H = np.eye(4)
        
        # Innovation calculation
        y = z - self.x_hat
        
        # Wrap the angular ERROR to take the shortest path [-π, π]
        y[0] = np.arctan2(np.sin(y[0]), np.cos(y[0]))
        y[2] = np.arctan2(np.sin(y[2]), np.cos(y[2]))
        
        S = H @ self.P @ H.T + self.R
        K = self.P @ H.T @ np.linalg.inv(S)
        
        self.x_hat = self.x_hat + K @ y
        self.P = (np.eye(4) - K @ H) @ self.P
        # NOTE: State wrapping removed here to ensure continuous state for NMPC
    
    def step(self, tau_applied, z_measurement):
        """Full EKF step: predict then update."""
        self.predict(tau_applied)
        self.update(z_measurement)
        return self.x_hat.copy()