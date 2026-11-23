import numpy as np

class SensorModel:
    """Adds realistic noise to the true state."""
    
    def __init__(self, sigma_theta=0.002, sigma_phi=0.003,
                 sigma_theta_dot=0.05, sigma_phi_dot=0.08):
        # sigmas in rad, rad, rad/s, rad/s
        self.sigma = np.array([sigma_theta, sigma_theta_dot,
                               sigma_phi, sigma_phi_dot])
    
    def measure(self, x_true):
        """Return noisy measurement z = x_true + noise."""
        noise = np.random.normal(0, self.sigma)
        z = np.zeros(4)
        
        # Apply noise first, then wrap the measured angles
        z[0] = np.arctan2(np.sin(x_true[0] + noise[0]), np.cos(x_true[0] + noise[0]))
        z[1] = x_true[1] + noise[1]
        z[2] = np.arctan2(np.sin(x_true[2] + noise[2]), np.cos(x_true[2] + noise[2]))
        z[3] = x_true[3] + noise[3]

        return z