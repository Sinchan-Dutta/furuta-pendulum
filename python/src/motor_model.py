def voltage_to_torque(V, theta_dot, motor_params):
  Kt = motor_params['Kt']
  Kb = motor_params['Kb']
  R = motor_params['R']
  N = motor_params['N']
  
  tau = N * Kt * (V - Kb * N * theta_dot) / R

  return tau


if __name__ == "__main__":
  motor_params = {
    "Kt" : 0.042,
    "Kb" : 0.042,
    "R" : 8.4,
    "N" : 8
  }

  print(voltage_to_torque(12, 0, motor_params=motor_params))