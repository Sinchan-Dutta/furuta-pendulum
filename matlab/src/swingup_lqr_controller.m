clc; clear; close all;

%% ── Physical parameters ──────────────────────────────────────────────────
la = 0.25;      % Arm length (m)
ma = 0.3;       % Arm mass (kg)
lp = 0.15;      % Pendulum length (m)
mp = 0.05;      % Pendulum mass (kg)
h  = 0.3;       % Base height (m)
J  = 1e-3;      % Base inertia (kg·m²)
g  = 9.81;      % Gravity (m/s²)
c  = 0.02;      % Arm damping (N·m·s/rad)
b  = 0.005;     % Pendulum damping (N·m·s/rad)

%% ── DC Motor parameters ──────────────────────────────────────────────────
R_m = 2.0;      % Armature resistance (Ω)
L_m = 5e-3;     % Armature inductance (H)
K_t = 0.05;     % Torque constant (N·m/A)
K_e = 0.05;     % Back-EMF constant (V·s/rad)  [= K_t for SI-consistent motor]
N   = 20;       % Gear ratio (motor → arm shaft)
V_max = 24;     % Supply voltage limit (V)

%  Torque delivered to arm shaft: tau = N * K_t * i
%  Back-EMF seen by motor:        e   = K_e * N * theta_dot

%% ── Inertia coefficients ─────────────────────────────────────────────────
alpha = J + (1/3)*ma*la^2 + mp*la^2;
beta  = (1/3)*mp*lp^2;
gamma_c = (1/2)*mp*la*lp;   % renamed to avoid clash with MATLAB gamma()
delta = (1/2)*mp*g*lp;

detM0 = alpha*beta - gamma_c^2;

%% ── Linearised state-space (5 states: θ, θ̇, φ, φ̇, i) ─────────────────
%
%  Mechanical equations (linearised about φ=π upright, θ=0):
%    detM0 * θ̈ =  -c*β*θ̇ + γ*δ*φ_err - b*γ*φ̇ + β*(N*K_t)*i
%    detM0 * φ̈ =  -c*γ*θ̇ + α*δ*φ_err - b*α*φ̇ + γ*(N*K_t)*i
%
%  Electrical equation:
%    L_m * i̇ = V_in - R_m*i - K_e*N*θ̇
%
%  State vector: x = [θ; θ̇; φ_err; φ̇; i]
%  Input:        u = V_in

NKt = N * K_t;   % effective torque gain seen at arm shaft
KeN = K_e * N;   % effective back-EMF coefficient

A = [0,                1,                  0,              0,              0;
     0,  -c*beta/detM0,  gamma_c*delta/detM0, -b*gamma_c/detM0, beta*NKt/detM0;
     0,                0,                  0,              1,              0;
     0, -c*gamma_c/detM0, alpha*delta/detM0, -b*alpha/detM0, gamma_c*NKt/detM0;
     0,         -KeN/L_m,                  0,              0,     -R_m/L_m];

B = [0; 0; 0; 0; 1/L_m];

%% ── Controllability check ────────────────────────────────────────────────
tol = 1e-8;
curly_C = zeros(5,5);
Bpow = B;
for i = 1:5
    curly_C(:,i) = Bpow;
    Bpow = A * Bpow;
end
if rank(curly_C, tol) == 5
    disp('System (with motor) is controllable.');
else
    warning('System may not be fully controllable — check parameters.');
end

%% ── LQR via Hamiltonian eigenvector method ───────────────────────────────
%  State weights: [θ, θ̇, φ_err, φ̇, i]
Q = diag([1, 1, 10, 1, 0.01]);
R_lqr = 0.5;   % voltage weight (higher → smoother, slower response)

H = [A, -B*(1/R_lqr)*B'; -Q, -A'];
[V_eig, D] = eig(H);
ev = diag(D);
neg_idx = find(real(ev) < 0);
M = V_eig(:, neg_idx);

X1 = M(1:5, :);
X2 = M(6:10, :);
P  = X2 / X1;
P  = 0.5*(P + P');
K  = real((1/R_lqr) * B' * P);   % 1×5 gain vector

fprintf('LQR gains K = [%.3f  %.3f  %.3f  %.3f  %.3f]\n', K);

%% ── Simulate ─────────────────────────────────────────────────────────────
% x = [theta; theta_dot; phi; phi_dot; i_motor]
x0    = [0; 0; 0.1; 0; 0];   % pendulum near bottom, zero current
tspan = [0 6];

options = odeset('RelTol',1e-6,'AbsTol',1e-7);
[t, X] = ode45(@(t,x) furuta_motor_dynamics(t, x, ...
                    alpha, beta, gamma_c, delta, c, b, ...
                    R_m, L_m, K_t, K_e, N, V_max, K), ...
               tspan, x0, options);

%% ── Plots ────────────────────────────────────────────────────────────────
figure('Color','w','Position',[100 100 900 700]);

subplot(3,2,1)
plot(t, X(:,3), 'LineWidth',1.5); yline(pi,'r--','LineWidth',1);
title('Pendulum Angle \phi'); ylabel('\phi (rad)'); grid on

subplot(3,2,2)
plot(t, X(:,1), 'LineWidth',1.5);
title('Arm Angle \theta'); ylabel('\theta (rad)'); grid on

subplot(3,2,3)
plot(t, X(:,4), 'LineWidth',1.5);
title('Pendulum Angular Velocity'); ylabel('\phi_{dot} (rad/s)'); grid on

subplot(3,2,4)
plot(t, X(:,2), 'LineWidth',1.5);
title('Arm Angular Velocity'); ylabel('\theta_{dot} (rad/s)'); grid on

subplot(3,2,5)
plot(t, X(:,5), 'LineWidth',1.5, 'Color',[0.2 0.6 0.2]);
title('Motor Armature Current'); ylabel('i (A)'); xlabel('Time (s)'); grid on

subplot(3,2,6)
% Reconstruct voltage command for plotting
V_log = zeros(length(t),1);
for k = 1:length(t)
    x_k = X(k,:)';
    phi_err = wrapToPi(x_k(3) - pi);
    if abs(phi_err) < deg2rad(40)
        x_ctrl = [x_k(1); x_k(2); phi_err; x_k(4); x_k(5)];
        V_log(k) = -K * x_ctrl;
    else
        % swing-up: estimate via back-calculation (motor eq.)
        V_log(k) = NaN;   % not LQR regime
    end
    V_log(k) = max(min(V_log(k), V_max), -V_max);
end
plot(t, V_log, 'LineWidth',1.5, 'Color',[0.8 0.2 0.2]);
title('Voltage Command V_{in}'); ylabel('V (V)'); xlabel('Time (s)'); grid on
yline(V_max,'k--'); yline(-V_max,'k--');

sgtitle('Furuta Pendulum — Swing-Up + LQR with DC Motor Model','FontWeight','bold');

%% ── Animation ────────────────────────────────────────────────────────────
figure('Color','w');
axis equal; axis([-0.4 0.4 -0.4 0.4 0 0.6]);
view(45,30); grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Furuta Pendulum: Swing-Up + LQR + DC Motor');
hold on;

th_circle = linspace(0,2*pi,100);
plot3(0.15*cos(th_circle), 0.15*sin(th_circle), zeros(1,100), 'k--');
plot3([0 0],[0 0],[0 h],'k','LineWidth',4);

arm1 = plot3([0 0],[0 0],[0 0],'b','LineWidth',3);
arm2 = plot3([0 0],[0 0],[0 0],'r','LineWidth',3);

fps = 60;
t_anim   = 0:(1/fps):max(t);
X_anim   = interp1(t, X, t_anim);

v = VideoWriter('furuta_motor.mp4','MPEG-4');
v.FrameRate = 40;
open(v);
pause(2);

for k = 1:length(t_anim)
    if ~isgraphics(arm1) || ~isgraphics(arm2)
        disp('Animation stopped by user.'); break;
    end

    theta = X_anim(k,1);
    phi   = X_anim(k,3);

    xb = la*cos(theta);
    yb = la*sin(theta);
    zb = h;

    xp = xb - lp*sin(phi)*sin(theta);
    yp = yb + lp*sin(phi)*cos(theta);
    zp = zb - lp*cos(phi);

    set(arm1,'XData',[0 xb],'YData',[0 yb],'ZData',[h zb]);
    set(arm2,'XData',[xb xp],'YData',[yb yp],'ZData',[zb zp]);
    pause(1/fps);
    drawnow;

    frame = getframe(gcf);
    writeVideo(v, frame);
end
close(v);

%% ── Dynamics function ────────────────────────────────────────────────────
function dx = furuta_motor_dynamics(~, x, ...
        alpha, beta, gamma_c, delta, c, b, ...
        R_m, L_m, K_t, K_e, N, V_max, K)

    theta     = x(1);
    theta_dot = x(2);
    phi       = x(3);
    phi_dot   = x(4);
    i_motor   = x(5);

    % Torque delivered to arm shaft by motor
    tau_motor = N * K_t * i_motor;

    phi_err = wrapToPi(phi - pi);

    %% ── Controller ───────────────────────────────────────────────────────
    if abs(phi_err) < deg2rad(40)
        % ── LQR stabilisation ──
        x_ctrl = [theta; theta_dot; phi_err; phi_dot; i_motor];
        V_in = -K * x_ctrl;
    else
        % ── Energy-based swing-up ──
        beta_  = (1/3)*(0.05)*(0.15)^2;   % reuse mp,lp inline
        delta_ = (1/2)*(0.05)*9.81*(0.15);
        E   = 0.5*beta_*phi_dot^2 + delta_*(1 - cos(phi));
        Ed  = 2*delta_;
        kE  = 30;
        eps = 0.1;

        tau_swingup = kE*(Ed - E)*tanh(phi_dot*cos(phi)/eps);
        tau_center  = -1.5*theta - 0.05*theta_dot;
        tau_desired = tau_swingup + tau_center;

        % Invert motor model to find required voltage:
        %   tau = N*K_t*i  →  i_desired = tau / (N*K_t)
        %   V = R_m*i + K_e*N*theta_dot  (quasi-static inversion for feed-fwd)
        i_desired = tau_desired / (N * K_t);
        V_in = R_m * i_desired + K_e * N * theta_dot;
    end

    % Voltage saturation
    V_in = max(min(V_in, V_max), -V_max);

    %% ── Nonlinear mechanical dynamics ────────────────────────────────────
    Den = alpha*beta + beta^2*sin(phi)^2 - (gamma_c*cos(phi))^2;

    U = -beta*sin(2*phi)*theta_dot*phi_dot ...
        + gamma_c*sin(phi)*phi_dot^2 ...
        - c*theta_dot + tau_motor;

    V_mech = 0.5*beta*sin(2*phi)*theta_dot^2 ...
             - delta*sin(phi) ...
             - b*phi_dot;

    theta_dd = (beta*U   - gamma_c*cos(phi)*V_mech) / Den;
    phi_dd   = (-gamma_c*cos(phi)*U + (alpha + beta*sin(phi)^2)*V_mech) / Den;

    %% ── Electrical dynamics (5th state) ──────────────────────────────────
    %  L_m * di/dt = V_in - R_m*i - K_e*(N*theta_dot)
    di_dt = (V_in - R_m*i_motor - K_e*N*theta_dot) / L_m;

    dx = [theta_dot;
          theta_dd;
          phi_dot;
          phi_dd;
          di_dt];
end