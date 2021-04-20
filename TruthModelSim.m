function [t_vec, X, Y] = TruthModelSim(SYSTEM, CONST, NOISE)
%TruthModelSim Simulates the orbital dynamics of multiple satellites to
%generate the truth model used for estimation
%
% Inputs:
%   
%
% Outputs:
%   t - [s] (1 x N) Time Vector
%   X - [km and km/s] (n x N x n_s) Matrix containing the n states of n_s
%       satellites for N time steps
%   Y - [km, km/s, and rad] (p x N x (n_s - 1)) Matrix containing
%       measurements of deputy satellites taken by the chief satellite

% System Parameters
n = SYSTEM.n; % Number of states
N = SYSTEM.N; % Number of time steps
dt = SYSTEM.dt; % [s] Time step intervals

k = 0:1:N-1; % Time step k vector
t_vec = dt*k; % [s] Time vector

% Initial Conditions
r0 = CONST.R_E + 2000; % [km]
x0_C = [r0; 0; 0; 9; 0; 0];
x0_D = [r0+100; 1; 0; 8.75; 0; 0];

% x0_C = [r0; 0; 0; 9]; % Chief satellite
% x0_D = [r0 + 50; 2.5; 0; 8.5]; % Deputy satellite

% x0_C = [r0; 0; 0; 9; 0; 0]; % Chief satellite
% x0_D = [r0 + 50; 2.5; 0; 8.5; 0; 0]; % Deputy satellite

% Noise
w_tilde = mvnrnd(zeros(SYSTEM.n/2, 1), NOISE.Q, N)'; % Process noise

% Numeric Integration
options = odeset('RelTol', 1e-12);

[~, x_chief] = ode45(@(t, x) TMFunc(t, x, SYSTEM, CONST, w_tilde, t_vec), t_vec, x0_C, options);
[~, x_deputy] = ode45(@(t, x) TMFunc(t, x, SYSTEM, CONST, w_tilde, t_vec), t_vec, x0_D, options);

X(:, :, 1) = x_chief';
X(:, :, 2) = x_deputy';

% Measurements
n_s = 2; % Number of satellites
X_chief = X(:, :, 1);

for i = 2:n_s
    X_deputy = X(:, :, i);
    
    Y = GetY(X_deputy, X_chief, NOISE.R);
end

end

