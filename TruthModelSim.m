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
N = SYSTEM.N; % Number of time steps
dt = SYSTEM.dt; % [s] Time step intervals

k = 0:1:N-1; % Time step k vector
t_vec = dt*k; % [s] Time vector

% Initial Conditions
r0 = CONST.R_E + 2000; % [km]
x0_C = [r0; 0; 0; 2.57; 0; 7.05];
% x0_C = [r0; 0; 0; 7.5; 0; 0];
x0_D = [r0+100; 1; 0; 8.75; 0; 0];

x0_E = [0; 7; r0-200; 0; r0/2; 1];

% x0_C = [r0; 0; 0; 9]; % Chief satellite
% x0_D = [r0 + 50; 2.5; 0; 8.5]; % Deputy satellite

% x0_C = [r0; 0; 0; 9; 0; 0]; % Chief satellite
% x0_D = [r0 + 50; 2.5; 0; 8.5; 0; 0]; % Deputy satellite

% Noise
w_tilde(:, :, 1) = mvnrnd(zeros(SYSTEM.n/2, 1), NOISE.Q, N)'; % Process noise
w_tilde(:, :, 2) = mvnrnd(zeros(SYSTEM.n/2, 1), NOISE.Q, N)'; % Process noise
w_tilde(:, :, 3) = mvnrnd(zeros(SYSTEM.n/2, 1), NOISE.Q, N)'; % Process noise

% Numeric Integration
options = odeset('RelTol', 1e-12);

[~, X_chief] = ode45(@(t, x) TMFunc(t, x, SYSTEM, CONST, w_tilde(:, :, 1), t_vec), t_vec, x0_C, options);
[~, X_deputy] = ode45(@(t, x) TMFunc(t, x, SYSTEM, CONST, w_tilde(:, :, 2), t_vec), t_vec, x0_D, options);
[~, X_eputy] = ode45(@(t, x) TMFunc(t, x, SYSTEM, CONST, w_tilde(:, :, 3), t_vec), t_vec, x0_E, options);

X_chief = X_chief';
X_deputy = X_deputy';
X_eputy = X_eputy';

X(:, :, 1) = X_chief;
X(:, :, 2) = X_deputy;
X(:, :, 3) = X_eputy;

% Measurements
Y(:,:,1) = measurementfcn(X_deputy, X_chief, NOISE.R);  
Y(:,:,2) = measurementfcn(X_eputy, X_chief, NOISE.R);

% Determine if the Line of Sight is obstructed
X1 = [X_chief(1, :); X_chief(3, :); X_chief(5, :)]; % [km] Position vector of the chief satellite

for i = 1:2
    
    if i == 1
        X2 = [X_deputy(1, :); X_deputy(3, :); X_deputy(5, :)]; % [km] Position vector of the deputy satellite
    else
        X2 = [X_eputy(1, :); X_eputy(3, :); X_eputy(5, :)]; % [km] Position vector of the deputy satellite
    end
    
    pt = [0; 0; 0]; % Point to find the distance to
    LOS = X2 - X1; % Line of sight vector

    V_pt_LOS = cross(cross(pt - X1, LOS), LOS); % Vector w/ minimum distance from pt to LOS
    % Normalize each vector if the magnitude of V_pt_LOS is non-zero
    V_pt_LOS(:, vecnorm(V_pt_LOS) ~= 0) = V_pt_LOS(:, vecnorm(V_pt_LOS) ~= 0)./vecnorm(V_pt_LOS(:, vecnorm(V_pt_LOS) ~= 0)); % Unit LOS vector

    % Calculate the distance between pt and LOS if the magnitude of V_pt_LOS is
    % non-zero
    d_pt_LOS = vecnorm(cross(X1, LOS))./vecnorm(LOS); % [km] Distance from pt to LOS

    pt_LOS = d_pt_LOS.*V_pt_LOS; % [km] Point in 3D space on LOS closest to pt

    % Calculate the value of t in the equation:
    %   pt_LOS = X1 + t * LOS
    % to determine whether the shortest distance occurs between the satellites
    % (Corresponds to t = [0, 1]
    t = (pt_LOS - X1)./LOS;

    bool_valid_rng_prm = zeros(1, size(Y, 2));
    % bool_valid_rng_prm(1, :) = abs(t(1, :) - t(2, :)) < 1e-12;
    % bool_valid_rng_prm(2, :) = abs(t(1, :) - t(3, :)) < 1e-12;
    bool_valid_rng_prm = d_pt_LOS > CONST.R_E | t(1, :) < 0 | t(1, :) > 1;

    valid_rng_prm = double(bool_valid_rng_prm);

    valid_rng_prm(valid_rng_prm == 0) = NaN;
    Y(:,:,i) = Y(:,:,i).*valid_rng_prm;
end
end

