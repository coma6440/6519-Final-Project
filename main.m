%% ASEN 6519 - Advanced State Estimation
%   Final Project
%   Jordan Abell and Conner Martin
%   
%   Final Project Main Script

%% Truth Model
clear; close all; clc;

% System parameters
SYSTEM.n = 4; % Number of states
SYSTEM.p = 3; % Number of measurements
SYSTEM.N = 5000; % Number of time steps
SYSTEM.dt = 50; % [s] Time step interval

% Constants
CONST.mu_E = 3.986e5; % [km^3/s^2] Earth's Gravitational Parameter
CONST.R_E = 6378; % [km] Earth's Radius

% Noise specifications
NOISE.Q = 1e-10 * eye(SYSTEM.n/2); % Process noise covariance
% NOISE.Q = zeros(SYSTEM.n/2);
NOISE.R = [ 1,      0,      0;
            0,      1,      0;
            0,      0,      0.01]; % Measurement noise covariance
% NOISE.R = zeros(3, 3);

rng(2021);

% Truth Model
[t, X, Y] = TruthModelSim(SYSTEM, CONST, NOISE);

fig = figure;
fig.WindowState = 'maximized';
PlotTrajectories(t, X, fig);
PlotMeasurements(t, Y, fig);

%% UKF
%Parameters for UKF
params.n = 2;
params.kappa = 0;
params.beta = 2;
params.alpha = 0.5;
params.lambda = params.alpha^2 * ...
                        (params.n + params.kappa) - ...
                        params.n;
params.wm = params.lambda/(params.n + params.lambda);
params.wc = params.wm + 1 - params.alpha^2 + params.beta;
params.wm(2:(2*params.n + 1)) = 1/(2*(params.n+params.lambda));
params.wc(2:(2*params.n + 1)) = 1/(2*(params.n+params.lambda));

xp = [0;0];
xm = [0;0];
Pp = 10*eye(2);
Pm = 10*eye(2);

for k = 1:150
    [xm, Pm] = predict_UKF(Pp, xp, u, params);
    [xp, Pp] = correct_UKF(Pm, xm, y, params);
end