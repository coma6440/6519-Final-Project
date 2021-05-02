%% ASEN 6519 - Advanced State Estimation
%   Final Project
%   Jordan Abell and Conner Martin
%   
%   Final Project Main Script

%% Truth Model
clear; close all; clc;

% System parameters
params.n = 6; % Number of states
params.N = 940; % Number of time steps
params.dt = 50; % [s] Time step interval

% Constants
CONST.mu_E = 3.986e5; % [km^3/s^2] Earth's Gravitational Parameter
CONST.R_E = 6378; % [km] Earth's Radius

% Noise specifications
NOISE.Q = 1e-10 * eye(params.n/2); % Process noise covariance
% NOISE.Q = zeros(params.n/2);
NOISE.R = [ 0.5,    0,      0,      0;
            0,      5e-6,   0,      0;
            0,      0,      0.001,   0;
            0,      0,      0,      0.001]; % Measurement noise covariance

% NOISE.R = zeros(4);

rng(04122021);

% Initial Conditions
r0 = CONST.R_E + 2000; % [km]
x0_C = [r0; 0; 0; 2.57; 0; 7.05];
% x0_C = [r0; 0; 0; 7.5; 0; 0];
x0_D = [r0+100; 1; 0; 8.75; 0; 0];
x0_E = [r0 - 200; 1; 0; 8.75; 0; 0];

% Truth Model
[t_vec, X, Y] = TruthModelSim(params, CONST, NOISE);

u = zeros(3,1);
x_chief = X(:,:,1);
x_deputy = X(:,:,2);
x_eputy = X(:,:,3);

params.p = size(Y, 1); % Number of measurements

fig = figure;
fig.WindowState = 'maximized';
PlotTrajectories(t_vec, X, fig);
PlotMeasurements(t_vec, Y, fig);
saveas(fig, 'TM.png')

%% Execute Filters
%Parameters
params.UKF.Q = 1e-7*[   0,      0,      0,      0,      0,      0;
                        0,      1,      0,      0,      0,      0;
                        0,      0,      0,      0,      0,      0;
                        0,      0,      0,      1,      0,      0;
                        0,      0,      0,      0,      0,      0;
                        0,      0,      0,      0,      0,      1];
                    
params.PF.Q = 1e-6*[    0,      0,      0,      0,      0,      0;
                        0,      1,      0,      0,      0,      0;
                        0,      0,      0,      0,      0,      0;
                        0,      0,      0,      1,      0,      0;
                        0,      0,      0,      0,      0,      0;
                        0,      0,      0,      0,      0,      1];

params.UKF.R = 10 * NOISE.R;
params.PF.R = 10 * NOISE.R;

params.Ns = 500; % Number of samples
[UKF, PF] = RunRBPF(x_chief, Y, params, CONST);


