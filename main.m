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

% Truth Model
[t_vec, X, Y] = TruthModelSim(params, CONST, NOISE);

u = zeros(3,1);
x_chief = X(:,:,1);
x_deputy = X(:,:,2);

params.p = size(Y, 1); % Number of measurements

fig = figure;
fig.WindowState = 'maximized';
PlotTrajectories(t_vec, X, fig);
PlotMeasurements(t_vec, Y, fig);
saveas(fig, 'TM.png')

%% Execute Filters
%Parameters
params.Q = 1e-6*[   0,      0,      0,      0,      0,      0;
                    0,      1,      0,      0,      0,      0;
                    0,      0,      0,      0,      0,      0;
                    0,      0,      0,      1,      0,      0;
                    0,      0,      0,      0,      0,      0;
                    0,      0,      0,      0,      0,      1];

params.R = 10 * NOISE.R;
% params.R = 10*[1,      0,      0,      0;
%             0,      1e-5,      0,      0;
%             0,      0,      0.01,    0; 
%             0,      0,      0,      0.01]; % Measurement noise covariance

[UKF, PF] = RunFilters(x_chief, Y, params, CONST);


%% State Estimation Errors for UKF

est_err_fig = PlotEstErr(t_vec, x_deputy, UKF.xp, UKF.Pp, 'UKF');
saveas(est_err_fig, 'EstErr_UKF.png')


%% Measurement Errors for UKF

resid_fig = PlotMeasInnov(t_vec, UKF.yres, 'UKF');
saveas(resid_fig, 'ResidualUKF.png')


%% State Estimation Errors for PF

PF_vec = VectorizeResults(PF);
est_err_fig = PlotEstErr(t_vec, x_deputy, PF_vec.x_mmse, PF_vec.P, 'Particle Filter, N_s = 1000');
saveas(est_err_fig, 'EstErr_PF_MMSE.png');

PF_vec = VectorizeResults(PF);
est_err_fig = PlotEstErr(t_vec, x_deputy, PF_vec.x_map, PF_vec.P, 'Particle Filter, N_s = 1000');
saveas(est_err_fig, 'EstErr_PF_MAP.png');

%% Measurement Errors for PF
%MMSE
resid_fig = PlotMeasInnov(t_vec, PF_vec.y_res_mmse, 'PF MMSE');
saveas(resid_fig, 'ResidualPF_MMSE.png')
%MAP
resid_fig = PlotMeasInnov(t_vec, PF_vec.y_res_map, 'PF MAP');
saveas(resid_fig, 'ResidualPF_MAP.png')