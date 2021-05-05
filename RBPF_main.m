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
% x0_E = [r0 - 200; 1; 0; 8.75; 0; 0];

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
                    
rough = 0.1e-6*[    0,      0,      0,      0,      0,      0;
                  0,      1,      0,      0,      0,      0;
                  0,      0,      0,      0,      0,      0;
                  0,      0,      0,      1,      0,      0;
                  0,      0,      0,      0,      0,      0;
                  0,      0,      0,      0,      0,      1];
params.PF.Q = params.PF.Q + rough;

params.UKF.R = 10 * NOISE.R;
params.PF.R = 10 * NOISE.R;
params.Ns = 100; % Number of samples
% x_meas = x_chief;
[RBPF] = RunRBPF(Y, params, CONST);
x_mmse = [RBPF(:).x_mmse];
x_map = [RBPF(:).x_map];
P = cat(3,RBPF(:).P);
xt_hat = [RBPF(:).xt_mmse];
y_res_mmse = [RBPF(:).y_res_mmse];
y_res_map = [RBPF(:).y_res_map];
Neff = [RBPF(:).Neff];

%% State Estimation Errors
% MMSE
est_err_fig = PlotEstErr(t_vec, x_deputy, x_mmse, P, sprintf('Particle Filter MAP, N_s = %d', params.Ns));
saveas(est_err_fig, 'EstErr_RBPF_MMSE.png');

% MAP
est_err_fig = PlotEstErr(t_vec, x_deputy, x_map, P, sprintf('Particle Filter MAP, N_s = %d', params.Ns));
saveas(est_err_fig, 'EstErr_RBPF_MAP.png');

% Chief Satellite
figure
tiledlayout('flow')
nexttile
plot(t_vec, xt_hat(1,:) - x_chief(1,:))
xlabel('Time [s]')
ylabel('X Position Error [km]')

nexttile
plot(t_vec, xt_hat(2,:) - x_chief(2,:))
xlabel('Time [s]')
ylabel('X Velocity Error [km/s]')

nexttile
plot(t_vec, xt_hat(3,:) - x_chief(3,:))
xlabel('Time [s]')
ylabel('Y Position Error [km]')

nexttile
plot(t_vec, xt_hat(4,:) - x_chief(4,:))
xlabel('Time [s]')
ylabel('Y Velocity Error [km/s]')

nexttile
plot(t_vec, xt_hat(5,:) - x_chief(5,:))
xlabel('Time [s]')
ylabel('Z Position Error [km]')

nexttile
plot(t_vec, xt_hat(6,:) - x_chief(6,:))
xlabel('Time [s]')
ylabel('Z Velocity Error [km/s]')

%% Measurement Residuals
% MMSE
resid_fig = PlotMeasInnov(t_vec, y_res_mmse, sprintf('PF MMSE, N_s = %d', params.Ns));
saveas(resid_fig, 'ResidualRBPF_MMSE.png')
% MAP
resid_fig = PlotMeasInnov(t_vec, y_res_map, sprintf('PF MAP, N_s = %d', params.Ns));
saveas(resid_fig, 'ResidualRBPF_MAP.png')

%% Effective Sample Size
Neff_fig = figure;
hold on
grid on
xlabel('Time [s]')
ylabel('Effective Sample Size, N_{eff}')
plot(t_vec, Neff, 'x', 'LineWidth', 2)
set(gca, 'FontSize', 14)
saveas(Neff_fig, 'Neff_RBPF.png')