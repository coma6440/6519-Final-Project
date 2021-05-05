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
x_deputy1 = X(:,:,2);
x_deputy2 = X(:,:,3);

params.p = size(Y, 1); % Number of measurements

fig = figure;
% fig.WindowState = 'maximized';
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
params.PF.Q = params.PF.Q;

params.UKF.R = 10 * NOISE.R;
params.PF.R = 10 * NOISE.R;
params.Ns = 100; % Number of samples
% x_meas = x_chief;
[RBPF] = RunRBPF(Y, params, CONST);
x_mmse1 = [RBPF(:).x_mmse1];
x_mmse2 = [RBPF(:).x_mmse2];
x_map1 = [RBPF(:).x_map1];
x_map2 = [RBPF(:).x_map2];
P1 = cat(3,RBPF(:).P1);
P2 = cat(3,RBPF(:).P2);
xt_hat = [RBPF(:).xt_mmse];
y_res_mmse1 = [RBPF(:).y_res_mmse1];
y_res_mmse2 = [RBPF(:).y_res_mmse2];
y_res_map1 = [RBPF(:).y_res_map1];
y_res_map2 = [RBPF(:).y_res_map2];
Neff = [RBPF(:).Neff];

%% State Estimation Errors
% MMSE
est_err_fig1 = PlotEstErr(t_vec, x_deputy1, x_mmse1, P1, sprintf('Particle Filter MAP, N_s = %d', params.Ns));
saveas(est_err_fig1, 'EstErr_RBPF_MMSE1.png');

est_err_fig2 = PlotEstErr(t_vec, x_deputy2, x_mmse2, P2, sprintf('Particle Filter MAP, N_s = %d', params.Ns));
saveas(est_err_fig2, 'EstErr_RBPF_MMSE2.png');

% MAP
est_err_fig1 = PlotEstErr(t_vec, x_deputy1, x_map1, P1, sprintf('Particle Filter MAP, N_s = %d', params.Ns));
saveas(est_err_fig1, 'EstErr_RBPF_MAP1.png');

est_err_fig2 = PlotEstErr(t_vec, x_deputy2, x_map2, P2, sprintf('Particle Filter MAP, N_s = %d', params.Ns));
saveas(est_err_fig2, 'EstErr_RBPF_MAP2.png');

% Chief Satellite
chief_fig = figure;
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
saveas(chief_fig, 'EstErr_RBPF_chief.png');
%% Measurement Residuals
% MMSE
resid_fig1 = PlotMeasInnov(t_vec, y_res_mmse1, sprintf('PF MMSE, N_s = %d', params.Ns));
saveas(resid_fig1, 'ResidualRBPF_MMSE1.png')

resid_fig2 = PlotMeasInnov(t_vec, y_res_mmse2, sprintf('PF MMSE, N_s = %d', params.Ns));
saveas(resid_fig2, 'ResidualRBPF_MMSE2.png')
% MAP
resid_fig1 = PlotMeasInnov(t_vec, y_res_map1, sprintf('PF MAP, N_s = %d', params.Ns));
saveas(resid_fig1, 'ResidualRBPF_MAP1.png')

resid_fig2 = PlotMeasInnov(t_vec, y_res_map2, sprintf('PF MAP, N_s = %d', params.Ns));
saveas(resid_fig2, 'ResidualRBPF_MAP2.png')
%% Effective Sample Size
Neff_fig = figure;
hold on
grid on
xlabel('Time [s]')
ylabel('Effective Sample Size, N_{eff}')
plot(t_vec, Neff, 'x', 'LineWidth', 2)
set(gca, 'FontSize', 14)
saveas(Neff_fig, 'Neff_RBPF_multi-sat.png')