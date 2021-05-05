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
x_hat = [RBPF(:).x_mmse];
P = cat(3,RBPF(:).P);
xt_hat = [RBPF(:).xt_mmse];

figure
tiledlayout('flow')
nexttile
hold on
plot(t_vec, x_hat(1,:) - x_deputy(1,:))
plot(t_vec, 2*sqrt(squeeze(P(1,1,:))),'--k');
plot(t_vec, -2*sqrt(squeeze(P(1,1,:))),'--k');
xlabel('time [s]')
ylabel('X Error')

nexttile
hold on
plot(t_vec, x_hat(3,:) - x_deputy(3,:))
plot(t_vec, 2*sqrt(squeeze(P(3,3,:))),'--k');
plot(t_vec, -2*sqrt(squeeze(P(3,3,:))),'--k');
xlabel('time [s]')
ylabel('Y Error')

nexttile
hold on
plot(t_vec, x_hat(5,:) - x_deputy(5,:))
plot(t_vec, 2*sqrt(squeeze(P(5,5,:))),'--k');
plot(t_vec, -2*sqrt(squeeze(P(5,5,:))),'--k');
xlabel('time [s]')
ylabel('Z Error')

figure
tiledlayout('flow')
nexttile
plot(t_vec, xt_hat(1,:) - x_chief(1,:))
xlabel('time [s]')
ylabel('X Error')
nexttile
plot(t_vec, xt_hat(3,:) - x_chief(3,:))
xlabel('time [s]')
ylabel('Y Error')
nexttile
plot(t_vec, xt_hat(5,:) - x_chief(5,:))
xlabel('time [s]')
ylabel('Z Error')




