%% ASEN 6519 - Advanced State Estimation
%   Final Project
%   Jordan Abell and Conner Martin
%   
%   Final Project Main Script

%% Truth Model
clear; close all; clc;

% System parameters
SYSTEM.n = 6; % Number of states
SYSTEM.N = 940; % Number of time steps
SYSTEM.dt = 50; % [s] Time step interval

% Constants
CONST.mu_E = 3.986e5; % [km^3/s^2] Earth's Gravitational Parameter
CONST.R_E = 6378; % [km] Earth's Radius

% Noise specifications
NOISE.Q = 1e-10 * eye(SYSTEM.n/2); % Process noise covariance
% NOISE.Q = zeros(SYSTEM.n/2);
NOISE.R = [ 0.1,    0,      0,      0;
            0,      0.1,   0,      0;
            0,      0,      0.01,   0;
            0,      0,      0,      0.01]; % Measurement noise covariance

% NOISE.R = zeros(3, 3);

rng(2021);

% Truth Model
[t, X, Y] = TruthModelSim(SYSTEM, CONST, NOISE);

SYSTEM.p = length(Y); % Number of measurements

fig = figure;
fig.WindowState = 'maximized';
PlotTrajectories(t, X, fig);
PlotMeasurements(t, Y, fig);
saveas(fig, 'TM.png')

%% UKF
%Parameters for UKF
params.n = 6;
params.p = 4;
params.Q = 2.5e-7*eye(6);
% params.Q(1,1) = 0;
% params.Q(3,3) = 0;
% params.Q(5,5) = 0;

% params.R = NOISE.R;
params.R = [1,      0,      0,      0;
            0,      1,      0,      0;
            0,      0,      0.1,    0;
            0,      0,      0,      0.1]; % Measurement noise covariance
params.dt = SYSTEM.dt;

r0 = CONST.R_E + 2000; % [km]
init_state = [r0+100; 1; 0; 8.75; 0; 0];
xmUKF = zeros(params.n,SYSTEM.N);
xpUKF = xmUKF;
xmUKF(:,1) = init_state;
xpUKF(:,1) = init_state;

xmPF = zeros(params.n,SYSTEM.N);
xpPF = xmUKF;
xmPF(:,1) = init_state;
xpPF(:,1) = init_state;

PpUKF = 1*eye(6);
PpUKF(2,2) = 0.0005;
PpUKF(4,4) = 0.0005;
PpUKF(6,6) = 0.0005;
PpPF = PpUKF;

PmUKF = PpUKF;
u = zeros(3,1);
x_chief = X(:,:,1);
x_deputy = X(:,:,2);
t_vec = t;
t = 0;
%UKF
filterUKF = trackingUKF(@transitionfcn, @GetY, init_state, ...
         'ProcessNoise', params.Q, 'MeasurementNoise', NOISE.R,...
         'StateCovariance', PpUKF,...
         'Alpha', 1e-3);
%Particle Filter  
filterPF = trackingPF(@particletransitionfcn, @particlemeasurementfcn, ...
           init_state, 'ProcessNoise', params.Q, 'MeasurementNoise', NOISE.R,...
           'StateCovariance', PpPF, 'NumParticles', 500);

for k = 1:(length(Y)-1)
    %Prediction step
    [xmUKF(:,k+1), PmUKF(:,:,k+1)] = predict(filterUKF,params.dt,CONST);
    [xmPF(:,k+1), PmPF(:,:,k+1)] = predict(filterPF,params.dt,CONST);
    %Measurement residual
    [yresUKF(:,k+1), yrescovUKF(:,:,k+1)] = residual(filterUKF, Y(:,k+1), {x_chief(:,k+1), zeros(params.p, params.p)});
    %[yresPF(:,k+1), yrescovPF(:,:,k+1)] = residual(filterPF, Y(:,k+1), {x_chief(:,k+1), zeros(params.p, params.p)});
    %Correction step
    [xpUKF(:,k+1), PpUKF(:,:,k+1)] = correct(filterUKF,Y(:,k+1),x_chief(:,k+1), zeros(params.p, params.p));
    [xpPF(:,k+1), PpPF(:,:,k+1)] = correct(filterPF,Y(:,k+1),x_chief(:,k+1), zeros(params.p, params.p));
    %Step time forward
    t = t + params.dt;
end


%% State Estimation Errors for UKF
version_flag = contains(version, '2020');

ex = (xpUKF - x_deputy);
est_err_fig = figure;
est_err_fig.Position = [268,98,885,601];
hold on
sgtitle('State Estimation Errors')

if version_flag
    tiledlayout('flow','TileSpacing','Compact')
    nexttile
else
    subplot(3, 2, 1)
end
hold on
plot(t_vec,ex(1,:))
plot(t_vec,2*sqrt(squeeze(PpUKF(1,1,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(PpUKF(1,1,:))),'--k')
ylabel({'X Position', 'Error [km]'})
grid on
set(gca, 'FontSize', 14)

if version_flag
    nexttile
else
    subplot(3, 2, 2)
end
hold on
plot(t_vec,ex(2,:))
plot(t_vec,2*sqrt(squeeze(PpUKF(2,2,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(PpUKF(2,2,:))),'--k')
ylabel({'X Velocity', 'Error [km/s]'})
grid on
set(gca, 'FontSize', 14)

if version_flag
    nexttile
else
    subplot(3, 2, 3)
end
hold on
plot(t_vec,ex(3,:))
plot(t_vec,2*sqrt(squeeze(PpUKF(3,3,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(PpUKF(3,3,:))),'--k')
ylabel({'Y Position', 'Error [km]'})
grid on
set(gca, 'FontSize', 14)

if version_flag
    nexttile
else
    subplot(3, 2, 4)
end
hold on
plot(t_vec,ex(4,:))
plot(t_vec,2*sqrt(squeeze(PpUKF(4,4,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(PpUKF(4,4,:))),'--k')
ylabel({'Y Velocity', 'Error [km/s]'});
grid on
set(gca, 'FontSize', 14)

if version_flag
    nexttile
else
    subplot(3, 2, 5)
end
hold on
plot(t_vec,ex(5,:))
plot(t_vec,2*sqrt(squeeze(PpUKF(5,5,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(PpUKF(5,5,:))),'--k')
xlabel('Time [s]')
ylabel({'Z Position', 'Error [km]'})
grid on
set(gca, 'FontSize', 14)

if version_flag
    nexttile
else
    subplot(3, 2, 6)
end
hold on
plot(t_vec,ex(6,:))
plot(t_vec,2*sqrt(squeeze(PpUKF(6,6,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(PpUKF(6,6,:))),'--k')
xlabel('Time [s]')
ylabel({'Z Velocity', 'Error [km/s]'});
grid on
set(gca, 'FontSize', 14)

saveas(est_err_fig, 'EstErr_UKF.png')

%% Measurement Errors for UKF
resid_fig = figure;
hold on
sgtitle('Measurement Residuals')
resid_fig.Position = [268,98,885,601];
if version_flag
    tiledlayout('flow','TileSpacing','Compact')
    nexttile
else
    subplot(2, 2, 1)
end
hold on
plot(t_vec,yresUKF(1,:))
ylabel(["Range Residual", "[km]"])
grid on
set(gca, 'FontSize', 14)

if version_flag
    nexttile
else
    subplot(2, 2, 2)
end
hold on
plot(t_vec,yresUKF(2,:))
ylabel(["Range Rate Residual", "[km/s]"])
grid on
set(gca, 'FontSize', 14)

if version_flag
    nexttile
else
    subplot(2, 2, 3)
end
hold on
plot(t_vec,yresUKF(3,:))
xlabel('Time [s]')
ylabel(["Azimuth Error", "[rad]"])
grid on
set(gca, 'FontSize', 14)

if version_flag
    nexttile
else
    subplot(2, 2, 4)
end
hold on
plot(t_vec,yresUKF(4,:))
xlabel('Time [s]')
ylabel(["Elevation Error", "[rad]"])
grid on
set(gca, 'FontSize', 14)

saveas(resid_fig, 'ResidualUKF.png')


%% State Estimation Errors for PF
version_flag = contains(version, '2020');

ex = (xpPF - x_deputy);
est_err_fig = figure;
est_err_fig.Position = [268,98,885,601];
hold on
sgtitle('State Estimation Errors')

if version_flag
    tiledlayout('flow','TileSpacing','Compact')
    nexttile
else
    subplot(3, 2, 1)
end
hold on
plot(t_vec,ex(1,:))
plot(t_vec,2*sqrt(squeeze(PpPF(1,1,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(PpPF(1,1,:))),'--k')
ylabel({'X Position', 'Error [km]'})
grid on
set(gca, 'FontSize', 14)

if version_flag
    nexttile
else
    subplot(3, 2, 2)
end
hold on
plot(t_vec,ex(2,:))
plot(t_vec,2*sqrt(squeeze(PpPF(2,2,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(PpPF(2,2,:))),'--k')
ylabel({'X Velocity', 'Error [km/s]'})
grid on
set(gca, 'FontSize', 14)

if version_flag
    nexttile
else
    subplot(3, 2, 3)
end
hold on
plot(t_vec,ex(3,:))
plot(t_vec,2*sqrt(squeeze(PpPF(3,3,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(PpPF(3,3,:))),'--k')
ylabel({'Y Position', 'Error [km]'})
grid on
set(gca, 'FontSize', 14)

if version_flag
    nexttile
else
    subplot(3, 2, 4)
end
hold on
plot(t_vec,ex(4,:))
plot(t_vec,2*sqrt(squeeze(PpPF(4,4,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(PpPF(4,4,:))),'--k')
ylabel({'Y Velocity', 'Error [km/s]'});
grid on
set(gca, 'FontSize', 14)

if version_flag
    nexttile
else
    subplot(3, 2, 5)
end
hold on
plot(t_vec,ex(5,:))
plot(t_vec,2*sqrt(squeeze(PpPF(5,5,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(PpPF(5,5,:))),'--k')
xlabel('Time [s]')
ylabel({'Z Position', 'Error [km]'})
grid on
set(gca, 'FontSize', 14)

if version_flag
    nexttile
else
    subplot(3, 2, 6)
end
hold on
plot(t_vec,ex(6,:))
plot(t_vec,2*sqrt(squeeze(PpPF(6,6,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(PpPF(6,6,:))),'--k')
xlabel('Time [s]')
ylabel({'Z Velocity', 'Error [km/s]'});
grid on
set(gca, 'FontSize', 14)

saveas(est_err_fig, 'EstErr_PF.png')
