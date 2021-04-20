%% ASEN 6519 - Advanced State Estimation
%   Final Project
%   Jordan Abell and Conner Martin
%   
%   Final Project Main Script

%% Truth Model
clear; close all; clc;

% System parameters
SYSTEM.n = 6; % Number of states
SYSTEM.N = 5000; % Number of time steps
SYSTEM.dt = 50; % [s] Time step interval

% Constants
CONST.mu_E = 3.986e5; % [km^3/s^2] Earth's Gravitational Parameter
CONST.R_E = 6378; % [km] Earth's Radius

% Noise specifications
NOISE.Q = 1e-10 * eye(SYSTEM.n/2); % Process noise covariance
% NOISE.Q = zeros(SYSTEM.n/2);
NOISE.R = [ 1,      0,      0,      0;
            0,      1,      0,      0;
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

%% UKF
%Parameters for UKF
params.n = 6;
params.p = 3;
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
params.Q = 1e-10*eye(6);
params.Q(1,1) = 0;
params.Q(3,3) = 0;
params.Q(5,5) = 0;

params.R = NOISE.R;
params.dt = SYSTEM.dt;

r0 = CONST.R_E + 2000; % [km]
init_state = [r0+100; 1; 0; 8.75; 0; 0];
xm = zeros(params.n,SYSTEM.N);
xp = xm;
xm(:,1) = init_state;
xp(:,1) = init_state;
Pp = 10*eye(6);
Pm = Pp;
u = zeros(3,1);
x_chief = X(:,:,1);
x_deputy = X(:,:,2);
t_vec = t;
t = 0;
filter = trackingUKF(@transitionfcn, @measurementfcn, init_state, ...
         'ProcessNoise', params.Q, 'MeasurementNoise', NOISE.R,...
         'StateCovariance', Pp,...
         'Alpha', 0.05);
for k = 1:(length(Y)-1)
    [xm(:,k+1), Pm(:,:,k+1)] = predict(filter,params.dt,CONST);
    ym(:,k+1) = measurementfcn(xm(:,k+1),x_chief(:,k+1));
    [xp(:,k+1), Pp(:,:,k+1)] = correct(filter,Y(:,k+1),x_chief(:,k+1));
    t = t + params.dt;
end

%% State Estimation Errors for UKF
ex = (xp - x_deputy);
figure
tiledlayout('flow','TileSpacing','Compact')
nexttile
hold on
plot(t_vec,ex(1,:))
plot(t_vec,2*sqrt(squeeze(Pp(1,1,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(Pp(1,1,:))),'--k')
ylabel('X Position Error')

nexttile
hold on
plot(t_vec,ex(2,:))
plot(t_vec,2*sqrt(squeeze(Pp(2,2,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(Pp(2,2,:))),'--k')
ylabel('X Velocity Error')

nexttile
hold on
plot(t_vec,ex(3,:))
plot(t_vec,2*sqrt(squeeze(Pp(3,3,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(Pp(3,3,:))),'--k')
ylabel('Y Position Error')

nexttile
hold on
plot(t_vec,ex(4,:))
plot(t_vec,2*sqrt(squeeze(Pp(4,4,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(Pp(4,4,:))),'--k')
ylabel('Y Velocity Error');

nexttile
hold on
plot(t_vec,ex(5,:))
plot(t_vec,2*sqrt(squeeze(Pp(5,5,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(Pp(5,5,:))),'--k')
ylabel('Z Position Error')

nexttile
hold on
plot(t_vec,ex(6,:))
plot(t_vec,2*sqrt(squeeze(Pp(6,6,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(Pp(6,6,:))),'--k')
ylabel('Z Velocity Error');

%% Measurement Errors for UKF
figure
tiledlayout('flow','TileSpacing','Compact')
nexttile
hold on
plot(t_vec,ym(1,:) - Y(1,:))
ylabel("Range Residual")

nexttile
hold on
plot(t_vec,ym(2,:) - Y(2,:))
ylabel("Range Rate Residual")

nexttile
hold on
plot(t_vec,ym(3,:) - Y(3,:))
ylabel("Azimuth Error")

nexttile
hold on
plot(t_vec,ym(4,:) - Y(4,:))
ylabel("Elevation Error")