function [fig] = PlotTrajectories(t, X, fig)
%PlotTrajectories Plots the orbital trajectories of each satellite from the
%given states X
%
% Inputs:
%   t - [s] (1 x N) Time vector
%   X - [km and km/s] (n x N x n_s) Matrix containing the n states of n_s satellites for
%       each time step from 0 to N
%   fig - Figure handle
%
% Outputs:
%   fig - Figure handle

[n, N, n_s] = size(X);

if n == 4
    X(5:6, :, :) = 0;
end

if nargin < 3
    fig = figure;
    % fig.WindowState = 'maximized';
    fig.Position = [1.8,41.8,766.4,740.8];
else
    figure(fig);
    subplot(3, 2, [1, 3, 5]);
end

earth_sphere;
hold on
grid on
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
title('Orbital Trajectories of Each Satellite')
quiver3(0, 0, 0, 10000, 0, 0, 'k', 'LineWidth', 3)
quiver3(0, 0, 0, 0, 10000, 0, 'k', 'LineWidth', 3)
quiver3(0, 0, 0, 0, 0, 10000, 'k', 'LineWidth', 3)
text(10200, 100, 100, 'X', 'FontSize', 12, 'FontWeight', 'bold')
text(100, 10200, 100, 'Y', 'FontSize', 12, 'FontWeight', 'bold')
text(100, 100, 10200, 'Z', 'FontSize', 12, 'FontWeight', 'bold')
set(gca, 'FontSize', 14)

p = [];

for i = n_s:-1:1 % For each satellite
    x = X(1, :, i); % [km] X-position
    y = X(3, :, i); % [km] Y-position
    z = X(5, :, i); % [km] Z-position
    
    if i == 1
        lgndStr = sprintf('Chief Satellite');
    else
        lgndStr = sprintf('Deputy Satellite %d', i-1);
    end
    
    p(end+1) = plot3(x, y, z, 'LineWidth', 2, 'DisplayName', lgndStr);
end

legend([p], 'location', 'northeast')
% plot3(x1_1, x3_1, x5_1, 'LineWidth', 2)
% plot3(x1_2, x3_2, x5_2, 'LineWidth', 2)
end

