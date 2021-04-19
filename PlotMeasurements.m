function [fig] = PlotMeasurements(t, Y, fig)
%PlotMeasurements Plots the time history of measurements of each satellite
%taken by the chief satellite

[p, N, n_sm1] = size(Y);

subplot_i = p; % Subplot row index

if nargin < 3
    fig = figure;
    fig.Position = [768.2,41.8,766.4,740.8];
    subplot_j = 1; % Subplot column index
else
    figure(fig);
    subplot_j = 2;
end

ttlStr = ["Range"; "Range Rate"; "Azimuth"; "Elevation"];
ylblStr = ["$\rho$ [km]"; "$\dot{\rho}$ [km/s]"; "$\theta$ [rad]"; "$\phi$ [rad]"];

for i = 1:p
    subplot(subplot_i, subplot_j, subplot_j*i);
    hold on
    grid on
    title(ttlStr(i));
    ylabel(ylblStr(i), 'Interpreter', 'latex');
    if i == p
        xlabel('Time [s]')
    end

    for j = 1:n_sm1
        subplt_sat = sprintf('Deputy Satellite %d', j);
        
        plot(t, Y(i, :, j), 'DisplayName', subplt_sat);
        
        if i == 1
            legend('location', 'best')
        end
    end
    set(gca, 'FontSize', 14);
end

end