function [resid_fig] = PlotMeasInnov(t_vec, yres, ttlStr)
%PlotMeasInnov Plots the measurement innovations

version_flag = contains(version, '2020');

resid_fig = figure;
hold on
sgtitle(sprintf('Measurement Residuals - %s', ttlStr))
resid_fig.Position = [268,98,885,601];
if version_flag
    tiledlayout('flow','TileSpacing','Compact')
    nexttile
else
    subplot(2, 2, 1)
end
hold on
plot(t_vec,yres(1,:))
ylabel(["Range Residual", "[km]"])
grid on
set(gca, 'FontSize', 14)

if version_flag
    nexttile
else
    subplot(2, 2, 2)
end
hold on
plot(t_vec,yres(2,:))
ylabel(["Range Rate Residual", "[km/s]"])
grid on
set(gca, 'FontSize', 14)

if version_flag
    nexttile
else
    subplot(2, 2, 3)
end
hold on
plot(t_vec,yres(3,:))
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
plot(t_vec,yres(4,:))
xlabel('Time [s]')
ylabel(["Elevation Error", "[rad]"])
grid on
set(gca, 'FontSize', 14)

end

