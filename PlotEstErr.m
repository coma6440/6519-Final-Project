function [est_err_fig] = PlotEstErr(t_vec, x_tm, xp, Pp, ttlStr)
%PlotEstErr Plots the estimation errors between the truth model states
%(x_tm) and the estimated states (xp) with the filter covariance matrices
%Pp

version_flag = contains(version, '2020');

ex = (xp - x_tm);
est_err_fig = figure;
est_err_fig.Position = [268,98,885,601];
hold on
sgtitle(sprintf('State Estimation Errors - %s', ttlStr))

if version_flag
    tiledlayout('flow','TileSpacing','Compact')
    nexttile
else
    subplot(3, 2, 1)
end
hold on
plot(t_vec,ex(1,:))
plot(t_vec,2*sqrt(squeeze(Pp(1,1,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(Pp(1,1,:))),'--k')
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
plot(t_vec,2*sqrt(squeeze(Pp(2,2,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(Pp(2,2,:))),'--k')
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
plot(t_vec,2*sqrt(squeeze(Pp(3,3,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(Pp(3,3,:))),'--k')
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
plot(t_vec,2*sqrt(squeeze(Pp(4,4,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(Pp(4,4,:))),'--k')
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
plot(t_vec,2*sqrt(squeeze(Pp(5,5,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(Pp(5,5,:))),'--k')
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
plot(t_vec,2*sqrt(squeeze(Pp(6,6,:))),'--k')
plot(t_vec,-2*sqrt(squeeze(Pp(6,6,:))),'--k')
xlabel('Time [s]')
ylabel({'Z Velocity', 'Error [km/s]'});
grid on
set(gca, 'FontSize', 14)

end

