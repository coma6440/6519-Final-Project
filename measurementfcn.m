function y = measurementfcn(x,x_chief)
rho = sqrt((x_chief(1) - x(1)).^2 + (x_chief(2) - x(2)).^2); % [km]
rhodot = ((x_chief(1) - x(1)).*(x_chief(2) - x(2)) + ...
    (x_chief(3) - x(3)).*(x_chief(4) - x(4)))./rho; % [km/s]
phi = atan2((x_chief(3) - x(3)), (x_chief(1) - x(1))); % [rad]
y = [rho;rhodot;phi];
end