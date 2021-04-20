function y = measurementfcn(x,x_chief)
    rho = sqrt((x_chief(1) - x(1)).^2 + (x_chief(3) - x(3)).^2 + (x_chief(5) - x(5)).^2); % [km]
    rhodot = ((x_chief(1) - x(1)).*(x_chief(2) - x(2)) + ...
        (x_chief(3) - x(3)).*(x_chief(4) - x(4)) + ...
        (x_chief(5) - x(5)).*(x_chief(6)- x(6)))./rho; % [km/s]
    az = atan2((x_chief(3) - x(3)), (x_chief(1) - x(1))); % [rad]
    el = atan((x_chief(5) - x(5))/(x_chief(1) - x(1))); % [rad]
    
    y = [rho; rhodot; az; el];
end