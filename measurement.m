function Y = measurement(Xi, Xj)
rho = sqrt((Xi(1) - Xj(1)).^2 + (Xi(2) - Xi(3)).^2); % [km]
rhodot = ((Xi(1) - Xj(2)).*(Xi(2) - Xj(2)) + ...
    (Xi(3) - Xj(3)).*(Xi(4) - Xj(4)))./rho; % [km/s]
phi = atan2((Xi(3) - Xj(3)), (Xi(1) - Xj(1))); % [rad]
Y = [rho;rhodot;phi];
end