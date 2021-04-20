function [Y] = GetY(X_deputy, X_chief, R)
%GetY Generates measurements given the state information X with measurement
%noise R
%
% Inputs:
%   X - [km and km/s] (n x N x n_s) Matrix containing the states of each satellite
%   R - Measurement noise covariance matrix
%
% Outputs:
%   Y - [km, km/s, and rad] (p x N x (n_s - 1)) Matrix containing the chief
%   satellite's measurements of each deputy satellite

[~, N] = size(X_chief);

p = 4; % Number of measurements
Y = zeros(p, N);

[Sv, pd] = chol(R); % Cholesky factorization and not positive definite flag
vk = randn(p, N);

rho = sqrt((X_chief(1, :) - X_deputy(1, :)).^2 + (X_chief(3, :) - X_deputy(3, :)).^2 + (X_chief(5, :) - X_deputy(5, :)).^2); % [km]
rhodot = ((X_chief(1, :) - X_deputy(1, :)).*(X_chief(2, :) - X_deputy(2, :)) + ...
    (X_chief(3, :) - X_deputy(3, :)).*(X_chief(4, :) - X_deputy(4, :)) + ...
    (X_chief(5, :) - X_deputy(5, :)).*(X_chief(6, :) - X_deputy(6, :)))./rho; % [km/s]
az = atan2((X_chief(3, :) - X_deputy(3, :)), (X_chief(1, :) - X_deputy(1, :))); % [rad]
el = atan((X_chief(5, :) - X_deputy(5, :))./(X_chief(1, :) - X_deputy(1, :))); % [rad]

Y = [rho; rhodot; az; el];
if pd == 0 % If R is positive definite
    Y = Y + Sv*vk;
end

end

