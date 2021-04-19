function [Y] = GetY(X, R)
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

[~, N, n_s] = size(X);

p = 4; % Number of measurements
Y = zeros(p, N, n_s-1);

[Sv, pd] = chol(R); % Cholesky factorization and not positive definite flag
vk = randn(p, N);

for i = 2:n_s
    rho = sqrt((X(1, :, 1) - X(1, :, i)).^2 + (X(3, :, 1) - X(3, :, i)).^2); % [km]
    rhodot = ((X(1, :, 1) - X(1, :, i)).*(X(2, :, i) - X(2, :, i)) + ...
        (X(3, :, 1) - X(3, :, i)).*(X(4, :, 1) - X(4, :, i)))./rho; % [km/s]
    az = atan2((X(3, :, 1) - X(3, :, i)), (X(1, :, 1) - X(1, :, i))); % [rad]
    el = atan2((X(5, :, 1) - X(5, :, i)), (X(1, :, 1) - X(1, :, i))); % [rad]
    
    Y = [rho; rhodot; az; el];
    
    if pd == 0 % If R is positive definite
        Y = Y + Sv*vk;
    end
end

end

