function [dXdt] = TMFunc(t, X, SYSTEM, CONST, w_tilde, t_vec)
%TMFunc Function used for numeric integration of the states to generate the
%Truth Model

x = X(1); % [km]
xdot = X(2); % [km/s]
y = X(3); % [km]
ydot = X(4); % [km/s]
z = X(5); % [km]
zdot = X(6); % [km/s]

% Process noise at time t
dt = SYSTEM.dt;
t_ind = find(t_vec <= t + dt & t_vec > t - dt);
w_til = w_tilde(:, t_ind(1));

% r = sqrt(x^2 + y^2); % [km]
r = sqrt(x^2 + y^2 + z^2); % [km]

dXdt = zeros(SYSTEM.n, 1);
dXdt(1) = xdot;
dXdt(2) = -CONST.mu_E/r^3 * x + w_til(1);
dXdt(3) = ydot;
dXdt(4) = -CONST.mu_E/r^3 * y + w_til(2);
dXdt(5) = zdot;
dXdt(6) = -CONST.mu_E/r^3 * z + w_til(3);

end

