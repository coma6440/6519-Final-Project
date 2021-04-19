function dx = ode_fun(t,x,CONST,w)
r = sqrt(x(1)^2 + x(3)^2); % [km]
dx = zeros(4, 1);
dx(1) = x(2);
dx(2) = -CONST.mu_E/r^3 * x(1) + w(1);
dx(3) = x(4);
dx(4) = -CONST.mu_E/r^3 * x(3) + w(2);
end