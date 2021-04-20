function xp = transitionfcn(x, dt, CONST)
options = odeset('RelTol', 1e-12);
[~,out] = ode45(@(t, x) ode_fun(t, x, CONST, [0,0,0]), [0, dt], x, options);
xp = out(end,:)';
end
