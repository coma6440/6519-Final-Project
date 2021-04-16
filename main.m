%% UKF
%Parameters for UKF
params.n = 2;
params.kappa = 0;
params.beta = 2;
params.alpha = 0.5;
params.lambda = params.alpha^2 * ...
                        (params.n + params.kappa) - ...
                        params.n;
params.wm = params.lambda/(params.n + params.lambda);
params.wc = params.wm + 1 - params.alpha^2 + params.beta;
params.wm(2:(2*params.n + 1)) = 1/(2*(params.n+params.lambda));
params.wc(2:(2*params.n + 1)) = 1/(2*(params.n+params.lambda));

xp = [0;0];
xm = [0;0];
Pp = 10*eye(2);
Pm = 10*eye(2);

for k = 1:150
    [xm, Pm] = predict_UKF(Pp, xp, u, params);
    [xp, Pp] = correct_UKF(Pm, xm, y, params);
end