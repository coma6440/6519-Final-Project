function [X_kp1] = SIS(X_k, y_kp1, x_chief, params, CONST)
%SIS Sequential Importance Sampling for the discrete time system

[~, Ns] = size(X_k.x);

X_kp1.x = zeros(params.n, Ns);
X_kp1.w = zeros(Ns, 1);

w_total = 0;
iter = 0;
max_iter = 100;

while w_total == 0 && iter < max_iter
    iter = iter + 1;
    
    for i = 1:Ns
        mu_q = transitionfcn(X_k.x(:, i), params.dt, CONST);
        X_kp1.x(:, i) = mvnrnd(mu_q, params.PF.Q)';

        if ~isnan(y_kp1) % If a valid measurement is taken at time k+1
            mu_p = measurementfcn(X_kp1.x(:, i), x_chief, zeros(params.p));
            X_kp1.w(i) = X_k.w(i)*mvnpdf(y_kp1, mu_p, params.PF.R);
        else
            X_kp1.w(i) = X_k.w(i);
        end
    end

    w_total = sum(X_kp1.w);

    if w_total ~= 0
        X_kp1.w = X_kp1.w./w_total; % Normalize
    else
        warning('Sum of the particle weights equals 0. All particles have died.');
        warning('Resampling on iteration %d.', iter);
    end
end

end

