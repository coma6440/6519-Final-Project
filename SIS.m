function [X_kp1] = SIS(X_k, y_kp1, x_chief, params, CONST)
%SIS Sequential Importance Sampling for the discrete time system

[~, Ns] = size(X_k.x);

X_kp1.x = zeros(params.n, Ns);
X_kp1.w = zeros(Ns, 1);

for i = 1:Ns
    mu_q = transitionfcn(X_k.x(:, i), params.dt, CONST);
    X_kp1.x(:, i) = mvnrnd(mu_q, params.Q)';
    
    mu_p = measurementfcn(X_kp1.x(:, i), x_chief, zeros(params.p));
    
%     X_kp1.w(i) = X_k.w(i);
%     for j = 1:params.p
%         pdf_ij = normpdf(y_kp1(j), mu_p(j), params.R(j, j));
%         X_kp1.w(i) = X_kp1.w(i) * pdf_ij;
%     end
    X_kp1.w(i) = X_k.w(i)*mvnpdf(y_kp1, mu_p, params.R);
    
end

w_total = sum(X_kp1.w);
X_kp1.w = X_kp1.w./w_total; % Normalize

end

