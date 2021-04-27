function [X_kp1] = Resample(X_kp1, w_fast, w_slow, nu)
%Resample Resamples the particles X_kp1 from the distribution of particles

[n, Ns] = size(X_kp1.x);

x_resamp = zeros(n, Ns);
w_cdf = cumsum(X_kp1.w);

for i = 1:Ns
    pSamp = rand;
    ii = find(pSamp < w_cdf, 1);
    x_resamp(:, i) = X_kp1.x(:, ii);
end

%Adaptive Injection
n_inject = round(Ns*max([0,1-nu*(w_fast/w_slow)]));
if n_inject > 0
    error("Time to implement adaptive injection in resample.m")
end

X_kp1.x = x_resamp;
X_kp1.w(:) = 1/Ns;

end

