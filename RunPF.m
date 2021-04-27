function [PF] = RunPF(x_chief, Y, params, CONST, Ns)
%RunPF Runs the Particle Filter to estimate the deputy satellite states
%over time given the chief satellite's states and the relative measurements

r0 = CONST.R_E + 2000; % [km]
init_state = [r0+100; 1; 0; 8.75; 0; 0];

pf0.x_mmse = NaN*zeros(params.n, 1);
pf0.y_res_mmse = NaN*zeros(params.p, 1);
pf0.x_map = NaN*zeros(params.n, 1);
pf0.y_res_map = NaN*zeros(params.p, 1);
pf0.P = 1*eye(params.n);
pf0.P(2,2) = 0.0005;
pf0.P(4,4) = 0.0005;
pf0.P(6,6) = 0.0005;
pf0.Neff = Ns;
PF(1) = pf0;

% Initial Particle Set
X_k.x = init_state.*ones(params.n, Ns); % Particles
X_k.w = (1/Ns)*ones(Ns, 1); % Weights

% figure;

%Moving averages
w_fast = 1;
w_slow = 1;
alpha_slow = 0.1;
alpha_fast = 0.001;
nu = 2;

for k = 1:params.N-1
    clc;
    fprintf('PF Progress:\t%3.0f %%\n', k/(params.N-1)*100);
    
%     pf_k = PF(k);
    y_kp1 = Y(:, k+1);
    
    X_kp1 = SIS(X_k, y_kp1, x_chief(:, k+1), params, CONST);
    
    pf_kp1.x_mmse = CalcMMSE(X_kp1);
    pf_kp1.y_res_mmse = y_kp1 - measurementfcn(pf_kp1.x_mmse, x_chief(:, k+1), zeros(params.p));
    pf_kp1.x_map = CalcMAP(X_kp1);
    pf_kp1.y_res_map = y_kp1 - measurementfcn(pf_kp1.x_map, x_chief(:, k+1), zeros(params.p));
    pf_kp1.P = cov(X_kp1.x');
    pf_kp1.Neff = 1/(sum((X_kp1.w).^2));
    PF(k+1) = pf_kp1;
    
%     plot(X_kp1.w, 'o', 'LineWidth', 2)
    
    X_kp1 = Resample(X_kp1, w_fast, w_slow, nu);
    X_k = X_kp1;
    
    %Update adaptive injection
    w_fast = w_fast + alpha_fast*(mean(X_k.w) - w_fast);
    w_slow = w_slow + alpha_slow*(mean(X_k.w) - w_slow);
end

end

