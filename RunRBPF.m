function [RBPF] = RunRBPF(Y, params, CONST)
%RunFilters Executes the UKF and PF for each time in t_vec with
%corresponding measurement in Y

t = 0;

RBPF = init_RBPF(CONST, params);

for k = 1:(length(Y)-1)
    % Print Progress
    clc;
    fprintf('RBPF Progress:\t%3.0f %%\n', k/(length(Y)-1)*100);
    
    %Prediction step
    y_e1 = Y(:, k+1, 1);
    y_e2 = Y(:, k+1, 2);
    RBPF = step_RBPF(RBPF, y_e1, y_e2, params, CONST);
    
    
    %Step time forward
    t = t + params.dt;
end

end

