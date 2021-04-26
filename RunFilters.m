function [UKF, PF] = RunFilters(x_chief, Y, params, CONST)
%RunFilters Executes the UKF and PF for each time in t_vec with
%corresponding measurement in Y

r0 = CONST.R_E + 2000; % [km]
init_state = [r0+100; 1; 0; 8.75; 0; 0];
UKF.xm = zeros(params.n, params.N);
UKF.xp = UKF.xm;
UKF.xm(:,1) = init_state;
UKF.xp(:,1) = init_state;

PF.xm = zeros(params.n, params.N);
PF.xp = UKF.xm;
PF.xm(:,1) = init_state;
PF.xp(:,1) = init_state;

UKF.Pp = 1*eye(6);
UKF.Pp(2,2) = 0.0005;
UKF.Pp(4,4) = 0.0005;
UKF.Pp(6,6) = 0.0005;
PF.Pp = UKF.Pp;

UKF.Pm = UKF.Pp;

t = 0;

%UKF
filterUKF = trackingUKF(@transitionfcn, @measurementfcn, init_state, ...
         'ProcessNoise', params.Q, 'MeasurementNoise', params.R,...
         'StateCovariance', UKF.Pp,...
         'Alpha', 1e-3);
%Particle Filter  
filterPF = trackingPF(@particletransitionfcn, @particlemeasurementfcn, ...
           init_state, 'ProcessNoise', params.Q, 'MeasurementNoise', params.R,...
           'StateCovariance', PF.Pp, 'NumParticles', 500);

for k = 1:(length(Y)-1)
    % Print Progress
    clc;
    fprintf('Progress:\t%3.0f %%\n', k/(length(Y)-1)*100);
    
    %Prediction step
    [UKF.xm(:,k+1), UKF.Pm(:,:,k+1)] = predict(filterUKF, params.dt, CONST);
    [PF.xm(:,k+1), PF.Pm(:,:,k+1)] = predict(filterPF, params.dt, CONST);
    
    %Measurement residual
    [UKF.yres(:,k+1), UKF.Sk(:,:,k+1)] = residual(filterUKF, Y(:,k+1), {x_chief(:,k+1), zeros(params.p, params.p)});
    %[PF.yres(:,k+1), PF.Sk(:,:,k+1)] = residual(filterPF, Y(:,k+1), {x_chief(:,k+1), zeros(params.p, params.p)});
    
    %Correction step
    [UKF.xp(:,k+1), UKF.Pp(:,:,k+1)] = correct(filterUKF,Y(:,k+1), x_chief(:,k+1), zeros(params.p, params.p));
    [PF.xp(:,k+1), PF.Pp(:,:,k+1)] = correct(filterPF,Y(:,k+1), x_chief(:,k+1), zeros(params.p, params.p));
    
    %Step time forward
    t = t + params.dt;
end

end

