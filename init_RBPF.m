function RBPF = init_RBPF(CONST, params)
Ns = 500;
r0 = CONST.R_E + 2000; % [km]
pf0.x_mmse = NaN*zeros(params.n, 1);
pf0.y_res_mmse = NaN*zeros(params.p, 1);
pf0.x_map = NaN*zeros(params.n, 1);
pf0.y_res_map = NaN*zeros(params.p, 1);
pf0.P = 1*eye(params.n);
pf0.P(2,2) = 0.0005;
pf0.P(4,4) = 0.0005;
pf0.P(6,6) = 0.0005;
pf0.Neff = Ns;
init_state_xt = [r0+100; 1; 0; 8.75; 0; 0];

init_state_xe = [r0-200; 1; 0; 8.75; 0; 0];
UKF.xm = zeros(params.n, params.N);
UKF.xp = UKF.xm;
UKF.xm(:,1) = init_state_xt;
UKF.xp(:,1) = init_state_xt;

Pp = 1*eye(6);
Pp(2,2) = 0.0005;
Pp(4,4) = 0.0005;
Pp(6,6) = 0.0005;
Pm = Pp;

filterUKF = trackingUKF(@transitionfcn, @measurementfcn, init_state_xe, ...
         'ProcessNoise', params.UKF.Q, 'MeasurementNoise', params.UKF.R,...
         'StateCovariance', Pp,...
         'Alpha', 1e-3);
for i = 1:Ns
    particles(i).x_t = init_state_xt;
    particles(i).w = 1/Ns;
    particles(i).xm = init_state_xe;
    particles(i).xp = particles(i).xm;
    particles(i).Pp = Pp;
    particles(i).Pm = Pm;
    particles(i).UKF = clone(filterUKF);
end
pf0.particles = particles;
RBPF = pf0;
end