function RBPF = init_RBPF(CONST, params)
Ns = 100;
r0 = CONST.R_E + 2000; % [km]
pf0.x_mmse1 = NaN*zeros(params.n, 1);
pf0.x_mmse2 = NaN*zeros(params.n, 1);
pf0.y_res_mmse1 = NaN*zeros(params.p, 1);
pf0.y_res_mmse2 = NaN*zeros(params.p, 1);
pf0.x_map1 = NaN*zeros(params.n, 1);
pf0.x_map2 = NaN*zeros(params.n, 1);
pf0.y_res_map1 = NaN*zeros(params.p, 1);
pf0.y_res_map2 = NaN*zeros(params.p, 1);
pf0.P1 = 1*eye(params.n);
pf0.P1(2,2) = 0.0005;
pf0.P1(4,4) = 0.0005;
pf0.P1(6,6) = 0.0005;
pf0.P2 = pf0.P1;
pf0.Neff = Ns;
init_state_xt = [r0; 0; 0; 2.57; 0; 7.05];
init_state_xe1 = [r0+100; 1; 0; 8.75; 0; 0];
init_state_xe2 = [r0-200; 1; 0; 8.75; 0; 0];


Pp = 1*eye(6);
Pp(2,2) = 0.0005;
Pp(4,4) = 0.0005;
Pp(6,6) = 0.0005;
Pm = Pp;

filterUKF1 = trackingUKF(@transitionfcn, @measurementfcn, init_state_xe1, ...
         'ProcessNoise', params.UKF.Q, 'MeasurementNoise', params.UKF.R,...
         'StateCovariance', Pp,...
         'Alpha', 1e-3);
filterUKF2 = trackingUKF(@transitionfcn, @measurementfcn, init_state_xe2, ...
         'ProcessNoise', params.UKF.Q, 'MeasurementNoise', params.UKF.R,...
         'StateCovariance', Pp,...
         'Alpha', 1e-3);
for i = 1:Ns
    particles(i).x_t = init_state_xt;
    particles(i).w = 1/Ns;
    particles(i).xm1 = init_state_xe1;
    particles(i).xm2 = init_state_xe2;
    particles(i).xp1 = particles(i).xm1;
    particles(i).xp2 = particles(i).xm2;
    particles(i).Pp1 = Pp;
    particles(i).Pp2 = Pp;
    particles(i).Pm1 = Pm;
    particles(i).Pm2 = Pm;
    particles(i).UKF1 = clone(filterUKF1);
    particles(i).UKF2 = clone(filterUKF2);
end
pf0.xt_mmse = init_state_xt;
pf0.xt_map = init_state_xt;
pf0.particles = particles;
RBPF = pf0;
end