function RBPF = step_RBPF(RBPF, y_e1, y_e2, params, CONST)
particles = SIS_RBPF(RBPF(end).particles, y_e1, y_e2, params, CONST);
rbpfkp1.particles = particles;
% MMSE
rbpfkp1.x_mmse1 = RBPF_MMSE(particles, 'xp1');
rbpfkp1.x_mmse2 = RBPF_MMSE(particles, 'xp2');
rbpfkp1.xt_mmse = RBPF_MMSE(particles, 'x_t');
% MAP
rbpfkp1.x_map1 = RBPF_MAP(particles, 'xp1');
rbpfkp1.x_map2 = RBPF_MAP(particles, 'xp2');
rbpfkp1.xt_map = RBPF_MAP(particles, 'x_t');
% Measurement Residueals
rbpfkp1.y_res_mmse1 = y_e1 - measurementfcn(rbpfkp1.x_mmse1, rbpfkp1.xt_mmse, zeros(params.p));
rbpfkp1.y_res_mmse2 = y_e2 - measurementfcn(rbpfkp1.x_mmse2, rbpfkp1.xt_mmse, zeros(params.p));
rbpfkp1.y_res_map1 = y_e1 - measurementfcn(rbpfkp1.x_map1, rbpfkp1.xt_map, zeros(params.p));
rbpfkp1.y_res_map2 = y_e2 - measurementfcn(rbpfkp1.x_map2, rbpfkp1.xt_map, zeros(params.p));
% Other Stuff
rbpfkp1.P1 = RBPF_cov(particles, rbpfkp1.x_mmse1, 'xp1');
rbpfkp1.P2 = RBPF_cov(particles, rbpfkp1.x_mmse2, 'xp2');
rbpfkp1.Neff = 1/sum([particles(:).w].^2);
if all(~isnan(y_e1)) && all(~isnan(y_e2))
    new_particles = RBPF_resample(particles);
    rbpfkp1.particles = new_particles;
end
RBPF(end+1) = rbpfkp1;
end