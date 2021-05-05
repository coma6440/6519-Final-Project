function RBPF = step_RBPF(RBPF, y_e, params, CONST)
particles = SIS_RBPF(RBPF(end).particles, y_e, params, CONST);
rbpfkp1.particles = particles;
% MMSE
rbpfkp1.x_mmse = RBPF_MMSE(particles, 'xp');
rbpfkp1.xt_mmse = RBPF_MMSE(particles, 'x_t');
rbpfkp1.y_res_mmse = y_e - measurementfcn(rbpfkp1.x_mmse, rbpfkp1.xt_mmse, zeros(params.p));
% MAP
rbpfkp1.x_map = RBPF_MAP(particles, 'xp');
rbpfkp1.xt_map = RBPF_MAP(particles, 'x_t');
rbpfkp1.y_res_map = y_e - measurementfcn(rbpfkp1.x_map, rbpfkp1.xt_map, zeros(params.p));
% Other Stuff
rbpfkp1.P = RBPF_cov(particles, rbpfkp1.x_mmse);
rbpfkp1.Neff = 1/sum([particles(:).w].^2);
if all(~isnan(y_e))
    new_particles = RBPF_resample(particles);
    rbpfkp1.particles = new_particles;
end
RBPF(end+1) = rbpfkp1;
end