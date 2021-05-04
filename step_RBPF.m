function RBPF = step_RBPF(RBPF, x_meas, y_t, y_e, params, CONST)
particles = SIS_RBPF(RBPF(end).particles, y_t, y_e, x_meas, params, CONST);
rbpfkp1.particles = particles;
% MMSE
rbpfkp1.x_mmse = RBPF_MMSE(particles);
rbpfkp1.y_res_mmse = y_e - measurementfcn(rbpfkp1.x_mmse, x_meas, zeros(params.p));
% MAP
rbpfkp1.x_map = RBPF_MAP(particles);
rbpfkp1.y_res_map = y_e - measurementfcn(rbpfkp1.x_map, x_meas, zeros(params.p));
% Other Stuff
rbpfkp1.P = RBPF_cov(particles);
rbpfkp1.Neff = 1/sum([particles(:).w].^2);
RBPF(end+1) = rbpfkp1;
end