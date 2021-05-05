function pp1 = SIS_RBPF(particles, y_e, params, CONST)
%SIS Sequential Importance Sampling for the discrete time system

[~, Ns] = size(particles);
dt = params.dt;

for i = 1:Ns
    pp1(i).x_t = zeros(params.n, 1);
    pp1(i).w = 0;
end

parfor i = 1:Ns
    mu_q = transitionfcn(particles(i).x_t, dt, CONST);
    pp1(i).x_t = mvnrnd(mu_q, params.PF.Q)';

%         if ~isnan(y_t) % If a valid measurement is taken at time k+1
%             mu_p = measurementfcn(pp1.x_t, x_meas, zeros(params.p));
%             pp1.w(i) = particles.w(i)*mvnpdf(y_t, mu_p, params.PF.R);
%         else
%             pp1.w(i) = particles.w(i);
%         end
    [pp1(i).xm, pp1(i).Pm] = predict(particles(i).UKF, dt, CONST);
    if all(~isnan(y_e))
        [pp1(i).xp, pp1(i).Pp] = correct(particles(i).UKF, y_e, pp1(i).x_t, zeros(params.p, params.p));
        [yres, Sk] = residual(particles(i).UKF, y_e, {pp1(i).x_t, zeros(params.p, params.p)});
        ym_hat = y_e - yres;
        pp1(i).w = particles(i).w*mvnpdf(y_e, ym_hat, Sk);
    else
        pp1(i).xp = pp1(i).xm;
        pp1(i).Pp = pp1(i).Pm;
        pp1(i).w = particles(i).w;
    end
    
    pp1(i).UKF = clone(particles(i).UKF);
end

w_total = sum([pp1(:).w]);

if w_total ~= 0
    for i = 1:Ns
        pp1(i).w = pp1(i).w/w_total;
    end
else
    error('Sum of the particle weights equals 0. All particles have died.');
end
end
