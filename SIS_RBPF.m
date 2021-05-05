function pp1 = SIS_RBPF(particles, y_e1, y_e2, params, CONST)
%SIS Sequential Importance Sampling for the discrete time system

[~, Ns] = size(particles);
dt = params.dt;

for i = 1:Ns
    pp1(i).x_t = zeros(params.n, 1);
    pp1(i).w = 0;
end

Q =  params.PF.Q;
p = params.p;
parfor i = 1:Ns
    mu_q = transitionfcn(particles(i).x_t, dt, CONST);
    pp1(i).x_t = mvnrnd(mu_q, Q)';

%         if ~isnan(y_t) % If a valid measurement is taken at time k+1
%             mu_p = measurementfcn(pp1.x_t, x_meas, zeros(params.p));
%             pp1.w(i) = particles.w(i)*mvnpdf(y_t, mu_p, params.PF.R);
%         else
%             pp1.w(i) = particles.w(i);
%         end
    [pp1(i).xm1, pp1(i).Pm1] = predict(particles(i).UKF1, dt, CONST);
    [pp1(i).xm2, pp1(i).Pm2] = predict(particles(i).UKF2, dt, CONST);
    if all(~isnan(y_e1)) && all(~isnan(y_e2))
        [pp1(i).xp1, pp1(i).Pp1] = correct(particles(i).UKF1, y_e1, pp1(i).x_t, zeros(p, p));
        [pp1(i).xp2, pp1(i).Pp2] = correct(particles(i).UKF2, y_e2, pp1(i).x_t, zeros(p, p));
                
        [yres1, Sk1] = residual(particles(i).UKF1, y_e1, {pp1(i).x_t, zeros(p, p)});
        [yres2, Sk2] = residual(particles(i).UKF2, y_e2, {pp1(i).x_t, zeros(p, p)});
        
        ym_hat1 = y_e1 - yres1;
        ym_hat2 = y_e2 - yres2;
        
%         test1 = mvnpdf(y_e1, ym_hat1, Sk1);
%         test2 = mvnpdf(y_e2, ym_hat2, Sk2);
%         
%         fprintf("pdf1 = %f, pdf2 = %f\n", test1, test2);
        pp1(i).w = particles(i).w*mvnpdf(y_e1, ym_hat1, Sk1)*mvnpdf(y_e2, ym_hat2, Sk2);
    else
        pp1(i).xp1 = pp1(i).xm1;
        pp1(i).xp2 = pp1(i).xm2;
        pp1(i).Pp1 = pp1(i).Pm1;
        pp1(i).Pp2 = pp1(i).Pm2;
        pp1(i).w = particles(i).w;
    end
    
    pp1(i).UKF1 = clone(particles(i).UKF1);
    pp1(i).UKF2 = clone(particles(i).UKF2);
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
