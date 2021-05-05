function out = RBPF_cov(particles, mean)
n = length(particles(1).xp);
P = zeros(n,n);
parfor i = 1:length(particles)
    P = P + particles(i).w*(particles(i).xp - mean)*(particles(i).xp - mean)';
end
out = P;
end