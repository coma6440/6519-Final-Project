function out = RBPF_cov(particles, mean, field_name)
n = length(particles(1).(field_name));
P = zeros(n,n);
parfor i = 1:length(particles)
    P = P + particles(i).w*(particles(i).(field_name) - mean)*(particles(i).(field_name) - mean)';
end
out = P;
end