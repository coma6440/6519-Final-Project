function out = RBPF_cov(particles)
x = zeros(length(particles(1).xp), length(particles));
for i = 1:length(particles)
    x(:,i) = particles(i).xp;
end
out = cov(x');
end