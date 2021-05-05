function x_map = RBPF_MAP(particles, field_name)
max_w = 0;
for i = 1:length(particles)
    if particles(i).w > max_w
        max_w = particles(i).w;
        x_map = particles(i).(field_name);
    end
end
end