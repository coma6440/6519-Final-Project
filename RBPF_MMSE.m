function x_mmse = RBPF_MMSE(particles, field_name)
Ns = length(particles);
x_mmse = 0;
for i = 1:Ns
    x_mmse = x_mmse + particles(i).w*particles(i).(field_name);
end
end