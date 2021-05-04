function x_mmse = RBPF_MMSE(particles)
Ns = length(particles);
x_mmse = 0;
for i = 1:Ns
    x_mmse = x_mmse + particles(i).w*particles(i).xp;
end
end