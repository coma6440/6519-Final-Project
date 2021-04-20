function predict_particles = particletransitionfcn(prev_particles,dt, CONST)
[~,c] = size(prev_particles);
%iterate over all particles
for i = 1:c
    predict_particles(:,i) = transitionfcn(prev_particles(:,i), dt, CONST);
end
end