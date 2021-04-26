function y_hat = particlemeasurementfcn(particles, X_chief, R)
[~,c] = size(particles);
%iterate over all particles
parfor i = 1:c
    y_hat(:,i) = measurementfcn(particles(:,i), X_chief, R);
end
end
