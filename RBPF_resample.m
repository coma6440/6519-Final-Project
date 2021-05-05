function new_particles = RBPF_resample(particles)
% Ns = length(particles);
% w_cdf = cumsum([particles(:).w]);
% new_particles = particles([]);
% for i = 1:Ns
%     pSamp = rand;
%     ii = find(pSamp <= w_cdf, 1);
%     
%     if ~isempty(ii)
%         new_particles(i) =  particles(ii);
%         new_particles(i).UKF = clone(particles(ii).UKF);
%     else
%         error("Error in resampling");
%     end
% end
% 
% for i = 1:Ns
%     new_particles(i).w = 1/Ns;
% end

% Cull dead particles
particles([particles(:).w] < 1e-5) = [];
w_total = sum([particles(:).w]);
for i = 1:length(particles)
    particles(i).w = particles(i).w/w_total;
end
new_particles = particles;
w_cdf = cumsum([particles(:).w]);
while length(new_particles) < 100
    pSamp = rand;
    ii = find(pSamp <= w_cdf, 1);
    new_particles(end + 1) = particles(ii);
    new_particles(end).UKF = clone(particles(ii).UKF);
end

for i = 1:length(new_particles)
    new_particles(i).w = 1/100;
end

end