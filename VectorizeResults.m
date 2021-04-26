function [KF_vec] = VectorizeResults(KF)
%VectorizeResults() takes the results of a Kalman Filter array of structs
%that contains a results struct for each time step k and converts it into a
%single vector of states and covariance matrices
%
% KF must have the following fields for each struct: .x_p, .P_p

fnames = fieldnames(KF);
N = length(KF);

for i = 1:length(fnames)
    [n, p] = size(KF(1).(fnames{i}));
    if p > 1
        KF_vec.(fnames{i}) = zeros(n, p, N);
        
        for k = 1:N
            KF_vec.(fnames{i})(:, :, k) = KF(k).(fnames{i});
        end
    else
        KF_vec.(fnames{i}) = zeros(n, N);
        for k = 1:N
            KF_vec.(fnames{i})(:, k) = KF(k).(fnames{i});
        end
    end
end

end

