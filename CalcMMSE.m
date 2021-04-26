function [x_mmse] = CalcMMSE(X_k)
%CalcMMSE Calculates the MMSE estimate using the states and weights of X_k

[~, Ns] = size(X_k.x);
mmse_sum = 0;
for i = 1:Ns
    mmse_sum = mmse_sum + X_k.w(i)*X_k.x(:, i);
end
x_mmse = mmse_sum;

end

