function [x_map] = CalcMAP(X_k)
%CalcMMSE Calculates the MAP estimate using the states and weights of X_k

[~, i_max] = max(X_k.w);
x_map = X_k.x(:, i_max);

end

