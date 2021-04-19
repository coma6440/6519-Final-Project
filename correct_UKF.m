function [xp, Pp] = correct_UKF(Pm, xm, t, y, x_chief, params, CONST)
n = params.n;
p = params.p;
lambda = params.lambda;
wm = params.wm;
wc = params.wc;
R = params.R;

S_b = chol(Pm);
ym = 0;
chi = zeros(length(xm),2*n+1);
for i = 1:(2*n+1)
    if i == 1
        chi(:,i) = xm;
    elseif i <= (n + 1)
        chi(:,i) = xm + sqrt(n + lambda)*S_b(i - 1,:)';
    else
        chi(:,i) = xm - sqrt(n + lambda)*S_b(i- (n+1),:)';
    end
    gam(:,i) = measurement(x_chief,chi(:,i));
    ym = ym + wm(i)*gam(:,i);
end
Si = zeros(length(y));
for i = 1:(2*n+1)
    Si = Si + wc(i)*(gam(:,i) - ym)*(gam(:,i) - ym)';
end
Si = Si + R;
C_xy = zeros(n,p);
for i = 1:(2*n+1)
    C_xy = C_xy + wc(i)*((chi(:,i) - xm)*(gam(:,i) - ym)');
end
K = C_xy*inv(Si);
xp = xm + K*(y - ym);
Pp = Pm - K*Si*K'; 
params.wm = wm;
params.wc = wc;
end