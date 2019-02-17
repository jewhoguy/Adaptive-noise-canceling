function [e,w,W] = NLMS(d, u, mu, M, a)
% Normalized Least Squares algorithm

% Inputs:
% d - desired signal
% u - input signal
% mu - µ - step-size
% M - filter taps
% a - NLMS constant

% Outputs:
% e - error vector
% w - last filter coeff
% W - saved filter coeff

% Init:
N_max = length(d);
u = [zeros(M-1,1);u];
w = zeros(M,1);
y = zeros(N_max,1);
e = zeros(N_max,1);

for n = 1:N_max
    uu = u(n+M-1:-1:n);
    y(n) = w'*uu;
    e(n) = d(n) - y(n);
    w = w + (mu/(a+uu'*uu))*e(n)*uu; %NLMS weight update rule
    W(:,n) = w;
end

end

