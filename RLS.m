function [e,w] = RLS(d,u,M,lambda,delta)
% Input:
% d - desired signal
% u - input signal (noise)
% M - number of taps
% lambda - forgetting factor
% delta - init param for P

% Output:
% e - error vector
% w - last filter coeff
% W - saved filter coeff

% Init:
w = zeros(M,1);
P = delta*eye(M);

N_max = length(d);
u = [zeros(M-1,1);u];
e = zeros(N_max,1);

for n = 1:N_max
    uu = u(n+M-1:-1:n);
    y = d(n);
    v = P*uu;
    k = v/(lambda+uu'*v);
    e(n) = y-w'*uu;
    
    w = w+k*e(n);
    P = P/lambda-k*v'/lambda;
end

% % Didn't work % %
% for n = (M+1):N_max
%     uu = u(n:-1:(n-M+1));
%     v = uu'*P;
%     gamma = lambda + v*uu;
%     k = v'/gamma;
%     e = d(n)-uu'*w;
%     w = w+k*e;
%     Pp = k*v;
%     P = (P-Pp)/lambda;
% end

end