function [rem,R]  = MosigMichalski( mu, k, s, omega, X, R)

R(k-1) = s;
% for j = 0 : k  % Change of limits from algo to prevent array access errors
%     d = X(k - j +1) - X(k -j);
%     eta = omega/(1 + mu * (j-1) * d/(X(k -j)));
%     R(k -j) = R(k-j+1) - eta*R(k-j)/(1-eta);
% end

for j  = 2 : k - 1
    d = X(k - j + 2) - X(k - j + 1);
    eta = omega / (1 + mu * (j-1) * d / (X(k - j + 1)));
    R(k - j + 1) = R(k-j+2) - eta * R(k - j + 1) / (1-eta);
end
rem = R(1);
end