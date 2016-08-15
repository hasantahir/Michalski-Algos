function [rem,R]  = MosigMichalski( mu, k, s, omega, X, R)

R(k) = s;
% for j = 0 : k  % Change of limits from algo to prevent array access errors
%     d = X(k - j +1) - X(k -j);
%     eta = omega/(1 + mu * (j-1) * d/(X(k -j)));
%     R(k -j) = R(k-j+1) - eta*R(k-j)/(1-eta);
% end

j = 2;
while ( j < k)
    d = X(abs(k - j) +1) - X(abs(k -j));
    eta = omega/(1 + mu * (j-1) * d/(X(abs(k -j))));
    R(abs(k -j)) = R(abs(k-j)+1) - eta*R(abs(k-j))/(1-eta);
    j = j +1;
end
rem = R(1);
end