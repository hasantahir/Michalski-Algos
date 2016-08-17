function [val, A, B] = LevinSidi (k,s, omega, X, A, B)

B(k-1) = 1/omega;
A(k-1) = s * B(k-1);

for j = 2 : k - 1
    d = 1/X(k) - 1/X(k - j + 1);
    A(k - j + 1) = (A( k - j + 2) - A(k - j + 1))/d;
    B(k - j + 1) = (B( k - j + 2) - B(k - j + 1))/d;
end
val = A(1)/B(1);