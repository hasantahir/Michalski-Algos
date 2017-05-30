function [val, A, B] = LevinSidi (k, S, w, X, A, B)

% The type of transformation based on omega is decided in the calling
% routine

% First elements of the remainder term
% Align index terms due to X(1) = a, X starts from index 2 (k = 2)
B(k - 1) = 1/w;
A(k - 1) = S * B(k - 1);


for j = 1 : k - 2 % Change of limits to align Matlab indices with algorithm
    
    d = 1/X(k) - 1/X(k - j);
    % The indices of A and B need to be fixed
    A(k - j - 1) = (A( k - j) - A(k - j - 1))/d;
    B(k - j - 1) = (B( k - j) - B(k - j - 1))/d;
%     sprintf('i = %d, j = %d, A(%d) = %d',j,k,k - j -1 ,A(k - j -1))
end
val = A(1)/B(1);